use pyo3::exceptions::socket;
use pyo3::prelude::*;
use pyo3::types::PyList;
use crate::filehandler::{read_bedfiles, chrombounds_from_bw, bwintervals, header_matrix, write_matrix};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;
use itertools::Itertools;
use crate::calc::{mean_float, median_float, max_float, min_float, sum_float};

#[pyfunction]
pub fn r_computematrix(
    mode: &str,
    bedlis: Py<PyList>,
    bwlis: Py<PyList>,
    sampleslabel: Py<PyList>,
    upstream: u32,
    downstream: u32,
    unscaled5prime: u32,
    unscaled3prime: u32,
    regionbodylength: u32,
    binsize: u32,
    missingdatazero: bool,
    nanafterend: bool, // end regions will treated as nans. Default is false.
    averagetypebins: &str, // operation to summarize values over bins. Default is mean.
    sortregions: &str, // either ascend, descend or keep. Default is keep (and ignores sortusing).
    sortusing: &str, // metric to sort on. Either mean median max min sum region_length. Default is mean.
    sortusingsamples: Py<PyList>, // list of samples to sort on. If empty, use all samples.
    referencepoint: &str,
    nproc: usize,
    verbose: bool,
    ofile: &str
) -> PyResult<()> {
    // Extract the bed and bigwig files from pyList to Vec.
    let mut bed_files: Vec<String> = Vec::new();
    let mut bw_files: Vec<String> = Vec::new();
    let mut samples_label: Vec<String> = Vec::new();
    let mut sort_using_samples: Vec<u32> = Vec::new();
    Python::with_gil(|py| {
        bed_files = bedlis.extract(py).expect("Failed to retrieve bed files.");
        bw_files = bwlis.extract(py).expect("Failed to retrieve bigwig filess.");
        samples_label = sampleslabel.extract(py).expect("Failed to retrieve samples label.");
        sort_using_samples = sortusingsamples.extract(py).expect("Failed to retrieve the samples to sort on.");
    });
    // Assert that samples_label equals bw_files, if samples_label is not empty.
    if !samples_label.is_empty() {
        assert_eq!(samples_label.len(), bw_files.len(), "Number of samplelabels do not match number of bigwig files.");
    }
    // Assert that sort_using_samples is smaller or equal to bw_files.
    assert!(sort_using_samples.len() <= bw_files.len(), "Number of samples to sort on is larger than number of bigwig files provided.");
    // Assert that no value in sort_using_samples is larger than bw_files.
    // assert!(sort_using_samples.max().unwrap() - 1 <= bw_files.len() as u32, "One or more samples to sort on are larger than the number of bigwig files provided.");

    // Get chromosome boundaries from first bigwig file.
    let chromsizes = chrombounds_from_bw(&bw_files.get(0).unwrap());
    // compute number of columns
    let bpsum = &upstream + &downstream + &unscaled5prime + &unscaled3prime + &regionbodylength;
    // Get the 'basepaths' of the bed files to use as labels later on.
    let mut bedlabels: Vec<String> = Vec::new();
    for bed in bed_files.iter() {
        let entryname = Path::new(bed)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        bedlabels.push(entryname);
    }
    if samples_label.is_empty() {
        // no samples labels provided via CLI, retrieve them from bigwig names.
        for bw in bw_files.iter() {
            let entryname = Path::new(bw)
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            samples_label.push(entryname);
        }
    }
    // Define the scaling regions in a struct
    let scale_regions = Scalingregions {
        upstream: upstream,
        downstream: downstream,
        unscaled5prime: unscaled5prime,
        unscaled3prime: unscaled3prime,
        regionbodylength: regionbodylength,
        binsize: binsize,
        cols_expected: ((bw_files.len() * bpsum as usize) / binsize as usize),
        missingdata_as_zero: missingdatazero,
        nan_after_end: nanafterend,
        referencepoint: referencepoint.to_string(),
        mode: mode.to_string(),
        bwfiles: bw_files.len(),
        avgtype: averagetypebins.to_string(),
        verbose: verbose,
        proc_number: nproc,
        bedlabels: bedlabels,
        bwlabels: samples_label
    };
    // Parse regions from bed files. Note that we retain the name of the bed file (in case there are more then 1)
    // Additionaly, score and strand are also retained, if it's a 3-column bed file we just fill in '.'
    let (mut regions, regionsizes) = read_bedfiles(&bed_files);
    let slopregions = slop_regions(
        &regions,
        &scale_regions,
        &chromsizes
    );
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let mut matrix: Vec<Vec<f32>> = pool.install(|| {
        bw_files.par_iter()
            .map(|i| bwintervals(&i, &regions, &slopregions, &scale_regions))
            .reduce(
                || vec![vec![]; regions.len()],
                |mut acc, vec_of_vecs| {
                    for (i, inner_vec) in vec_of_vecs.into_iter().enumerate() {
                        acc[i].extend(inner_vec);
                    }
                    acc
                },
            )
    });
    // Resort the matrix, if this is requested.
    if sortregions != "keep" {
        println!("Sorting output matrix with settings: sortRegions: {}, sortUsing {}", sortregions, sortusing);
        // If sortusingsamples is set, we need a vector to subset the columns of interest
        let mut cols_of_interest: Vec<usize> = Vec::new();
        if !sort_using_samples.is_empty() {
            for sample_ix in sort_using_samples.iter() {
                // Note that sort_using_samples is assumed to be 1-index. Hence we need to subtract 1.
                let start = (sample_ix - 1) * bpsum;
                let end = start + bpsum;
                cols_of_interest.extend(start as usize..end as usize);
            }
        }
        let mut regionslices: Vec<(usize, usize)> = Vec::new();
        let mut rstart: usize = 0;
        let mut rend: usize = 0;
        let mut lastregion = &regions[0].5;
        for (ix, region) in regions.iter().enumerate() {
            if &region.5 != lastregion {
                regionslices.push((rstart, rend));
                rstart = ix;
                rend = ix;
                lastregion = &region.5;
            }
            rend = ix;
        }
        regionslices.push((rstart, rend));
        let mut sortedix: Vec<usize> = Vec::new();
        if sortusing == "region_length" {
            if !sort_using_samples.is_empty() {
                println!("Sort using samples is set ({:?}), but is not used when sorting on region_length. It is thus ignored.", sort_using_samples);
            }
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
                    println!("region length sorting. Start = {}, end = {}", start, end);
                    let rslice = &regions[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, vals)| {
                            let regionlength = vals.2 - vals.1;
                            (ix + *start, regionlength)
                        })
                        .collect::<Vec<_>>()
                        .iter()
                        .sorted_by(|ix, metric| ix.1.partial_cmp(&metric.1).unwrap())
                        .map(|(ix, _)| *ix)
                        .collect::<Vec<usize>>();
                    match sortregions {
                        "ascend" => tix,
                        "descend" => tix.into_iter().rev().collect(),
                        _ => panic!("If sortRegions is not keep, it should be either ascend or descend. Not {}", sortregions),
                    }
                })
                .collect();
            println!("Sortedix = {:?}", sortedix);
        } else {
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
                    println!("Matrix sorting. Start = {}, end = {}", start, end);
                    let rslice = &matrix[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, vals)| {
                            let subset: Vec<_> = if cols_of_interest.is_empty() {
                                vals.iter().collect()
                            } else {
                                cols_of_interest
                                    .iter()
                                    .filter_map(|&index| vals.get(index))  // `vec.get(index)` returns Option<&T>
                                    .collect()
                            };
                            let metric = match sortusing {
                                "mean" => mean_float(subset),
                                "median" => median_float(subset),
                                "max" => max_float(subset),
                                "min" => min_float(subset),
                                "sum" => sum_float(subset),
                                _ => panic!("Sortusing should be either mean, median, max, min, sum or region_length. Not {}", sortusing),
                            };
                            (ix + *start, metric)
                        })
                        .collect::<Vec<_>>()
                        .iter()
                        .sorted_by(|ix, metric| ix.1.partial_cmp(&metric.1).unwrap())
                        .map(|(ix, _)| *ix)
                        .collect::<Vec<usize>>();
                    match sortregions {
                        "ascend" => tix,
                        "descend" => tix.into_iter().rev().collect(),
                        _ => panic!("If sortRegions is not keep, it should be either ascend or descend. Not {}", sortregions),
                    }
                })
                .collect();
            }
        // assert sorted ix length == matrix length == regions length
        assert_eq!(
            sortedix.len(),
            matrix.len(),
            "Length of sorted indices does not match matrix length: {} != {}", sortedix.len(), matrix.len()
        );
        assert_eq!(
            sortedix.len(),
            regions.len(),
            "Length of sorted indices does not match regions length: {} ! = {}", sortedix.len(), regions.len()
        );
        let mut sortedmatrix: Vec<Vec<f32>> = Vec::new();
        let mut sortedregions: Vec<(String, u32, u32, String, String, String)> = Vec::new();
        // Reorder matrix & regions
        sortedmatrix = sortedix
            .iter()
            .map(|ix| matrix[*ix].clone())
            .collect();
        sortedregions = sortedix
            .into_iter()
            .map(|ix| regions[ix].clone())
            .collect();
        write_matrix(
            header_matrix(&scale_regions, regionsizes),
            sortedmatrix,
            ofile,
            sortedregions
        );
    } else {
        write_matrix(
            header_matrix(&scale_regions, regionsizes),
            matrix,
            ofile,
            regions
        );
    }
    
    Ok(())
}

fn slop_regions(
    regions: &Vec<(String, u32, u32, String, String, String)>,
    scale_regions: &Scalingregions,
    chromsizes: &HashMap<String, u32>
) -> Vec<Vec<(u32, u32)>> {
    let mut regionranges: Vec<Vec<(u32, u32)>> = Vec::new();
    // Idea is to create a vector of tuples with start and end of every bin (binsize passed by computeMatrix).
    // The number of columns per region needs to be fixed per region.
    // Note that the before / after could mean that we run out of chromosome. 
    // Invalid regions (later to be encoded as NA or 0), will be pushed as (0,0) tuples.
    // Note that if nan_after_end is set to true, we will push (0,0) tuples after the end of the region.
    let col_expected = scale_regions.cols_expected / scale_regions.bwfiles;
    for region in regions.iter() {
        // To implement: 
        // // scale-regions + unscaled 5 and 3
        // // + and - encodings
        let chromend: u32 = *chromsizes.get(&region.0).unwrap();
        assert!(region.2 <= chromend, "Region end goes beyond chromosome boundary. Fix your bed files. {:?} > {}", region, chromend);
        assert!(region.1 <= chromend, "Region start goes beyond chromosome boundary. Fix your bed files. {:?} > {}", region, chromend);

        let anchorpoint = match scale_regions.referencepoint.as_str() {
            "TSS" => region.1,
            "TES" => region.2,
            "center" => (region.1 + region.2) / 2,
            _ => panic!("Reference should either be TSS, TES or center. {:?} is not supported.", scale_regions.referencepoint),
        };

        let mut regionsizes: Vec<(u32, u32)> = Vec::new();
        let mut absstart: i32 = anchorpoint as i32 - scale_regions.upstream as i32;
        let absstop: i32 = anchorpoint as i32 + scale_regions.downstream as i32;
        while absstart < absstop {
            let bin = absstart + scale_regions.binsize as i32;
            if absstart < 0  || bin > chromend as i32 {
                regionsizes.push((0,0));
            } else {
                // If we reached end of region, and nan_after_end is true, push (0,0)
                if scale_regions.nan_after_end && absstart as u32 >= region.2 {
                    regionsizes.push((0,0))
                } else {
                    regionsizes.push((absstart as u32, bin as u32))
                }
            }
            absstart = bin;
        }
        assert!(
            regionsizes.len() == col_expected,
            "Number of bins does not match expected number of columns: (CHROMLEN = {}) {:?} \n \n {:?}, \n \n {} != {}",
            *chromsizes.get(&region.0).unwrap(),
            region,
            regionsizes,
            regionsizes.len(),
            col_expected,
        );
        regionranges.push(regionsizes);
    }
    return regionranges;
}

pub struct Scalingregions {
    pub upstream: u32,
    pub downstream: u32,
    pub unscaled5prime: u32,
    pub unscaled3prime: u32,
    pub regionbodylength: u32,
    pub binsize: u32,
    pub cols_expected: usize,
    pub missingdata_as_zero: bool,
    pub nan_after_end: bool,
    pub referencepoint: String,
    pub mode: String,
    pub bwfiles: usize,
    pub avgtype: String,
    pub verbose: bool,
    pub proc_number: usize,
    pub bedlabels: Vec<String>,
    pub bwlabels: Vec<String>
}