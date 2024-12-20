use pyo3::prelude::*;
use pyo3::types::PyList;
use crate::filehandler::{read_bedfile, read_gtffile, chrombounds_from_bw, bwintervals, header_matrix, write_matrix};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;
use std::fmt;
use itertools::Itertools;
use crate::calc::{mean_float, median_float, max_float, min_float, sum_float};

#[pyfunction]
pub fn r_computematrix(
    mode: &str,
    regionlis: Py<PyList>,
    bwlis: Py<PyList>,
    sampleslabel: Py<PyList>,
    upstream: u32,
    downstream: u32,
    unscaled5prime: u32,
    unscaled3prime: u32,
    regionbodylength: u32,
    binsize: u32,
    missingdatazero: bool,
    metagene: bool,
    txnid: &str,
    exonid: &str,
    txniddesignator: &str,
    scale: f32, // scaling factor for writing out values. default is 1.0 (no scaling)
    nanafterend: bool, // end regions will treated as nans. Default is false.
    skipzeros: bool, // skip regions with all zeros. Default is false.
    minthresh: f32, // minimum threshold to keep a region. If not set it will equal 0.0
    maxthresh: f32, // maximum threshold to keep a region. if not set it will equal 0.0
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
    let mut region_files: Vec<String> = Vec::new();
    let mut bw_files: Vec<String> = Vec::new();
    let mut samples_label: Vec<String> = Vec::new();
    let mut sort_using_samples: Vec<u32> = Vec::new();
    Python::with_gil(|py| {
        region_files = regionlis.extract(py).expect("Failed to retrieve bed files.");
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
    // Get chromosome boundaries from first bigwig file.
    let chromsizes = chrombounds_from_bw(&bw_files.get(0).unwrap());
    // compute number of columns
    let bpsum = &upstream + &downstream + &unscaled5prime + &unscaled3prime + &regionbodylength;
    // Get the 'basepaths' of the bed files to use as labels later on.
    let mut regionlabels: Vec<String> = Vec::new();
    for bed in region_files.iter() {
        let entryname = Path::new(bed)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        regionlabels.push(entryname);
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
        scale: scale,
        nan_after_end: nanafterend,
        skipzero: skipzeros,
        minthresh: minthresh,
        maxthresh: maxthresh,
        referencepoint: referencepoint.to_string(),
        mode: mode.to_string(),
        bwfiles: bw_files.len(),
        avgtype: averagetypebins.to_string(),
        verbose: verbose,
        proc_number: nproc,
        regionlabels: regionlabels,
        bwlabels: samples_label
    };
    if verbose {
        println!("Region files: {:?}", &region_files);
        println!("Bigwig files: {:?}", &bw_files);
        println!("Samples labels: {:?}", scale_regions.bwlabels);
        println!("Sort using samples: {:?}", &sort_using_samples);
    }
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    
    // Parse regions from bed files. Note that we retain the name of the bed file (in case there are more then 1)
    // Additionaly, score and strand are also retained, if it's a 3-column bed file we just fill in '.'
    let mut regions: Vec<Region> = Vec::new();
    let mut regionsizes: HashMap<String, u32> = HashMap::new();

    let gtfparse = Gtfparse {
        metagene: metagene,
        txnid: txnid.to_string(),
        exonid: exonid.to_string(),
        txniddesignator: txniddesignator.to_string(),
    };

    region_files.iter()
        .map(|r| {
            let ext = Path::new(r)
                .extension()
                .and_then(|e| e.to_str())
                .map(|e| e.to_ascii_lowercase());

            match ext {
                Some(v) if v == "gtf".to_string() => read_gtffile(r, &gtfparse, chromsizes.keys().collect()),
                Some(v) if v == "bed".to_string() => read_bedfile(r, metagene, chromsizes.keys().collect()),
                _ => panic!("Only .bed and .gtf files are allowed as regions. File = {}, Extension = {:?}", r, ext),
            }
        })
        .for_each(|(reg, regsize)| {
            regions.extend(reg);
            regionsizes.insert(regsize.0, regsize.1);
        });
    
    // let (mut regions, regionsizes) = read_bedfiles(&bed_files, metagene);
    // Slop the regions
    
    let slopregions = pool.install(|| {
        regions.par_iter()
            .map(|region| slop_region(&region, &scale_regions, &chromsizes))
            .collect::<Vec<_>>()
    });
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
        if verbose {
            println!("Sorting output matrix with settings: sortRegions: {}, sortUsing {}", sortregions, sortusing);
        }
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
        let mut lastregion = &regions[0].name;
        for (ix, region) in regions.iter().enumerate() {
            if region.name != *lastregion {
                regionslices.push((rstart, rend));
                rstart = ix;
                rend = ix;
                lastregion = &region.name;
            }
            rend = ix;
        }
        regionslices.push((rstart, rend));
        let mut sortedix: Vec<usize> = Vec::new();
        if sortusing == "region_length" {
            if !sort_using_samples.is_empty() && verbose {
                println!("Sort using samples is set ({:?}), but is not used when sorting on region_length. It is thus ignored.", sort_using_samples);
            }
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &regions[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, region)| {
                            (ix + *start, region.regionlength)
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
        } else {
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
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
        let mut sortedregions: Vec<Region> = Vec::new();
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
            sortedregions,
            &scale_regions
        );
    } else {
        write_matrix(
            header_matrix(&scale_regions, regionsizes),
            matrix,
            ofile,
            regions,
            &scale_regions
        );
    }
    
    Ok(())
}

fn slop_region(
    region: &Region,
    scale_regions: &Scalingregions,
    chromsizes: &HashMap<String, u32>
) -> Vec<Bin> {
    // Idea is to create a vector of tuples with start and end of every bin (binsize passed by computeMatrix).
    // The number of columns per region needs to be fixed per region.
    // Note that the before / after could mean that we run out of chromosome. 
    // Invalid regions (later to be encoded as NA or 0), will be pushed as (0,0) tuples.
    // Note that if nan_after_end is set to true, we will push (0,0) tuples after the end of the region.

    // chromosome end
    let chromend: u32 = *chromsizes.get(&region.chrom).unwrap();

    // Assert region stays within chromosome boundary.
    region.assert_end(chromend);
    // Get the anchorpoint per region.
    let anchorpoint = region.get_anchorpoint(&scale_regions.referencepoint);
    // Create Vector of bins wrt. the anchorpoint.
    region.get_anchor_bins(anchorpoint, scale_regions, chromend)
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
    pub scale: f32,
    pub nan_after_end: bool,
    pub skipzero: bool,
    pub minthresh: f32,
    pub maxthresh: f32,
    pub referencepoint: String,
    pub mode: String,
    pub bwfiles: usize,
    pub avgtype: String,
    pub verbose: bool,
    pub proc_number: usize,
    pub regionlabels: Vec<String>,
    pub bwlabels: Vec<String>
}

#[derive(Clone)]
pub struct Gtfparse {
    pub metagene: bool,
    pub txnid: String,
    pub exonid: String,
    pub txniddesignator: String,
}

#[derive(Clone)]
pub enum Revalue {
    U(u32),
    V(Vec<u32>),
}

impl Revalue {
    pub fn rewrites(&self) -> String {
        match self {
            Revalue::U(v) => format!("{}", v),
            Revalue::V(vs) => vs.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

impl fmt::Debug for Revalue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Revalue::U(value) => write!(f, "U({})", value),
            Revalue::V(values) => write!(f, "V({:?})", values),
        }
    }
}

impl fmt::Display for Revalue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Revalue::U(value) => write!(f, "U({})", value),
            Revalue::V(values) => write!(f, "V({})", values.iter().map(|v| v.to_string()).collect::<Vec<String>>().join(", ")),
        }
    }
}

#[derive(Debug)]
pub enum Bin {
    Conbin(u32, u32),
    Catbin(Vec<(u32, u32)>),
}

#[derive(Clone, Debug)]
pub struct Region {
    pub chrom: String,
    pub start: Revalue,
    pub end: Revalue,
    pub score: String,
    pub strand: String,
    pub name: String,
    pub regionlength: u32
}

impl Region {
    pub fn assert_end(&self, chromend: u32) {
        match &self.end {
            Revalue::U(end) => {
                assert!(
                    *end <= chromend,
                    "Region end goes beyond chromosome boundary. Fix {}. {} {} {} (chr end = {})", self.name, self.chrom, self.start, self.end, chromend
                );
            },
            Revalue::V(ends) => {
                for end in ends.iter() {
                    assert!(
                        *end <= chromend,
                        "Region end goes beyond chromosome boundary. Fix {}. {} {} {} (chr end = {})", self.name, self.chrom, self.start, end, chromend
                    );
                }
            }
        }
    }

    pub fn get_anchorpoint(&self, referencepoint: &str) -> u32 {
        // reference-point mode.
        // What is exactly returned depends on a couple of parameters
        // what is the referencepoint : TSS, center, TES
        // what is the strand: +, -, . (note . is assumed to be +)
        // depending on if we have exon blocks (start / end are Revalue V == Vectors) or not (start / end are Revalue U == u32's)
        match referencepoint {
            "TSS" => {
                match self.strand.as_str() {
                    "+" | "." => match &self.start {Revalue::U(start) => *start, Revalue::V(starts) => starts[0]},
                    "-" => match &self.end {Revalue::U(end) => *end, Revalue::V(ends) => *ends.last().unwrap()},
                    _ => panic!("Strand should either be + or - or . {:?} is not supported.", self.strand),
                }
            },
            "TES" => {
                match self.strand.as_str() {
                    "+" | "." => match &self.end {Revalue::U(end) => *end, Revalue::V(ends) => *ends.last().unwrap()},
                    "-" => match &self.start {Revalue::U(start) => *start, Revalue::V(starts) => starts[0]},
                    _ => panic!("Strand should either be + or - or . {:?} is not supported.", self.strand),
                }
            },
            "center" => {
                // Here + or - doesn't matter. It is important though if we have 'metagenes' or not.
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        (*start + *end) / 2
                    },
                    (Revalue::V(starts), Revalue::V(ends)) => {
                        let exonlength: u32 = starts.iter().zip(ends.iter()).map(|(s, e)| e - s).sum();
                        let middle = exonlength / 2;
                        let mut cumsum: u32 = 0;
                        for (s, e) in starts.iter().zip(ends.iter()) {
                            cumsum += e - s;
                            if cumsum >= middle {
                                return s + (middle - (cumsum - (e - s)));
                            }
                        }
                    panic!(
                        "Middle of region not found. Fix {}. {}:{}-{}",
                        self.name, self.chrom, self.start, self.end
                    )
                    },
                    _ => panic!(
                        "Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                        self.name, self.chrom, self.start, self.end
                    ),
                }
            },
            _ => panic!(
                "Reference should either be TSS, TES or center. {:?} is not supported.",
                referencepoint
            ),
        }
    }

    pub fn get_anchor_bins(&self, anchorpoint: u32, scale_regions: &Scalingregions, chromend: u32) -> Vec<Bin> {
        // Given an anchorpoint, return a vector, start, end , middle
        // The order of the vector is always 5' -> 3', meaning 'increasing' for +/. regions, and 'decreasing' for - regions.
        // At this stage, two situations are possible:
        // - self.start / self.end are Revalue::U, meaning we are in 'non metagene' mode.
        // - self.start / self.end are Revalue::V, meaning we are in 'metagene' mode, and the bins returned are exon-aware.
        // We need a notion of bins that don't make sense (i.e. beyond chromosome boundaries). These are encoded as (0,0)
        let mut bins: Vec<Bin> = Vec::new();

        match self.strand.as_str() {
            "+" | "." => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // simplest scenario,
                        let mut absstart: i32 = anchorpoint as i32 - scale_regions.upstream as i32;
                        let absstop: i32 = anchorpoint as i32 + scale_regions.downstream as i32;
                        while absstart < absstop {
                            let bin = absstart + scale_regions.binsize as i32;
                            if absstart < 0 || absstart as u32 >= chromend || bin as u32 >= chromend {
                                bins.push(Bin::Conbin(0,0));
                            } else {
                                // If we reached end of region, and nan_after_end is true, push (0,0)
                                if scale_regions.nan_after_end && absstart as u32 >= *end {
                                    bins.push(Bin::Conbin(0,0));
                                } else {
                                    bins.push(Bin::Conbin(absstart as u32, bin as u32));
                                }
                            }
                            absstart = bin;
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorpoint;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.downstream {
                            if lastanchor >= chromend {
                                rightbins.push(Bin::Conbin(0,0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true
                                );
                                rightbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Left side.
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorpoint;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.upstream {
                            if lastanchor == 0 {
                                leftbins.push(Bin::Conbin(0,0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    false,
                                    false
                                );
                                leftbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Now we need to reverse the leftbins, as we walked backwards.
                        leftbins.reverse();
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    _ => panic!("Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                     self.name, self.chrom, self.start, self.end),
                }
            },
            "-" => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // Still simple scenario, we just gotta walk the other way.
                        let mut absstart: i32 = anchorpoint as i32 + scale_regions.upstream as i32;
                        let absstop: i32 = anchorpoint as i32 - scale_regions.downstream as i32;
                        while absstart > absstop {
                            let bin = absstart - scale_regions.binsize as i32;
                            if absstart as u32 > chromend || bin < 0 {
                                bins.push(Bin::Conbin(0,0));
                            } else {
                                if scale_regions.nan_after_end && absstart as u32 <= *end {
                                    bins.push(Bin::Conbin(0,0));
                                } else {
                                    // Push in the opposite direction so we walk backwards, but the interval is positive.
                                    bins.push(Bin::Conbin(bin as u32, absstart as u32));
                                }
                            }
                            absstart = bin;
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorpoint;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.upstream {
                            if lastanchor >= chromend {
                                rightbins.push(Bin::Conbin(0,0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    false,
                                    true
                                );
                                rightbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Left side.
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorpoint;
                        let mut walked_bps: u32 = 0;
                        while walked_bps < scale_regions.downstream {
                            if lastanchor == 0 {
                                leftbins.push(Bin::Conbin(0,0));
                                walked_bps += scale_regions.binsize;
                            } else {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false
                                );
                                leftbins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        // Note that now we need to go the exact opposite way as for the + strand as the 'highest position' is the 'starting point'.
                        rightbins.reverse();
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    _ => panic!("Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}. {}:{}-{}",
                     self.name, self.chrom, self.start, self.end),
                }
            },
            _ => panic!("Strand should either be + or - or . {:?} is not supported.", self.strand),
        };
        assert_eq!(
            bins.len(),
            scale_regions.cols_expected / scale_regions.bwfiles,
            "Number of bins does not match expected number of columns: {} != {}",
            bins.len(),
            scale_regions.cols_expected / scale_regions.bwfiles
        );
        bins
    }
}

fn refpoint_exonwalker(exons: &Vec<(u32, u32)>, anchor: u32, binsize: u32, chromend: u32, nan_after_end: bool, forward: bool) -> (Bin, u32) {
    // Basic function that walks over exons, and returns a Bin (Either Conbin or Catbin) and the last anchorpoint.
    let mut anchorix: Option<usize> = None;

    for (ix, exon) in exons.iter().enumerate() {
        if anchor >= exon.0 && anchor <= exon.1 {
            anchorix = Some(ix);
        }
    }
    if forward {
        // Walk forward (downstream, towards chromosome end)
        match anchorix {
            Some(i) => {
                // anchor sits in an exon. Check if anchor + binsize is also in same exon.
                if anchor + binsize <= exons[i].1 {
                    (Bin::Conbin(anchor, anchor + binsize), anchor + binsize)
                } else {
                    // anchor + binsize is not in same exon. We need a Catbin.
                    // Things are a bit more difficult here as well, as we need to walk exons.
                    let mut start_end_vec: Vec<(u32, u32)> = Vec::new();
                    start_end_vec.push( (anchor, exons[i].1) );
                    
                    let mut remainingbin: u32 = binsize - (exons[i].1 - anchor);
                    let mut lastix: usize = i;
                    let mut lastanchor: u32 = exons[i].1;
    
                    while remainingbin != 0 {
                        if lastix + 1 < exons.len() {
                            // next exon is available.
                            // Two options here:
                            // the remainder fits in the lastix + 1 exon. We are done.
                            // the remainder doesn't fit in the lastix + 1 exon. We need to walk further.
                            if exons[lastix+1].1 - exons[lastix+1].0 >= remainingbin {
                                // remainder fits in next exon.
                                start_end_vec.push( (exons[lastix+1].0, exons[lastix+1].0 + remainingbin) );
                                lastanchor = exons[lastix+1].0 + remainingbin;
                                remainingbin = 0;
                            } else {
                                // Remainder is larger then our exon. We need another walk.
                                start_end_vec.push( (exons[lastix+1].0, exons[lastix+1].1) );
                                remainingbin -= exons[lastix+1].1 - exons[lastix+1].0;
                                lastix += 1;
                                lastanchor = exons[lastix].1;
                            }
                        } else {
                            // No next exon available. Remainder can just be genomic.
                            // The last entry here can be changed to include the last part.
                            if nan_after_end {
                                start_end_vec.push((0,0));
                            } else {
                                let last = start_end_vec.last_mut().unwrap();
                                assert_eq!(
                                    last.1,
                                    lastanchor,
                                    "In the exon - genomic walk, our coordinates are not contiguous"
                                );
                                // Check we don't fall of the chromosome.
                                if lastanchor + remainingbin > chromend {
                                    last.1 = chromend;
                                    lastanchor = chromend;
                                } else {
                                    last.1 = lastanchor + remainingbin;
                                    lastanchor = lastanchor + remainingbin;
                                }
                            }
                            remainingbin = 0;
                        }
                    }
                    // We now have a Vec of start - end, we can construct a CatBin.
                    // Note that CatBins are (absstart, absstop, ((intstart1, intstart2), ...))
                    // This seems weird, but makes sure we need to slice the bigwig file only once per bin.
                    if start_end_vec.len() == 1 {
                        (Bin::Conbin(anchor, lastanchor), lastanchor)
                    } else {
                        (Bin::Catbin(start_end_vec), lastanchor)
                    }
                }
            },
            None => {
                // our anchor doesn't sit in exons. We just return the anchor + binsize as Bin
                if anchor + binsize > chromend {
                    (Bin::Conbin(anchor, chromend), chromend)
                } else {
                    (Bin::Conbin(anchor, anchor + binsize), anchor + binsize)
                }
            }
        }
    } else {
        // Walk backwards (upstream, towards chromosome start)
        match anchorix {
            Some(i) => {
                // anchor sits in an exon. Check if anchor - binsize is also in same exon.
                if anchor - binsize >= exons[i].0 {
                    (Bin::Conbin(anchor - binsize, anchor), anchor - binsize)
                } else {
                    // anchor + binsize is not in same exon. We need a Catbin.
                    // Things are a bit more difficult here as well, as we need to walk exons.
                    let mut start_end_vec: Vec<(u32, u32)> = Vec::new();
                    start_end_vec.push( (exons[i].0, anchor) );
                    
                    let mut remainingbin: u32 = binsize - (anchor - exons[i].0);
                    let mut lastix: usize = i;
                    let mut lastanchor: u32 = exons[i].0;

                    while remainingbin != 0 {
                        if lastix >= 1 {
                            // previous exon is available.
                            // Two options here:
                            // the remainder fits in the previous exon. We are done.
                            // the remainder doesn't fit in the previous exon. We need to walk further.
                            if exons[lastix-1].1 - exons[lastix-1].0 >= remainingbin {
                                // remainder fits in next exon.
                                start_end_vec.push( (exons[lastix-1].1 - remainingbin, exons[lastix-1].1) );
                                lastanchor = exons[lastix-1].1 - remainingbin;
                                remainingbin = 0;
                            } else {
                                // Remainder is larger then our exon. We need another walk.
                                start_end_vec.push( (exons[lastix-1].0, exons[lastix-1].1 ) );
                                remainingbin -= exons[lastix-1].1 - exons[lastix-1].0;
                                lastix -= 1;
                                lastanchor = exons[lastix].0;
                            }
                        } else {
                            // No previous exon available. Remainder can just be genomic.
                            // The last entry here can be changed to include the last part.
                            if nan_after_end {
                                start_end_vec.push( (0,0) );
                            } else {
                                let last = start_end_vec.last_mut().unwrap();
                                assert_eq!(
                                    last.0,
                                    lastanchor,
                                    "In the exon - genomic walk (reverse), our coordinates are not contiguous"
                                );
                                // Check we don't go in the negative.
                                if lastanchor < remainingbin {
                                    last.0 = 0;
                                    lastanchor = 0;
                                } else {
                                    last.0 = lastanchor - remainingbin;
                                    lastanchor = lastanchor - remainingbin;
                                }
                            }
                            remainingbin = 0;
                        }
                    }
                    // We now have a Vec of start - end, we can construct a CatBin.
                    // Note that CatBins are (absstart, absstop, ((intstart1, intstart2), ...))
                    // This seems weird, but makes sure we need to slice the bigwig file only once per bin.
                    if start_end_vec.len() == 1 {
                        (Bin::Conbin(lastanchor, anchor), lastanchor)
                    } else {
                        start_end_vec.reverse();
                        (Bin::Catbin(start_end_vec), lastanchor)
                    }
                }
            },
            None => {
                // our anchor doesn't sit in exons. We just return the anchor - binsize as Bin
                if anchor < binsize {
                    (Bin::Conbin(0, anchor), 0)
                } else {
                    (Bin::Conbin(anchor - binsize, anchor), anchor - binsize)
                }
            }
        }
    }
}