use pyo3::prelude::*;
use pyo3::types::PyList;
use crate::filehandler::{read_bedfiles, chrombounds_from_bw, bwintervals, header_matrix, write_matrix};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;

#[pyfunction]
pub fn r_computematrix(
    mode: &str,
    bedlis: Py<PyList>,
    bwlis: Py<PyList>,
    upstream: u64,
    downstream: u64,
    unscaled5prime: u64,
    unscaled3prime: u64,
    regionbodylength: u64,
    binsize: u64,
    missingdatazero: bool,
    referencepoint: &str,
    nproc: usize,
    verbose: bool,
    ofile: &str
) -> PyResult<()> {

    // Extract the bed and bigwig files from pyList to Vec.
    let mut bed_files: Vec<String> = Vec::new();
    let mut bw_files: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        bed_files = bedlis.extract(py).expect("Failed to retrieve bed file strings.");
        bw_files = bwlis.extract(py).expect("Failed to retrieve bed file strings.");
    });
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
    // Define the scaling regions in a struct
    let scale_regions = Scalingregions {
        upstream: upstream,
        downstream: downstream,
        unscaled5prime: unscaled5prime,
        unscaled3prime: unscaled3prime,
        regionbodylength: regionbodylength,
        binsize: binsize,
        cols_expected: ((bw_files.len() as u64 * bpsum) / binsize) as usize,
        missingdata_as_zero: missingdatazero,
        referencepoint: referencepoint.to_string(),
        mode: mode.to_string(),
        bwfiles: bw_files.len(),
        avgtype: "mean".to_string(),
        verbose: verbose,
        proc_number: nproc,
        bedlabels: bedlabels
    };
    // Parse regions from bed files. Note that we retain the name of the bed file (in case there are more then 1)
    // Additionaly, score and strand are also retained, if it's a 3-column bed file we just fill in '.'
    let (regions, regionsizes) = read_bedfiles(&bed_files);
    let slopregions = slop_regions(
        &regions,
        &scale_regions,
        &chromsizes
    );
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let matrix: Vec<Vec<f64>> = pool.install(|| {
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
    write_matrix(
        header_matrix(&scale_regions, regionsizes),
        matrix,
        ofile,
        regions
    );
    Ok(())
}

fn slop_regions(
    regions: &Vec<(String, u64, u64, String, String, String)>,
    scale_regions: &Scalingregions,
    chromsizes: &HashMap<String, u64>
) -> Vec<Vec<(u64, u64)>> {
    let mut regionranges: Vec<Vec<(u64, u64)>> = Vec::new();
    // Idea is to create a vector of tuples with start and end of every bin (binsize passed by computeMatrix).
    // The number of columns per region needs to be fixed per region.
    // Note that the before / after could mean that we run out of chromosome. 
    // Invalid regions (later to be encoded as NA or 0), will be pushed as (0,0) tuples.
    let col_expected = scale_regions.cols_expected / scale_regions.bwfiles;
    for region in regions.iter() {
        // To implement: 
        // // scale-regions + unscaled 5 and 3
        // // + and - encodings
        let chromend: u64 = *chromsizes.get(&region.0).unwrap();
        assert!(region.2 <= chromend, "Region end goes beyond chromosome boundary. Fix your bed files. {:?} > {}", region, chromend);
        assert!(region.1 <= chromend, "Region start goes beyond chromosome boundary. Fix your bed files. {:?} > {}", region, chromend);

        let anchorpoint = match scale_regions.referencepoint.as_str() {
            "TSS" => region.1,
            "TES" => region.2,
            "center" => (region.1 + region.2) / 2,
            _ => panic!("Reference should either be TSS, TES or center. {:?} is not supported.", scale_regions.referencepoint),
        };

        let mut regionsizes: Vec<(u64, u64)> = Vec::new();
        let mut absstart: i64 = anchorpoint as i64 - scale_regions.upstream as i64;
        let absstop: i64 = anchorpoint as i64 + scale_regions.downstream as i64;
        while absstart < absstop {
            let bin = absstart + scale_regions.binsize as i64;
            if absstart < 0  || bin > chromend as i64 {
                regionsizes.push((0,0));
            } else {
                regionsizes.push((absstart as u64, bin as u64))
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
    pub upstream: u64,
    pub downstream: u64,
    pub unscaled5prime: u64,
    pub unscaled3prime: u64,
    pub regionbodylength: u64,
    pub binsize: u64,
    pub cols_expected: usize,
    pub missingdata_as_zero: bool,
    pub referencepoint: String,
    pub mode: String,
    pub bwfiles: usize,
    pub avgtype: String,
    pub verbose: bool,
    pub proc_number: usize,
    pub bedlabels: Vec<String>
}