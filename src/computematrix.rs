use pyo3::prelude::*;
use pyo3::types::PyList;
use crate::filehandler::{read_bedfile, read_gtffile, chrombounds_from_bw, bwintervals, header_matrix, write_matrix};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;
use std::fmt;
use itertools::Itertools;
use ndarray::Array1;
use crate::calc::{mean_float, median_float, max_float, min_float, sum_float};


#[pyfunction]
pub fn r_computematrix(
    mode: &str, // reference-point or scale-regions
    regionlis: Py<PyList>, // python list of region files (bed or gtf)
    bwlis: Py<PyList>, // python list of bigwig files
    sampleslabel: Py<PyList>, // python list of sample labels, if empty, use bigwig file names.
    upstream: u32, // upstream region to consider
    downstream: u32, // downstream region to consider
    unscaled5prime: u32, // unscaled region 5' of the anchorpoint, only used in scale-regions mode.
    unscaled3prime: u32, // unscaled region 3' of the anchorpoint, only used in scale-regions mode.
    regionbodylength: u32, // length of the region body (after scaling), only used in scale-regions mode.
    binsize: u32, // binsize to use for the matrix
    missingdatazero: bool, // Encode missing data as 0. Default is false (and will be encoded as NA).
    metagene: bool, // If set, 'exons' are stitched together to form a metagene
    txnid: &str, // transcript id to use when parsing GTF file
    exonid: &str, // exon id to use when parsing GTF file
    txniddesignator: &str, // designator to use when parsing GTF file
    scale: f32, // scaling factor for writing out values. default is 1.0 (no scaling)
    nanafterend: bool, // end regions will treated as nans. Default is false.
    skipzeros: bool, // skip regions with all zeros. Default is false.
    minthresh: f32, // minimum threshold to keep a region. If not set it will equal 0.0
    maxthresh: f32, // maximum threshold to keep a region. if not set it will equal 0.0
    averagetypebins: &str, // operation to summarize values over bins. Default is mean.
    sortregions: &str, // either ascend, descend or keep. Default is keep (and ignores sortusing).
    sortusing: &str, // metric to sort on. Either mean median max min sum region_length. Default is mean.
    sortusingsamples: Py<PyList>, // list of samples to sort on. If empty, use all samples.
    referencepoint: &str, // reference point to use. Either TSS, TES or center. Default is TSS. Only used in reference-point mode.
    nproc: usize, // number of threads.
    verbose: bool, // verbose output.
    ofile: &str // npz file to write to.
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
        bpsum: bpsum,
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
    let gtfparse = Gtfparse {
        metagene: metagene,
        txnid: txnid.to_string(),
        exonid: exonid.to_string(),
        txniddesignator: txniddesignator.to_string(),
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
    // Define slop regions, which contain the actual 'bins' to query the bigwig files.
    let slopregions = pool.install(|| {
        regions.par_iter()
            .map(|region| slop_region(&region, &scale_regions, &chromsizes))
            .collect::<Vec<_>>()
    });

    // Discriminate between reference-point and scale-regions mode.

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
    matrix_dump(
        sortregions,
        sortusing,
        sort_using_samples,
        regions,
        matrix, 
        scale_regions,
        regionsizes,
        ofile,
        verbose
    );

    Ok(())
}

fn slop_region(
    region: &Region,
    scale_regions: &Scalingregions,
    chromsizes: &HashMap<String, u32>
) -> Vec<Bin> {
    // Idea is to create a vector Bins (Conbin or Catbin) which encodes start and end of every bin (binsize passed by computeMatrix).
    // Catbin takes care of the situation where one needs metagenes, and thus multiple start/end per bin are possible.
    // The number of columns is predetermined
    // Note that the before / after could mean that we run out of chromosome. 
    // Invalid regions (later to be encoded as NA or 0), will be pushed as (0,0) tuples.
    // Note that if nan_after_end is set to true, we will push (0,0) tuples after the end of the region.
    
    // Get the chromosome end for a specific region, and assert that the region stays within the chromosome boundary.
    // Note that only a right check is needed, as positions are u32.
    // Note that we know &region.chrom is inside chromsizes already, since this filtering is done at the region reading stage.
    let chromend: u32 = *chromsizes.get(&region.chrom).unwrap();
    region.assert_end(chromend);
    region.get_anchor_bins(scale_regions, chromend)
}

#[derive(Debug)]
pub struct Scalingregions {
    pub upstream: u32,
    pub downstream: u32,
    pub unscaled5prime: u32,
    pub unscaled3prime: u32,
    pub regionbodylength: u32,
    pub binsize: u32,
    pub cols_expected: usize,
    pub bpsum: u32,
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

#[derive(Clone, Debug)]
pub enum Bin {
    Conbin(u32, u32),
    Catbin(Vec<(u32, u32)>),
}

impl Bin {
    pub fn get_start(&self) -> u32 {
        match self {
            Bin::Conbin(start, _) => *start,
            Bin::Catbin(starts) => starts.first().unwrap().0,
        }
    }
    pub fn get_end(&self) -> u32 {
        match self {
            Bin::Conbin(_, end) => *end,
            Bin::Catbin(ends) => ends.last().unwrap().1,
        }
    }
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

    pub fn get_anchor_bins(&self, scale_regions: &Scalingregions, chromend: u32) -> Vec<Bin> {
        // Given an anchorpoint, return a vector, start, end , middle
        // The order of the vector is always 5' -> 3', meaning 'increasing' for +/. regions, and 'decreasing' for - regions.
        // At this stage, two situations are possible:
        // - self.start / self.end are Revalue::U, meaning we are in 'non metagene' mode.
        // - self.start / self.end are Revalue::V, meaning we are in 'metagene' mode, and the bins returned are exon-aware.
        // We need a notion of bins that don't make sense (i.e. beyond chromosome boundaries). These are encoded as (0,0)

        let mut bins: Vec<Bin> = Vec::new();
        let mut bodybins: Vec<Bin> = Vec::new();

        // Define anchorpoints
        let anchorstart;
        let anchorstop;
        match scale_regions.mode.as_str() {
            "reference-point" => {
                anchorstart = self.get_anchorpoint(&scale_regions.referencepoint);
                anchorstop = anchorstart;
            },
            "scale-regions" => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        anchorstart = *start;
                        anchorstop = *end;
                    },
                    (Revalue::V(start), Revalue::V(end)) => {
                        anchorstart = *start.first().unwrap();
                        anchorstop = *end.last().unwrap();
                    },
                    _ => panic!("Start and End are not either both u32, or Vecs. This means your regions file is ill-defined. Fix {}.",self.name),
                }
            },
            _ => panic!("Mode should either be reference-point or scale-regions. {} is not supported.", scale_regions.mode),
        }
        if scale_regions.mode != "reference-point" {
            // scale-regions mode. Assert
            assert!(scale_regions.regionbodylength != 0, "scale-regions mode, but regionbodylength is 0.");
            if self.regionlength < (scale_regions.unscaled5prime + scale_regions.unscaled3prime) ||
               self.regionlength - (scale_regions.unscaled5prime + scale_regions.unscaled3prime) < scale_regions.binsize {
                println!("Warning ! Region {} is shorter than the binsize (potentially after unscaled regions taken into account. Whole region encoded as 0 or NA", self.name);
                let nbin = scale_regions.cols_expected / scale_regions.bwfiles;
                for _ in 0..nbin {
                    bins.push(Bin::Conbin(0,0));
                }
                return bins;
            } else {
                bodybins.extend(self.scale_regionbody(scale_regions, chromend));
            }
        }

        // Get flanking regions.
        // Note that we still need to deal with exon - non-exon as reference-point mode could require metagene walks.
        match self.strand.as_str() {
            "+" | "." => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut rightbins: Vec<Bin> = Vec::new();

                        let mut absstart: i64 = anchorstart as i64 - scale_regions.upstream as i64;
                        let absstop: i64 = anchorstop as i64 + scale_regions.downstream as i64;

                        for binix in (absstart..anchorstart as i64).step_by(scale_regions.binsize as usize) {
                            if binix < 0 || binix as u32 >= chromend || (binix + scale_regions.binsize as i64) as u32 >= chromend {
                                leftbins.push(Bin::Conbin(0,0));
                            } else if scale_regions.nan_after_end && binix as u32 <= *start  {
                                leftbins.push(Bin::Conbin(0,0));
                            } else {
                                leftbins.push(Bin::Conbin(binix as u32, (binix as u32) + scale_regions.binsize));
                            }
                        }

                        for binix in (anchorstop as i64..absstop).step_by(scale_regions.binsize as usize) {
                            if binix < 0 || binix as u32 >= chromend || (binix + scale_regions.binsize as i64) as u32 >= chromend {
                                rightbins.push(Bin::Conbin(0,0));
                            } else if scale_regions.nan_after_end && binix as u32 >= *end {
                                rightbins.push(Bin::Conbin(0,0));
                            } else {
                                rightbins.push(Bin::Conbin(binix as u32, (binix as u32) + scale_regions.binsize));
                            }
                        }
                        
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                        // If we have bodybins, they should be squeezed in here.
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstop;
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
                        let mut lastanchor: u32 = anchorstart;
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
                                    scale_regions.nan_after_end,
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
                        // If we have bodybins, they should be squeezed in here.
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
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
                        let mut leftbins: Vec<Bin> = Vec::new();
                        let mut rightbins: Vec<Bin> = Vec::new();

                        let mut absstart: i64 = anchorstop as i64 + scale_regions.upstream as i64;
                        let absstop: i64 = anchorstart as i64 - scale_regions.downstream as i64;

                        let steps: Vec<_> = (anchorstop as i64..absstart)
                            .step_by(scale_regions.binsize as usize)
                            .collect();
                        for binix in steps.into_iter().rev() {
                            if binix as u32 > chromend || (binix + scale_regions.binsize as i64) as u32 > chromend {
                                rightbins.push(Bin::Conbin(0,0));
                            } else if scale_regions.nan_after_end && binix as u32 >= *end {
                                leftbins.push(Bin::Conbin(0,0));
                            } else {
                                leftbins.push(Bin::Conbin(binix as u32, (binix as u32) + scale_regions.binsize));
                            }
                        }
                        
                        let steps: Vec<_> = (absstop..anchorstart as i64)
                            .step_by(scale_regions.binsize as usize)
                            .collect();
                        for binix in steps.into_iter().rev() {
                            if binix < 0 {
                                leftbins.push(Bin::Conbin(0,0));
                            } else if scale_regions.nan_after_end && binix as u32 + scale_regions.binsize <= *start {
                                rightbins.push(Bin::Conbin(0,0));
                            } else {
                                rightbins.push(Bin::Conbin(binix as u32, (binix as u32) + scale_regions.binsize));
                            }
                        }
                        
                        for bin in rightbins.into_iter() {
                            bins.push(bin);
                        }
                        // If we have bodybins, they should be squeezed in here.
                        bodybins.reverse();
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        // Reset bodybins, as they are consumed.
                        bodybins = Vec::new();
                        for bin in leftbins.into_iter() {
                            bins.push(bin);
                        }
                    }
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        // Right side.
                        let mut rightbins: Vec<Bin> = Vec::new();
                        let mut lastanchor: u32 = anchorstop;
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
                        let mut lastanchor: u32 = anchorstart;
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
                        bodybins.reverse();
                        for bin in bodybins.into_iter() {
                            bins.push(bin);
                        }
                        bodybins = Vec::new();
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

    pub fn scale_regionbody(&self, scale_regions: &Scalingregions, chromend: u32) -> Vec<Bin> {
        // A vector of bins needs to be constructed for regionbody.
        // Since this is scaling mode, 'linspace' functionality is reproduced.
        match self.strand.as_str() {
            "+" | "." => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // No exons, forward strand. divide start  - end as such:
                        // |---un5prime---|---bodylength---|---un3prime---|
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        let mut innerbins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            un5bins.extend((0..scale_regions.unscaled5prime)
                                .step_by(scale_regions.binsize as usize)
                                .map(|i| Bin::Conbin(*start + i, *start + i + scale_regions.binsize))
                                .collect::<Vec<Bin>>());
                        }
                        
                        if scale_regions.unscaled3prime > 0 {
                            un3bins.extend( (0..scale_regions.unscaled3prime)
                                .step_by(scale_regions.binsize as usize)
                                .rev()
                                .map(|i| Bin::Conbin(*end - i - scale_regions.binsize, *end - i))
                                .collect::<Vec<Bin>>() );
                        }
                        let bodystart = *start + scale_regions.unscaled5prime;
                        let bodyend = *end - scale_regions.unscaled3prime;
                        
                        // Get the bins over the body length. These need to be scaled, so similar to deeptools < 4, linspace is used.
                        let neededbins = (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        // There's multiple options here:
                        // transcriptlength >= regionbodylength -> linspace
                        // regionbodylength / binsize > transcriptlength <= regionbodylength -> 1 >= binsize > binsize.
                        // transcriptlength <= regionbodylength / binsize -> index repetitions with binsize of one.
                        let scaledbinsize = std::cmp::min(std::cmp::max((bodyend - bodystart) / neededbins as u32, 1), scale_regions.binsize);

                        innerbins.extend( Array1::linspace(bodystart as f32, (bodyend - scaledbinsize) as f32, neededbins)
                            .mapv(|x| x as u32)
                            .map(|x| Bin::Conbin(*x, *x + scaledbinsize))
                            .into_iter()
                            .collect::<Vec<_>>() );

                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        return combined_bins;
                    },
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = start[0];

                            while walked_bps < scale_regions.unscaled5prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true
                                );
                                un5bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        if scale_regions.unscaled3prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = *end.last().unwrap();

                            while walked_bps < scale_regions.unscaled3prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false
                                );
                                un3bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        un3bins.reverse();
                        
                        let bodystart: u32;
                        let bodyend: u32;
                        if scale_regions.unscaled5prime > 0 {
                            bodystart = un5bins.last().unwrap().get_end();
                        } else {
                            bodystart = *start.first().unwrap();
                        }
                        if scale_regions.unscaled3prime > 0 {
                            bodyend = un3bins.first().unwrap().get_start();
                        } else {
                            bodyend = *end.last().unwrap();
                        }
                        let truebodylength = self.regionlength - scale_regions.unscaled5prime - scale_regions.unscaled3prime;
                        let neededbins = (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        let scaledbinsize = std::cmp::min(std::cmp::max(truebodylength / neededbins as u32, 1), scale_regions.binsize);
                        // Things are a bit tricky now, as we can do a linspace over the region, but we don't have a notion of the exons.
                        // I think easiest is to just pull a hashmap over the entire region, get linspace from hashmap to vec, and be done with it.
                        // technically we fetch a bunch of regions we don't need, but this operation is not too expensive.

                        let mut binmap: HashMap<u32, Bin> = HashMap::new();
                        let mut lastanchor: u32 = bodystart;

                        for ix in 0..((truebodylength/scaledbinsize)+1) {
                            let (bin, anchor) = refpoint_exonwalker(
                                &exons,
                                lastanchor,
                                scaledbinsize,
                                chromend,
                                scale_regions.nan_after_end,
                                true
                            );
                            lastanchor = anchor;
                            match bin {
                                Bin::Conbin(start, end) => {
                                    if end > bodyend {
                                        binmap.insert(ix, Bin::Conbin(start, bodyend));
                                    } else {
                                        binmap.insert(ix, Bin::Conbin(start, end));
                                    }
                                },
                                Bin::Catbin(bins) => {
                                    if bins.last().unwrap().1 > bodyend {
                                        let mut newbins: Vec<(u32, u32)> = Vec::new();
                                        for (s, e) in bins.iter() {
                                            if *e > bodyend {
                                                newbins.push((*s, bodyend));
                                            } else {
                                                newbins.push((*s, *e));
                                            }
                                        }
                                        binmap.insert(ix, Bin::Catbin(newbins));
                                    } else {
                                        binmap.insert(ix, Bin::Catbin(bins));
                                    }
                                },
                            }
                        }

                        let innerbins = Array1::linspace(0 as f32, ((truebodylength)/scaledbinsize) as f32, neededbins)
                            .mapv(|x| x as u32)
                            .map(|x| binmap.get(&x).unwrap().clone())
                            .into_iter()
                            .collect::<Vec<Bin>>();

                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        return combined_bins;
                    },
                    _ => panic!("Start and End are not either both u32, or Vecs. This means your regions file is ill-defined."),
                }
            },
            "-" => {
                match (&self.start, &self.end) {
                    (Revalue::U(start), Revalue::U(end)) => {
                        // No exons, negative strand. divide start  - end as such:
                        // |---un3prime---|---bodylength---|---un5prime---|
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        let mut innerbins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            un5bins.extend((0..scale_regions.unscaled5prime)
                                .step_by(scale_regions.binsize as usize)
                                .map(|i| Bin::Conbin(*end - i - scale_regions.binsize, *end - i))
                                .collect::<Vec<Bin>>());
                        }
                        
                        if scale_regions.unscaled3prime > 0 {
                            un3bins.extend( (0..scale_regions.unscaled3prime)
                                .step_by(scale_regions.binsize as usize)
                                .rev()
                                .map(|i| Bin::Conbin(*start + scale_regions.unscaled3prime - i - scale_regions.binsize, *start + scale_regions.unscaled3prime - i))
                                .collect::<Vec<Bin>>() );
                        }
                        let bodystart = *start + scale_regions.unscaled3prime;
                        let bodyend = *end - scale_regions.unscaled5prime;

                        // Get the bins over the body length. These need to be scaled, so similar to deeptools < 4, linspace is used.
                        let neededbins = (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        // There's multiple options here:
                        // transcriptlength >= regionbodylength -> linspace
                        // regionbodylength / binsize > transcriptlength <= regionbodylength -> 1 >= binsize > binsize.
                        // transcriptlength <= regionbodylength / binsize -> index repetitions with binsize of one.
                        let mut scaledbinsize = (bodyend - bodystart)/neededbins as u32;
                        if scaledbinsize == 0 {
                            scaledbinsize = 1;
                        }
                        if scaledbinsize > scale_regions.binsize {
                            scaledbinsize = scale_regions.binsize;
                        }

                        innerbins.extend( Array1::linspace(bodystart as f32, (bodyend - scaledbinsize) as f32, neededbins)
                            .mapv(|x| x as u32)
                            .map(|x| Bin::Conbin(*x, *x + scaledbinsize))
                            .into_iter()
                            .collect::<Vec<_>>() );
                        
                        // Reverse innerbins to go from 3' -> 5'
                        innerbins.reverse();
                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        return combined_bins;
                    },
                    (Revalue::V(start), Revalue::V(end)) => {
                        let exons: Vec<(u32, u32)> = start.iter().zip(end.iter())
                            .map(|(&s, &e)| (s, e)) 
                            .collect();
                        let mut un5bins: Vec<Bin> = Vec::new();
                        let mut un3bins: Vec<Bin> = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = *end.last().unwrap();

                            while walked_bps < scale_regions.unscaled5prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    false
                                );
                                un5bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        un5bins.reverse();

                        if scale_regions.unscaled3prime > 0 {
                            let mut walked_bps: u32 = 0;
                            let mut lastanchor: u32 = start[0];

                            while walked_bps < scale_regions.unscaled3prime {
                                let (bin, retanch) = refpoint_exonwalker(
                                    &exons,
                                    lastanchor,
                                    scale_regions.binsize,
                                    chromend,
                                    scale_regions.nan_after_end,
                                    true
                                );
                                un3bins.push(bin);
                                walked_bps += scale_regions.binsize;
                                lastanchor = retanch;
                            }
                        }
                        
                        let bodystart: u32;
                        let bodyend: u32;
                        if scale_regions.unscaled3prime > 0 {
                            bodystart = un3bins.last().unwrap().get_end();
                        } else {
                            bodystart = *start.first().unwrap();
                        }
                        if scale_regions.unscaled5prime > 0 {
                            bodyend = un5bins.first().unwrap().get_start();
                        } else {
                            bodyend = *end.last().unwrap();
                        }
                        let truebodylength = self.regionlength - scale_regions.unscaled5prime - scale_regions.unscaled3prime;
                        let neededbins = (scale_regions.regionbodylength / scale_regions.binsize) as usize;
                        let scaledbinsize = std::cmp::min(std::cmp::max(truebodylength / neededbins as u32, 1), scale_regions.binsize);
                        // Things are a bit tricky now, as we can do a linspace over the region, but we don't have a notion of the exons.
                        // I think easiest is to just pull a hashmap over the entire region, get linspace from hashmap to vec, and be done with it.
                        // technically we fetch a bunch of regions we don't need, but this operation is not too expensive.

                        let mut binmap: HashMap<u32, Bin> = HashMap::new();
                        let mut lastanchor: u32 = bodystart;

                        for ix in 0..((truebodylength/scaledbinsize)+1) {
                            let (bin, anchor) = refpoint_exonwalker(
                                &exons,
                                lastanchor,
                                scaledbinsize,
                                chromend,
                                scale_regions.nan_after_end,
                                true
                            );
                            lastanchor = anchor;
                            match bin {
                                Bin::Conbin(start, end) => {
                                    if end > bodyend {
                                        binmap.insert(ix, Bin::Conbin(start, bodyend));
                                    } else {
                                        binmap.insert(ix, Bin::Conbin(start, end));
                                    }
                                },
                                Bin::Catbin(bins) => {
                                    if bins.last().unwrap().1 > bodyend {
                                        let mut newbins: Vec<(u32, u32)> = Vec::new();
                                        for (s, e) in bins.iter() {
                                            if *e > bodyend {
                                                newbins.push((*s, bodyend));
                                            } else {
                                                newbins.push((*s, *e));
                                            }
                                        }
                                        binmap.insert(ix, Bin::Catbin(newbins));
                                    } else {
                                        binmap.insert(ix, Bin::Catbin(bins));
                                    }
                                },
                            }
                        }

                        let innerbins = Array1::linspace(0 as f32, ((truebodylength)/scaledbinsize) as f32, neededbins)
                            .mapv(|x| x as u32)
                            .map(|x| binmap.get(&x).unwrap().clone())
                            .into_iter()
                            .collect::<Vec<Bin>>();
                        // Combine the vectors and return
                        let mut combined_bins = Vec::new();
                        if scale_regions.unscaled5prime > 0 {
                            combined_bins.extend(un5bins.into_iter());
                        }
                        combined_bins.extend(innerbins.into_iter());
                        if scale_regions.unscaled3prime > 0 {
                            combined_bins.extend(un3bins.into_iter());
                        }
                        return combined_bins;
                    },
                    _ => panic!("Start and End are not either both u32, or Vecs. This means your regions file is ill-defined."),
                }
            },
            _ => panic!("Strand should either be + or - or . {:?} is not supported.", self.strand),
        }
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

fn matrix_dump(
    sortregions: &str,
    sortusing: &str,
    sort_using_samples: Vec<u32>,
    regions: Vec<Region>,
    matrix: Vec<Vec<f32>>,
    scale_regions: Scalingregions,
    regionsizes: HashMap<String, u32>,
    ofile: &str,
    verbose: bool
) {
    // Takes a pre-computed matrix, resorts it if requested, and writes it to file.
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
                let start = (sample_ix - 1) * scale_regions.bpsum;
                let end = start + scale_regions.bpsum;
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

        // Reorder matrix & regions
        let sortedmatrix = sortedix
            .iter()
            .map(|ix| matrix[*ix].clone())
            .collect();
        let sortedregions = sortedix
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
}