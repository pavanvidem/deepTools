use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::filehandler::bam_stats;
use crate::normalization::scale_factor_bamcompare;
use crate::covcalc::{bam_pileup, parse_regions};

#[pyfunction]
pub fn r_bamcompare(
    bam_ifile: &str,
    bam_ifile2: &str,
    _ofile: &str,
    _ofiletype: &str,
    _norm: &str,
    scalefactorsmethod: &str,
    effective_genome_size: u64,
    nproc: usize,
    binsize: u32,
    regions: Vec<(String, u64, u64)>,
    verbose: bool
) -> PyResult<()> {
    // put statistics into scope, this should probably be rewritten. (can't we always assume at least 2 threads ? )
    // will need to be revisited for multiBamsummaries / computeMatrix.
    let mut total_reads1: u64 = 0;
    let mut total_reads2: u64 = 0;
    let mut mapped_reads1: u64 = 0;
    let mut mapped_reads2: u64 = 0;
    let mut unmapped_reads1: u64 = 0;
    let mut unmapped_reads2: u64 = 0;
    let mut readlen1: f32 = 0.0;
    let mut readlen2: f32 = 0.0;
    let mut fraglen1: f32 = 0.0;
    let mut fraglen2: f32 = 0.0;
    // Get statistics of bam file.
    if nproc > 1 {
        let pool2 = ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let bamstatvec: Vec<_> = pool2.install(|| {
            vec![
                (bam_ifile, &verbose),
                (bam_ifile2, &verbose)
            ]
            .par_iter()
            .map(|(bam_ifile, verbose)| bam_stats(bam_ifile, verbose))
            .collect()
        });
        let (_total_reads1, _mapped_reads1, _unmapped_reads1, _readlen1, _fraglen1) = bamstatvec[0];
        let (_total_reads2, _mapped_reads2, _unmapped_reads2, _readlen2, _fraglen2) = bamstatvec[1];
        total_reads1 = _total_reads1;
        total_reads2 = _total_reads2;
        mapped_reads1 = _mapped_reads1;
        mapped_reads2 = _mapped_reads2;
        unmapped_reads1 = _unmapped_reads1;
        unmapped_reads2 = _unmapped_reads2;
        readlen1 = _readlen1;
        readlen2 = _readlen2;
        fraglen1 = _fraglen1;
        fraglen2 = _fraglen2;

    } else {
        let (_total_reads1, _mapped_reads1, _unmapped_reads1, _readlen1, _fraglen1) = bam_stats(bam_ifile, &verbose);
        let (_total_reads2, _mapped_reads2, _unmapped_reads2, _readlen2, _fraglen2) = bam_stats(bam_ifile2, &verbose);
        total_reads1 = _total_reads1;
        total_reads2 = _total_reads2;
        mapped_reads1 = _mapped_reads1;
        mapped_reads2 = _mapped_reads2;
        unmapped_reads1 = _unmapped_reads1;
        unmapped_reads2 = _unmapped_reads2;
        readlen1 = _readlen1;
        readlen2 = _readlen2;
        fraglen1 = _fraglen1;
        fraglen2 = _fraglen2;
    }

    // Calculate scale factors
    let (scale_factor1, scale_factor2) = scale_factor_bamcompare(
        scalefactorsmethod,
        mapped_reads1,
        mapped_reads2,
        binsize as u64,
        effective_genome_size
    );
    println!("scalefactors = {} and {}", scale_factor1, scale_factor2);
    let (regions, chromsizes)  = parse_regions(&regions, bam_ifile);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let _bg1: Vec<(String, u64, u64, f64)> = pool.install(|| {
        regions.par_iter()
            .flat_map(|i| bam_pileup(bam_ifile, &i, &binsize, scale_factor1))
            .collect()
    });
    let _bg2: Vec<(String, u64, u64, f64)> = pool.install(|| {
        regions.par_iter()
            .flat_map(|i| bam_pileup(bam_ifile2, &i, &binsize, scale_factor2))
            .collect()
    });
    Ok(())
}