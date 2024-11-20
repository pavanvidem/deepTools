use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::covcalc::{bam_pileup, parse_regions};
use crate::filehandler::{bam_stats, write_file};
use crate::normalization::scale_factor;

#[pyfunction]
pub fn r_bamcoverage(
    bam_ifile: &str,
    ofile: &str,
    ofiletype: &str,
    norm: &str,
    effective_genome_size: u64,
    nproc: usize,
    binsize: u32,
    regions: Vec<(String, u64, u64)>,
    _verbose: bool
) -> PyResult<()> {

    // Get statistics of bam file
    let (total_reads, mapped_reads, unmapped_reads, fraglen) = bam_stats(bam_ifile);
    println!("Total reads: {}", total_reads);
    println!("Mapped reads: {}", mapped_reads);
    println!("Unmapped reads: {}", unmapped_reads);
    println!("Fragment length: {}", fraglen);

    let scale_factor = scale_factor(norm, mapped_reads, binsize as u64, effective_genome_size);
    println!("Scale factor: {} calculated for norm: {}", scale_factor, norm);

    // Parse regions & calculate coverage
    let (regions, chromsizes)  = parse_regions(&regions, bam_ifile);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let _bg: Vec<(String, u64, u64, f64)> = pool.install(|| {
        regions.par_iter()
            .flat_map(|i| bam_pileup(bam_ifile, &i, &binsize, scale_factor))
            .collect()
    });

    // Create output
    write_file(ofile, ofiletype, _bg, chromsizes).unwrap();
    Ok(())
}