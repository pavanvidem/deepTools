use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::covcalc::{bam_pileup, parse_regions};
use crate::bamhandler::bam_stats;
use crate::normalization::scale_factor;

#[pyfunction]
pub fn r_bamcoverage(
    bam_ifile: &str,
    bedgraph_ofile: &str,
    norm: &str,
    effective_genome_size: u64,
    nproc: usize,
    binsize: u32,
    regions: Vec<(String, u64, u64)>,
    _verbose: bool
) -> PyResult<()> {
    // Get statistics of bam file
    let (total_reads, mapped_reads, unmapped_reads) = bam_stats(bam_ifile);
    println!("Total reads: {}", total_reads);
    println!("Mapped reads: {}", mapped_reads);
    println!("Unmapped reads: {}", unmapped_reads);

    // Calculate scale factor
    let scale_factor = scale_factor(norm, mapped_reads, binsize as u64, effective_genome_size);
    println!("Scale factor: {} calculated for norm: {}", scale_factor, norm);

    // Parse regions & calculate coverage
    let regions = parse_regions(&regions, bam_ifile);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let _bg: Vec<(String, u64, u64, f64)> = pool.install(|| {
        regions.par_iter()
            .flat_map(|i| bam_pileup(bam_ifile, &i, &binsize, scale_factor))
            .collect()
    });
    // write bedgraph
    let mut writer = BufWriter::new(File::create(bedgraph_ofile).unwrap());
    for i in _bg {
        writeln!(writer, "{}\t{}\t{}\t{}", i.0, i.1, i.2, i.3).unwrap();
    }
    Ok(())
}