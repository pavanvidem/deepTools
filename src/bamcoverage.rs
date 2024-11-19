use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{BufWriter, Write};
use bigtools::{BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use crate::covcalc::{bam_pileup, parse_regions};
use crate::filehandler::bam_stats;
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
    if ofiletype == "bedgraph" {
        // write output file, bedgraph
        let mut writer = BufWriter::new(File::create(ofile).unwrap());
        for i in _bg {
            writeln!(writer, "{}\t{}\t{}\t{}", i.0, i.1, i.2, i.3).unwrap();
        }
    } else {
        // write output file, bigwig
        let vals = BedParserStreamingIterator::wrap_infallible_iter(
            _bg.iter().map(
                |(chr, start, end, cov)| {(chr.as_str(), Value {start: *start as u32, end: *end as u32, value: *cov as f32 } )}
            ),
            true
        );
        // Theoretically one could add more threads here too, but this would require rewrite of the _bg iter upstream.
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()
            .expect("Unable to create tokio runtime for bw writing.");
        let writer = BigWigWrite::create_file(ofile, chromsizes).unwrap();
        let _ = writer.write(vals, runtime);
    }
    Ok(())
}