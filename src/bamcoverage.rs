use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::compute::{bam_pileup, parse_regions};

#[pyfunction]
pub fn r_bamcoverage(
    bam_ifile: &str,
    bedgraph_ofile: &str,
    nproc: usize,
    binsize: u32,
    regions: Vec<(String, u64, u64)>,
    _verbose: bool
) -> PyResult<()> {
    
    let regions = parse_regions(&regions, bam_ifile);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let _bg: Vec<(String, u64, u64, f64)> = pool.install(|| {
        regions.par_iter()
            .flat_map(|i| bam_pileup(bam_ifile, &i, &binsize))
            .collect()
    });
    // write bedgraph
    let mut writer = BufWriter::new(File::create(bedgraph_ofile).unwrap());
    for i in _bg {
        writeln!(writer, "{}\t{}\t{}\t{}", i.0, i.1, i.2, i.3).unwrap();
    }
    Ok(())
}