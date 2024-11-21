use pyo3::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::covcalc::{bam_pileup, parse_regions, collapse_bgvec};
use crate::filehandler::{bam_ispaired, write_file};
use crate::normalization::scale_factor;
use crate::calc::median;

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
    verbose: bool
) -> PyResult<()> {
    let ispe = bam_ispaired(bam_ifile);
    if verbose {
        println!("Sample: {} is-paired: {}", bam_ifile, ispe);
    }
    // Parse regions & calculate coverage
    let (regions, chromsizes)  = parse_regions(&regions, bam_ifile);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let (bg, mapped, _unmapped, readlen, fraglen) = pool.install(|| {
        regions.par_iter()
            .map(|i| bam_pileup(bam_ifile, &i, &binsize, &ispe))
            .reduce(
                || (vec![], 0, 0, vec![], vec![]),
                |(mut _bg, mut _mapped, mut _unmapped, mut _readlen, mut _fraglen), (bg, mapped, unmapped, readlen, fraglen)| {
                    _bg.extend(bg);
                    _readlen.extend(readlen);
                    _fraglen.extend(fraglen);
                    _mapped += mapped;
                    _unmapped += unmapped;
                    (_bg, _mapped, _unmapped, _readlen, _fraglen)
                }
            )
    });
    let readlen = median(readlen);
    let fraglen = median(fraglen);
    let sf = scale_factor(
        norm, 
        mapped,
        binsize,
        effective_genome_size,
        readlen,
        fraglen,
        &verbose
    );
    let bg_scaled = collapse_bgvec(bg, sf);
    // Create output
    write_file(ofile, ofiletype, bg_scaled, chromsizes);
    Ok(())
}