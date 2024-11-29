use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::filehandler::{bam_ispaired, write_file};
use crate::covcalc::{bam_pileup, parse_regions, alignmentfilters};
use crate::normalization::scale_factor_bamcompare;
use crate::calc::median;

#[pyfunction]
pub fn r_bamcompare(
    // input and output
    bam_ifile1: &str, // input bamfile 1
    bam_ifile2: &str, // input bamfile 2
    ofile: &str, // output file
    ofiletype: &str, // ouput file type, bedgraph or bigwig
    // norm options
    norm: &str,
    effective_genome_size: u64,
    scalefactorsmethod: &str,
    operation: &str,
    pseudocount: f32,
    // filtering options
    ignoreduplicates: bool,
    minmappingquality: u8, // 
    samflaginclude: u16,
    samflagexclude: u16,
    minfraglen: u32,
    maxfraglen: u32,
    nproc: usize,
    _ignorechr: Py<PyList>,
    binsize: u32,
    regions: Vec<(String, u32, u32)>,
    verbose: bool
) -> PyResult<()> {
    let ispe1 = bam_ispaired(bam_ifile1);
    let ispe2 = bam_ispaired(bam_ifile2);

    if verbose {
        println!("Sample1: {} is-paired: {}", bam_ifile1, ispe1);
        println!("Sample2: {} is-paired: {}", bam_ifile2, ispe2);
    }
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        ignorechr = _ignorechr.extract(py).expect("Failed to retrieve ignorechr.");
    });
    // Set alignment filters
    let filters = alignmentfilters {
        minmappingquality: minmappingquality,
        samflaginclude: samflaginclude,
        samflagexclude: samflagexclude,
        minfraglen: minfraglen,
        maxfraglen: maxfraglen
    };

    // Parse regions & calculate coverage
    let (regions, chromsizes)  = parse_regions(&regions, bam_ifile1);
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    
    // // Parse first bamfile
    // let (bg1, mapped1, _unmapped1, readlen1, fraglen1) = pool.install(|| {
    //     regions.par_iter()
    //         .map(|i| bam_pileup(bam_ifile1, &i, &binsize, &ispe1, &ignorechr, &filters))
    //         .reduce(
    //             || (vec![], 0, 0, vec![], vec![]),
    //             |(mut _bg, mut _mapped, mut _unmapped, mut _readlen, mut _fraglen), (bg, mapped, unmapped, readlen, fraglen)| {
    //                 _bg.extend(bg);
    //                 _readlen.extend(readlen);
    //                 _fraglen.extend(fraglen);
    //                 _mapped += mapped;
    //                 _unmapped += unmapped;
    //                 (_bg, _mapped, _unmapped, _readlen, _fraglen)
    //             }
    //         )
    // });
    // let _readlen1 = median(readlen1);
    // let _fraglen1 = median(fraglen1);

    // // Parse second bamfile
    // let (bg2, mapped2, _unmapped2, readlen2, fraglen2) = pool.install(|| {
    //     regions.par_iter()
    //         .map(|i| bam_pileup(bam_ifile2, &i, &binsize, &ispe2, &ignorechr, &filters))
    //         .reduce(
    //             || (vec![], 0, 0, vec![], vec![]),
    //             |(mut _bg, mut _mapped, mut _unmapped, mut _readlen, mut _fraglen), (bg, mapped, unmapped, readlen, fraglen)| {
    //                 _bg.extend(bg);
    //                 _readlen.extend(readlen);
    //                 _fraglen.extend(fraglen);
    //                 _mapped += mapped;
    //                 _unmapped += unmapped;
    //                 (_bg, _mapped, _unmapped, _readlen, _fraglen)
    //             }
    //         )
    // });
    // let _readlen2 = median(readlen2);
    // let _fraglen2 = median(fraglen2);
    // let (sf1, sf2) = scale_factor_bamcompare(scalefactorsmethod, mapped1, mapped2, binsize, effective_genome_size, norm);
    // println!("scale factor1 = {}, scale factor2 = {}", sf1, sf2);
    // let bge = collapse_bgvecs(bg1, bg2, sf1, sf2, pseudocount, operation);
    // write_file(ofile, ofiletype, bge, chromsizes);
    Ok(())
}