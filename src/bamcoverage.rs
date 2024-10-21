use pyo3::prelude::*;
use rust_htslib::{bam, bam::Read, bam::IndexedReader};
use std::fs::File;
use std::io::{BufWriter, Write};


#[pyfunction]
pub fn r_bamcoverage(bam_ifile: &str, bedgraph_ofile: &str) -> PyResult<()> {
    // files
    let mut bam = bam::IndexedReader::from_path(bam_ifile)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM file: {}", e)))?;
    let output_file = File::create(bedgraph_ofile)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output file: {}", e)))?;
    let mut writer = BufWriter::new(output_file);

    let mut l_chr = String::new();
    let mut l_start: u64 = 0;
    let mut l_end: u64 = 0;
    let mut l_cov: u64 = 0;
    let mut l_chromi: u32 = 0;
    let header = bam.header().clone();

    for pileup in bam.pileup() {

        let pileup = pileup 
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error in pileup: {}", e)))?;
        let chrom = String::from_utf8(header.tid2name(pileup.tid() as u32).to_owned())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid UTF-8 in chromosome name: {}", e)))?;
        let chromi = pileup.tid() as u32;
        let pos = pileup.pos() as u64;
        let cov = pileup.depth() as u64;

        // edge cases
        if chrom != l_chr {
            // catch case that last chromosome was 0 at end
            // This is no longer needed when parallel chunking over chromosomes
            if !l_chr.is_empty() {
                let chrend = header.target_len(l_chromi)
                    .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Chromosome {} not found in the BAM file header.", l_chr)))?;
                if l_end + 1 < chrend {
                    writeln!(writer, "{}\t{}\t{}\t{}", l_chr, l_end + 1, chrend, 0)?;
                }
            }
            // catch case that beginning of chromosome is 0
            if pos > 0 {
                writeln!(writer, "{}\t{}\t{}\t{}", chrom, 0, pos, 0)?;
            }
            l_chromi = chromi;
            l_chr = chrom;
            l_start = pos;
            l_cov = cov;
        } else {
            if pos != l_end + 1 {
                writeln!(writer, "{}\t{}\t{}\t{}", l_chr, l_start, l_end + 1, l_cov)?;
                writeln!(writer, "{}\t{}\t{}\t{}", l_chr, l_end + 1, pos, 0)?;
                l_start = pos;
                l_cov = cov;
            } else if l_cov != cov {
                writeln!(writer, "{}\t{}\t{}\t{}", l_chr, l_start, pos, l_cov)?;
                l_start = pos;
            }
        }
        l_end = pos;
        l_cov = cov;
    }
    Ok(())
}