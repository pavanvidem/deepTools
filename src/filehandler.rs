use rust_htslib::bam::{Read, Reader};
use std::io::{BufWriter, Write};
use std::fs::File;
use bigtools::{BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use std::collections::HashMap;

use crate::calc::median;

pub fn bam_stats(bam_ifile: &str) -> (u64, u64, u64, f32) {
    let mut fraglens: Vec<u32> = Vec::new();
    let ispe = bam_ispaired(bam_ifile);
    let mut bam = Reader::from_path(bam_ifile).unwrap();
    let mut total_reads: u64 = 0;
    let mut mapped_reads: u64 = 0;
    let mut unmapped_reads: u64 = 0;
    for record in bam.records() {
        total_reads += 1;
        let record = record.expect("Error parsing record.");
        if record.is_unmapped() {
            unmapped_reads += 1;
        } else {
            mapped_reads += 1;
            if ispe {
                if record.is_paired() && record.is_proper_pair() {
                    let flen = record.insert_size() as u32;
                    if flen > 0 {
                        fraglens.push(flen);
                    }
                }
            } else {
                fraglens.push(record.seq().len() as u32);
            }
        }
    }
    return (total_reads, mapped_reads, unmapped_reads, median(fraglens));
}

fn bam_ispaired(bam_ifile: &str) -> bool {
    let mut bam = Reader::from_path(bam_ifile).unwrap();
    for record in bam.records() {
        let record = record.expect("Error parsing record.");
        if record.is_paired() {
            return true;
        }
    }
    return false;
}

pub fn write_file(ofile: &str, filetype: &str, bg: Vec<(String, u64, u64, f64)>, chromsizes: HashMap<String, u32>) -> Result<(), std::string::String> {
    if filetype == "bedgraph" {
        // write output file, bedgraph
        let mut writer = BufWriter::new(File::create(ofile).unwrap());
        for i in bg {
            writeln!(writer, "{}\t{}\t{}\t{}", i.0, i.1, i.2, i.3).unwrap();
        }
    } else {
        // write output file, bigwig
        let vals = BedParserStreamingIterator::wrap_infallible_iter(
            bg.iter().map(
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