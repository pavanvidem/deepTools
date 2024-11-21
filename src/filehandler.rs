use rust_htslib::bam::{Read, Reader};
use std::io::{BufWriter, Write};
use std::fs::File;
use bigtools::{BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use std::collections::HashMap;

pub fn bam_ispaired(bam_ifile: &str) -> bool {
    let mut bam = Reader::from_path(bam_ifile).unwrap();
    for record in bam.records() {
        let record = record.expect("Error parsing record.");
        if record.is_paired() {
            return true;
        }
    }
    return false;
}

pub fn write_file(ofile: &str, filetype: &str, bg: Vec<(String, u64, u64, f64)>, chromsizes: HashMap<String, u32>) {
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
            false
        );
        // Theoretically one could add more threads here too, but this would require rewrite of the _bg iter upstream.
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()
            .expect("Unable to create tokio runtime for bw writing.");
        println!("Init writer");
        let writer = BigWigWrite::create_file(ofile, chromsizes).unwrap();
        println!("Start writer");
        let _ = writer.write(vals, runtime);
    }
}