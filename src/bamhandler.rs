use rust_htslib::bam::{Read, Reader};
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