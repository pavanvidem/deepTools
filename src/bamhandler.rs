use rust_htslib::bam::{Read, Reader};

pub fn bam_stats(bam_ifile: &str) -> (u64, u64, u64, u32) {
    let mut fraglens: Vec<u32> = Vec::new();
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
        }
    }
    return (total_reads, mapped_reads, unmapped_reads)
}