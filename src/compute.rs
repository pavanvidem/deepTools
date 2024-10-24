use rust_htslib::bam::{Read, IndexedReader};

pub fn parse_regions(regions: &Vec<(String, u64, u64)>, bam_ifile: &str) -> Vec<(String, u64, u64)> {
    // Takes a vector of regions, and a bam reference
    // returns a vector of regions, with all chromosomes and full lengths if original regions was empty

    let bam = IndexedReader::from_path(bam_ifile).unwrap();
    let header = bam.header().clone();

    let mut chromregions: Vec<(String, u64, u64)> = Vec::new();
    if regions.is_empty() {
        // if regions is empty, we default to all chromosomes, full length
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header.target_len(tid)
                .expect("Error retrieving length for chromosome");
            chromregions.push((chromname, 0, chromlen));
        }
    } else {
        let validchroms: Vec<String> = header
            .target_names()
            .iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap())
            .collect();
        
        for region in regions {
            let chromname = &region.0;
            assert!(region.1 < region.2, "Region start must be strictly less than region end.");
            assert!(validchroms.contains(chromname), "Chromosome {} not found in bam header", chromname);
            chromregions.push((chromname.clone(), region.1, region.2));
        }
    }
    return chromregions;
}

// unused as values are only used in the next iteration of the loop.
#[allow(unused_assignments)]
pub fn bam_pileup(bam_ifile: &str, region: &(String, u64, u64)) -> Vec<(String, u64, u64, u64)> {

    let mut bam = IndexedReader::from_path(bam_ifile).unwrap();
    bam.fetch((region.0.as_str(), region.1, region.2))
        .expect(&format!("Error fetching region: {:?}", region));

    let mut bg: Vec<(String, u64, u64, u64)> = Vec::new();
    let mut l_start: u64 = region.1;
    let mut l_end: u64 = region.1;
    let mut l_cov: u64 = 0;
    // let chrlen: u64 = bam.header().target_len(bam.header().tid(region.0.as_bytes()).unwrap()).unwrap();
    let mut pileup_start: bool = true;

    for p in bam.pileup() {
        // Per default pileups count deletions in cigar string too.
        // For consistency with previous deepTools functionality, we ignore them.
        // to be fair I think they shouldn't be counted anyhow, but who am I ?
        // Note that coverages can be 0 now.
        let pileup = p.expect("Error parsing pileup.");
        let mut cov: u64 = 0;
        for _a in pileup.alignments() {
            if !_a.is_del() {
                cov += 1;
            }
        }
        let pos = pileup.pos() as u64;
        if pileup_start {
            // if the first pileup is not at the start of the region, write 0 coverage
            if pos > l_start {
                bg.push((region.0.clone(), l_start, pos, 0));
            }
            pileup_start = false;
            l_start = pos;
            l_cov = cov;
        } else {
            if pos != l_end + 1 {
                bg.push((region.0.clone(), l_start, l_end + 1, l_cov));
                bg.push((region.0.clone(), l_end + 1, pos, 0));
                l_start = pos;
                l_cov = cov;
            } else if l_cov != cov {
                bg.push((region.0.clone(), l_start, pos, l_cov));
                l_start = pos;
            }
        }
        l_end = pos;
        l_cov = cov;
        }
    // if bg is empty, whole region is 0 coverage
    println!("Reached end of chrom {}, l_start = {}, l_end = {}, cov = {}, l_cov = {}, regionend = {}", region.0, l_start, l_end, l_cov, l_cov, region.2);
    if bg.is_empty() {
        bg.push((region.0.clone(), l_start, region.2, 0));
    } else {
        // Still need to write the last pileup(s)
        bg.push((region.0.clone(), l_start, l_end + 1, l_cov));
        // Make sure that if we didn't reach end of chromosome, we still write 0 cov.
        if l_end + 1 < region.2 {
            bg.push((region.0.clone(), l_end + 1, region.2, 0));
        }
    }
    return bg;
}