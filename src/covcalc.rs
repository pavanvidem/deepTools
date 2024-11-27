use rust_htslib::bam::{self, Read, IndexedReader, record::Cigar};
use std::collections::HashMap;
use itertools::Itertools;

pub fn parse_regions(regions: &Vec<(String, u64, u64)>, bam_ifile: &str) -> (Vec<(String, u64, u64)>, HashMap<String, u32>) {
    // Takes a vector of regions, and a bam reference
    // returns a vector of regions, with all chromosomes and full lengths if original regions was empty
    // Else it validates the regions against the information from the bam header
    // Finally, a Vec with chromsizes is returned as well.

    let bam = IndexedReader::from_path(bam_ifile).unwrap();
    let header = bam.header().clone();
    let mut chromregions: Vec<(String, u64, u64)> = Vec::new();
    let mut chromsizes = HashMap::new();
    if regions.is_empty() {
        // if regions is empty, we default to all chromosomes, full length
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header.target_len(tid)
                .expect("Error retrieving length for chromosome");
            chromregions.push((chromname.clone(), 0, chromlen));
            chromsizes.insert(chromname.to_string(), chromlen as u32);
        }
    } else {
        // populate chromsizes
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header.target_len(tid)
                .expect("Error retrieving length for chromosome");
            chromsizes.insert(chromname, chromlen as u32);
        }
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
    // Sort regions to make our live easier down the line
    chromregions.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    return (chromregions, chromsizes);
}

/// Main workhorse for bamCoverage and bamCompare
/// Calculates coverage either per bp (bs = 1) or over bins (bs > 1)
#[allow(unused_assignments)]
pub fn bam_pileup(
    bam_ifile: &str,
    region: &(String, u64, u64),
    binsize: &u32,
    ispe: &bool
) -> (
    Vec<(String, u64, u64, f64)>, // bedgraph vec
    u64, // mapped reads
    u64, // unmapped reads
    Vec<u32>, // read lengths
    Vec<u32>, // fragment lengths
)  {
    
    // constant to check if read is first in pair (if relevant later)
    const FREAD: u16 = 0x40;

    // open bam file and fetch proper chrom
    let mut bam = IndexedReader::from_path(bam_ifile).unwrap();
    bam.fetch((region.0.as_str(), region.1, region.2))
        .expect(&format!("Error fetching region: {:?}", region));
    
    // init variables for mapping statistics and lengths
    let mut mapped_reads: u64 = 0;
    let mut unmapped_reads: u64 = 0;
    let mut readlens: Vec<u32> = Vec::new();
    let mut fraglens: Vec<u32> = Vec::new();
    
    // Create the output vector
    let mut bg: Vec<(String, u64, u64, f64)> = Vec::new();

    // Two cases: either the binsize is 1, or it is > 1.
    if binsize > &1 {
        // Construct a hashmap of bins with their counts
        let mut bin_counts: HashMap<u64, (u64, u64, u64)> = HashMap::new();
        
        // populate hashmap
        let mut binstart = region.1;
        let mut binix: u64 = 0;
        while binstart < region.2 {
            let mut binend = binstart + *binsize as u64;
            if binend > region.2 {
                binend = region.2;
            }
            bin_counts.insert(binix, (binstart, binend, 0));
            binstart = binend;
            binix += 1;
        }

        for record in bam.records() {
            let record = record.expect("Error parsing record.");
            if record.is_unmapped() {
                unmapped_reads += 1;
            } else {
                mapped_reads += 1;
                if *ispe {
                    if record.is_paired() && record.is_proper_pair() && (record.flags() & FREAD != 0) {
                        fraglens.push(record.insert_size().abs() as u32);
                    }
                }
                readlens.push(record.seq_len() as u32);
            }

            let blocks = bam_blocks(record);
            // keep a record of bin indices that have been updated already
            let mut changedbins: Vec<u64> = Vec::new();
            // split off first entry
            for block in blocks.into_iter() {
                // Don't count blocks that exceed the chromosome
                if block.0 as u64 >= region.2 {
                    continue;
                }
                binix = block.0 as u64 / *binsize as u64;
                if !changedbins.contains(&binix) {
                    bin_counts.get_mut(&binix).expect(
                        &format!("Bin index {} not in hashmap. Blocks {:?}. Region = {}. Bamfile = {}", &binix, block, region.0, bam_ifile)
                    ).2 += 1;
                    changedbins.push(binix);
                }
                // Don't count blocks that exceed the chromosome
                if block.1 as u64 >= region.2 {
                    continue;
                }
                // if block.1 is at the end of a region, it should be counted in the block before (if different from first block)
                if block.1 as u64 == region.2 {
                    binix = (block.1 as u64 - 1) / *binsize as u64;
                } else {
                    binix = block.1 as u64 / *binsize as u64;
                }
                if !changedbins.contains(&binix) {
                    bin_counts.get_mut(&binix).expect(
                        &format!("Bin index {} not in hashmap. Blocks {:?}. Region = {}. Bamfile = {}", &binix, block, region.0, bam_ifile)
                    ).2 += 1;
                    changedbins.push(binix);
                }
            }
            // if changedbins contains non-continuous blocks (perhaps read length spans multiple bins), update the inbetweens
            if let (Some(&min), Some(&max)) = (changedbins.iter().min(), changedbins.iter().max()) {
                for ix in min..=max {
                    if !changedbins.contains(&ix) {
                        bin_counts.get_mut(&ix).unwrap().2 += 1;
                    }
                }
            }
        }
        bg = bin_counts
            .iter()
            .sorted_by_key(|(&k, _)| k)
            .map(|(_k, &(binstart, binend, count))| (region.0.clone(), binstart, binend, count as f64))
            .collect();

    } else {
        let mut l_start: u64 = region.1;
        let mut l_end: u64 = region.1;
        let mut l_cov: u64 = 0;
        let mut pileup_start: bool = true;

        for record in bam.records() {
            let record = record.expect("Error parsing record.");
            if record.is_unmapped() {
                unmapped_reads += 1;
            } else {
                mapped_reads += 1;
                if *ispe {
                    if record.is_paired() && record.is_proper_pair() && (record.flags() & FREAD != 0) {
                        fraglens.push(record.insert_size().abs() as u32);
                    }
                }
                readlens.push(record.seq_len() as u32);
            }
        }
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
                    bg.push((region.0.clone(), l_start, pos, 0 as f64));
                }
                pileup_start = false;
                l_start = pos;
                l_cov = cov;
            } else {
                if pos != l_end + 1 {
                    bg.push((region.0.clone(), l_start, l_end + 1, l_cov as f64));
                    bg.push((region.0.clone(), l_end + 1, pos, 0 as f64));
                    l_start = pos;
                    l_cov = cov;
                } else if l_cov != cov {
                    bg.push((region.0.clone(), l_start, pos, l_cov as f64));
                    l_start = pos;
                }
            }
            l_end = pos;
            l_cov = cov;
            }
        // if bg is empty, whole region is 0 coverage
        if bg.is_empty() {
            bg.push((region.0.clone(), l_start, region.2, 0 as f64));
        } else {
            // Still need to write the last pileup(s)
            bg.push((region.0.clone(), l_start, l_end + 1, l_cov as f64));
            // Make sure that if we didn't reach end of chromosome, we still write 0 cov.
            if l_end + 1 < region.2 {
                bg.push((region.0.clone(), l_end + 1, region.2, 0 as f64));
            }
        }
    }
    // Collect median read lengths and fragment lengths if needed
    return (bg, mapped_reads, unmapped_reads, readlens, fraglens);
}

/// Extract contigious blocks from a bam record
/// Blocks are split per insertion, padding, deletion and ref skip
fn bam_blocks(rec: bam::Record) -> Vec<(i64, i64)> {
    let mut blocks: Vec<(i64, i64)> = Vec::new();
    let mut pos = rec.pos();

    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let start = pos;
                let end = pos + *len as i64;
                blocks.push((start, end));
                pos = end;
            }
            Cigar::Ins(len) | Cigar::Pad(len) => {
                pos += *len as i64;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                pos += *len as i64;
            }
            _ => (),
        }
    }
    return blocks;
}

/// Takes a bedgraph vector, and combines adjacent blocks with equal coverage
#[allow(unused_assignments)]
pub fn collapse_bgvec(mut bg: Vec<(String, u64, u64, f64)>, scale_factor: f64) -> Vec<(String, u64, u64, f64)> {
    let mut cvec: Vec<(String, u64, u64, f64)> = Vec::new();
    // initialize values
    let (mut lchrom, mut lstart, mut lend, mut lcov) = bg.remove(0);
    for (chrom, start, end, cov) in bg.into_iter() {
        if chrom != lchrom {
            cvec.push((lchrom, lstart, lend, lcov * scale_factor));
            lchrom = chrom;
            lstart = start;
            lend = end;
            lcov = cov;
        } else if cov != lcov {
            cvec.push((lchrom, lstart, lend, lcov * scale_factor));
            lchrom = chrom;
            lstart = start;
            lend = end;
            lcov = cov;
        }
        lend = end;
    }
    cvec.push((lchrom, lstart, lend, lcov * scale_factor));
    return cvec;
}

/// Takes two bedgraph vectors, and combines adjacent blocks with equal coverage
#[allow(unused_assignments)]
pub fn collapse_bgvecs(
    mut bg1: Vec<(String, u64, u64, f64)>,
    mut bg2: Vec<(String, u64, u64, f64)>,
    sf1: f64,
    sf2: f64,
    pseudocount: f64,
    operation: &str
) -> Vec<(String, u64, u64, f64)> {
    assert_eq!(bg1.len(), bg2.len(), "Bedgraph vectors must be of equal length.");
    let mut cvec: Vec<(String, u64, u64, f64)> = Vec::new();
    // initialize values, assert equal bins
    let (mut lchrom, mut lstart, mut lend, lcov1) = bg1.remove(0);
    let (lchrom2, lstart2, lend2, lcov2) = bg2.remove(0);
    assert_eq!(lchrom, lchrom2,"Chromosomes in bedgraph vectors must be equal.");
    assert_eq!(lstart, lstart2,"Start positions in bedgraph vectors must be equal.");
    assert_eq!(lend, lend2,"End positions in bedgraph vectors must be equal.");

    let mut lcov: f64 = calc_ratio(lcov1, lcov2, &sf1, &sf2, &pseudocount, &operation);
    for (
        (chrom, start, end, cov1),
        (chrom2, start2, end2, cov2)
    ) in bg1.into_iter().zip(bg2.into_iter()) {
        // assert equal bins.
        assert_eq!(chrom, chrom2,"Chromosomes in bedgraph vectors must be equal.");
        assert_eq!(start, start2,"Start positions in bedgraph vectors must be equal.");
        assert_eq!(end, end2,"End positions in bedgraph vectors must be equal.");
        // Calculate coverage, depending on what operation is fed
        let cov: f64 = calc_ratio(cov1, cov2, &sf1, &sf2, &pseudocount, &operation);
        // Collapse bg vector properly
        if chrom != lchrom {
            cvec.push((lchrom, lstart, lend, lcov));
            lchrom = chrom;
            lstart = start;
            lend = end;
            lcov = cov;
        } else if cov != lcov {
            cvec.push((lchrom, lstart, lend, lcov));
            lchrom = chrom;
            lstart = start;
            lend = end;
            lcov = cov;
        }
        lend = end;
    }
    cvec.push((lchrom, lstart, lend, lcov));
    return cvec;
}

fn calc_ratio(
    cov1: f64,
    cov2: f64,
    sf1: &f64,
    sf2: &f64,
    pseudocount: &f64,
    operation: &str
) -> f64 {
    // Pseudocounts are only used in log2 and ratio operations
    // First scale factor is applied, then pseudocount, if applicable.
    match operation {
        "log2" => {
            let num: f64 = (cov1 * *sf1) + *pseudocount;
            let den: f64 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
        "ratio" => {
            let num: f64 = (cov1 * *sf1) + *pseudocount;
            let den: f64 = (cov2 * *sf2) + *pseudocount;
            return num / den;
        }
        "reciprocal_ratio" => {
            let num: f64 = (cov1 * *sf1) + *pseudocount;
            let den: f64 = (cov2 * *sf2) + *pseudocount;
            let ratio: f64 = num / den;
            if ratio >= 1.0 {
                return den / num;
            } else {
                return -num / den;
            }
        }
        _ => {
            // No operation is never allowed (on the py arg level, so just default to log2)
            let num: f64 = (cov1 * *sf1) + *pseudocount;
            let den: f64 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
    }
}