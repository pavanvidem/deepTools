use rust_htslib::bam::{self, Read, IndexedReader, record::Cigar};
use std::collections::HashMap;
use std::fs::File;
use itertools::Itertools;
use tempfile::{Builder, NamedTempFile, TempPath};
use std::io::{BufReader, BufWriter, Write};

pub fn parse_regions(regions: &Vec<(String, u32, u32)>, bam_ifile: &str) -> (Vec<(String, u32, u32)>, HashMap<String, u32>) {
    // Takes a vector of regions, and a bam reference
    // returns a vector of regions, with all chromosomes and full lengths if original regions was empty
    // Else it validates the regions against the information from the bam header
    // Finally, a Vec with chromsizes is returned as well.

    let bam = IndexedReader::from_path(bam_ifile).unwrap();
    let header = bam.header().clone();
    let mut chromregions: Vec<(String, u32, u32)> = Vec::new();
    let mut chromsizes = HashMap::new();
    if regions.is_empty() {
        // if regions is empty, we default to all chromosomes, full length
        for tid in 0..header.target_count() {
            let chromname = String::from_utf8(header.tid2name(tid).to_vec())
                .expect("Invalid UTF-8 in chromosome name");
            let chromlen = header.target_len(tid)
                .expect("Error retrieving length for chromosome");
            chromregions.push((chromname.clone(), 0, chromlen as u32));
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
pub fn bam_pileup<'a>(
    bam_ifile: &str,
    region: &'a (String, u32, u32),
    binsize: &u32,
    ispe: &bool,
    ignorechr: &Vec<String>,
    filters: &alignmentfilters
) -> (
    Vec<TempPath>, // temp bedgraph file.
    u32, // mapped reads
    u32, // unmapped reads
    Vec<u32>, // read lengths
    Vec<u32>, // fragment lengths
)  {
    // constant to check if read is first in pair (relevant later)
    const FREAD: u16 = 0x40;

    // open bam file and fetch proper chrom
    let mut bam = IndexedReader::from_path(bam_ifile).unwrap();
    bam.fetch((region.0.as_str(), region.1, region.2))
        .expect(&format!("Error fetching region: {:?}", region));
    
    // init variables for mapping statistics and lengths
    let mut mapped_reads: u32 = 0;
    let mut unmapped_reads: u32 = 0;
    let mut readlens: Vec<u32> = Vec::new();
    let mut fraglens: Vec<u32> = Vec::new();
    
    // Create the output vector
    // let mut bg: Vec<(&str, u32, u32, f32)> = Vec::new();
    let bg = Builder::new()
        .prefix("deeptoolstmp_")
        .suffix(".bedgraph")
        .rand_bytes(12)
        .tempfile()
        .expect("Failed to create temporary file.");
    
    // Two cases: either the binsize is 1, or it is > 1.
    // Counting between the two modes is different. In binsize == 1 we compute pileups
    // for binsize > 1, we count the number of reads that overlap a bin.
    if binsize > &1 {
        // populate the bg vector with 0 counts over all bins
        let mut counts: Vec<f32> = vec![0.0; 1+((region.2 - region.1) / *binsize) as usize];
        let mut binstart = region.1;
        let mut binix: u32 = 0;
        
        for record in bam.records() {
            let record = record.expect("Error parsing record.");
            if !ignorechr.contains(&region.0) {
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

            let blocks = bam_blocks(record);
            // keep a record of bin indices that have been updated already
            let mut changedbins: Vec<u32> = Vec::new();
            for block in blocks.into_iter() {
                // Don't count blocks that exceed the chromosome
                if block.0 >= region.2 {
                    continue;
                }
                binix = block.0 / *binsize;
                if !changedbins.contains(&binix) {
                    counts[binix as usize] += 1.0;
                    changedbins.push(binix);
                }
                // Don't count blocks that exceed the chromosome
                if block.1 >= region.2 {
                    continue;
                }
                // if block.1 is at the end of a region, it should be counted in the block before (if different from first block)
                if block.1 == region.2 {
                    binix = (block.1 - 1) / *binsize;
                } else {
                    binix = block.1 / *binsize;
                }
                if !changedbins.contains(&binix) {
                    counts[binix as usize] += 1.0;
                    changedbins.push(binix);
                }
            }
            // if changedbins contains non-continuous blocks (perhaps read length spans multiple bins), update the inbetweens
            if let (Some(&min), Some(&max)) = (changedbins.iter().min(), changedbins.iter().max()) {
                for ix in min..=max {
                    if !changedbins.contains(&ix) {
                        counts[binix as usize] += 1.0;
                    }
                }
            }
        }
        // collaps bincounts into bedgraph format. 
        // merge bins with same coverage
        let mut lcov = counts[0];
        let mut lstart = region.1;
        let mut lend = region.1 + *binsize;
        let mut writer = BufWriter::new(&bg);
        for (ix, count) in counts.into_iter().skip(1).enumerate() {
            if count != lcov {
                //bg.push((&region.0, lstart, lend, lcov));
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, lstart, lend, lcov).unwrap();
                lstart = region.1 + (ix as u32 * binsize);
                lcov = count;
            }
            if region.1 + ((ix + 1) as u32 * binsize) > region.2 {
                lend = region.2;
            } else {
                lend = region.1 + ((ix + 1) as u32 * binsize);
            }
        }
        // write last entry
        writeln!(writer, "{}\t{}\t{}\t{}", region.0, lstart, lend, lcov).unwrap();
    } else {
        let mut l_start: u32 = region.1;
        let mut l_end: u32 = region.1;
        let mut l_cov: u32 = 0;
        let mut pileup_start: bool = true;
        let mut writer = BufWriter::new(&bg);
        let mut written: bool = false;
        for record in bam.records() {
            let record = record.expect("Error parsing record.");
            if !ignorechr.contains(&region.0) {
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
        }
        // Refetch bam file to reset iterator.
        bam.fetch((region.0.as_str(), region.1, region.2))
            .expect(&format!("Error fetching region: {:?}", region));
        for p in bam.pileup() {
            // Per default pileups count deletions in cigar string too.
            // For consistency with previous deepTools functionality, we ignore them.
            // to be fair I think they shouldn't be counted anyhow, but who am I ?
            // Note that coverages can be 0 now.
            let pileup = p.expect("Error parsing pileup.");
            let mut cov: u32 = 0;
            for _a in pileup.alignments() {
                if !_a.is_del() {
                    cov += 1;
                }
            }
            let pos = pileup.pos();
            if pileup_start {
                // if the first pileup is not at the start of the region, write 0 coverage
                if pos > l_start {
                    //bg.push((&region.0, l_start, pos, 0 as f32));
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, pos, 0 as f32).unwrap();
                    written = true;
                }
                pileup_start = false;
                l_start = pos;
                l_cov = cov;
            } else {
                if pos != l_end + 1 {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, l_end + 1, l_cov as f32).unwrap();
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_end + 1, pos, 0 as f32).unwrap();
                    written = true;
                    //bg.push((&region.0, l_start, l_end + 1, l_cov as f32));
                    //bg.push((&region.0, l_end + 1, pos, 0 as f32));
                    l_start = pos;
                    l_cov = cov;
                } else if l_cov != cov {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, pos, l_cov as f32).unwrap();
                    written = true;
                    // bg.push((&region.0, l_start, pos, l_cov as f32));
                    l_start = pos;
                }
            }
            l_end = pos;
            l_cov = cov;
            }
        // if bg is empty, whole region is 0 coverage
        if !written {
            writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, region.2, 0 as f32).unwrap();
            //bg.push((&region.0, l_start, region.2, 0 as f32));
        } else {
            // Still need to write the last pileup(s)
            writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, l_end + 1, l_cov as f32).unwrap();
            //bg.push((&region.0, l_start, l_end + 1, l_cov as f32));
            // Make sure that if we didn't reach end of chromosome, we still write 0 cov.
            if l_end + 1 < region.2 {
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_end + 1, region.2, 0 as f32).unwrap();
                //bg.push((&region.0, l_end + 1, region.2, 0 as f32));
            }
        }
    }
    let bgpath = bg.into_temp_path();
    let tmpvec   = vec![bgpath];
    return (tmpvec, mapped_reads, unmapped_reads, readlens, fraglens);
}

/// Extract contigious blocks from a bam record
/// Blocks are split per insertion, padding, deletion and ref skip
fn bam_blocks(rec: bam::Record) -> Vec<(u32, u32)> {
    let mut blocks: Vec<(u32, u32)> = Vec::new();
    let mut pos = rec.pos();

    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let start = pos;
                let end = pos + *len as i64;
                blocks.push((start as u32, end as u32));
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
pub fn collapse_bgvec(mut bg: Vec<(&str, u32, u32, f32)>, scale_factor: f32) -> Vec<(&str, u32, u32, f32)> {
    let mut cvec: Vec<(&str, u32, u32, f32)> = Vec::new();
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

// /// Takes two bedgraph vectors, and combines adjacent blocks with equal coverage
// #[allow(unused_assignments)]
// pub fn collapse_bgvecs(
//     mut bg1: Vec<(&str, u32, u32, f32)>,
//     mut bg2: Vec<(&str, u32, u32, f32)>,
//     sf1: f32,
//     sf2: f32,
//     pseudocount: f32,
//     operation: &str
// ) -> Vec<(&str, u32, u32, f32)> {
//     assert_eq!(bg1.len(), bg2.len(), "Bedgraph vectors must be of equal length.");
//     let mut cvec: Vec<(&str, u32, u32, f32)> = Vec::new();
//     // initialize values, assert equal bins
//     let (mut lchrom, mut lstart, mut lend, lcov1) = bg1.remove(0);
//     let (lchrom2, lstart2, lend2, lcov2) = bg2.remove(0);
//     assert_eq!(lchrom, lchrom2,"Chromosomes in bedgraph vectors must be equal.");
//     assert_eq!(lstart, lstart2,"Start positions in bedgraph vectors must be equal.");
//     assert_eq!(lend, lend2,"End positions in bedgraph vectors must be equal.");

//     let mut lcov: f32 = calc_ratio(lcov1, lcov2, &sf1, &sf2, &pseudocount, &operation);
//     for (
//         (chrom, start, end, cov1),
//         (chrom2, start2, end2, cov2)
//     ) in bg1.into_iter().zip(bg2.into_iter()) {
//         // assert equal bins.
//         assert_eq!(chrom, chrom2,"Chromosomes in bedgraph vectors must be equal.");
//         assert_eq!(start, start2,"Start positions in bedgraph vectors must be equal.");
//         assert_eq!(end, end2,"End positions in bedgraph vectors must be equal.");
//         // Calculate coverage, depending on what operation is fed
//         let cov: f32 = calc_ratio(cov1, cov2, &sf1, &sf2, &pseudocount, &operation);
//         // Collapse bg vector properly
//         if chrom != lchrom {
//             cvec.push((lchrom, lstart, lend, lcov));
//             lchrom = chrom;
//             lstart = start;
//             lend = end;
//             lcov = cov;
//         } else if cov != lcov {
//             cvec.push((lchrom, lstart, lend, lcov));
//             lchrom = chrom;
//             lstart = start;
//             lend = end;
//             lcov = cov;
//         }
//         lend = end;
//     }
//     cvec.push((lchrom, lstart, lend, lcov));
//     return cvec;
// }

fn calc_ratio(
    cov1: f32,
    cov2: f32,
    sf1: &f32,
    sf2: &f32,
    pseudocount: &f32,
    operation: &str
) -> f32 {
    // Pseudocounts are only used in log2 and ratio operations
    // First scale factor is applied, then pseudocount, if applicable.
    match operation {
        "log2" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
        "ratio" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return num / den;
        }
        "reciprocal_ratio" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            let ratio: f32 = num / den;
            if ratio >= 1.0 {
                return den / num;
            } else {
                return -num / den;
            }
        }
        _ => {
            // No operation is never allowed (on the py arg level, so just default to log2)
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
    }
}
pub struct alignmentfilters {
    pub minmappingquality: u8,
    pub samflaginclude: u16,
    pub samflagexclude: u16,
    pub minfraglen: u32,
    pub maxfraglen: u32,
}