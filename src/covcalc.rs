use rust_htslib::bam::{self, Read, IndexedReader, record::Cigar};
use std::collections::HashMap;
use tempfile::{Builder, TempPath};
use std::io::{BufWriter, Write};
use std::cmp::min;

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
    _filters: &Alignmentfilters,
    collapse: bool,
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
        let mut counts: Vec<f32> = vec![0.0; (region.2 - region.1).div_ceil(*binsize) as usize];
        // let mut binstart = region.1;
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
                binix = block.0 / binsize;
                if !changedbins.contains(&binix) {
                    counts[binix as usize] += 1.0;
                    changedbins.push(binix);
                }
                // Don't count blocks that exceed the chromosome
                if block.1 >= region.2 {
                    continue;
                }
                // Note that our blocks are strictly less then region ends, hence no problem with bin spewing at end:
                // binix = [a,b) where b == region.2 would divide into binix+1 (doesn't exist).
                binix = block.1 / binsize;
                if !changedbins.contains(&binix) {
                    counts[binix as usize] += 1.0;
                    changedbins.push(binix);
                }
            }
            // if changedbins contains non-continuous blocks (perhaps read length spans multiple bins), update the inbetweens
            if let (Some(&min), Some(&max)) = (changedbins.iter().min(), changedbins.iter().max()) {
                for ix in min..=max {
                    if !changedbins.contains(&ix) {
                        counts[ix as usize] += 1.0;
                        changedbins.push(ix);
                    }
                }
            }
        }
        let mut writer = BufWriter::new(&bg);
        // There are two scenarios: 
        // bamCoverage mode -> we can collapse bins with same coverage (collapse = true)
        // bamCompare & others -> We cannot collapse the bins, yet. (collapse = false)
        if counts.len() == 1 {
            writeln!(writer, "{}\t{}\t{}\t{}", region.0, region.1, region.2, counts[0]).unwrap();
        } else {
            if collapse {
                let mut lcov = counts[0];
                let mut lstart = region.1;
                let mut lend = region.1 + binsize;
                let mut start = lstart;
                let mut end = lend;
                let mut bin: u32 = 0;
    
                for (ix, count) in counts.into_iter().skip(1).enumerate() {
                    bin = (ix + 1) as u32; // offset of 1 due to skip(1)
                    start = (bin * binsize) + region.1;
                    end = min(start + binsize, region.2);
                    if count != lcov {
                        //bg.push((&region.0, lstart, lend, lcov));
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, lstart, lend, lcov).unwrap();
                        lstart = lend;
                        lcov = count;
                    }
                    lend = end;
                }
                // write last entry
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, lstart, lend, lcov).unwrap();
            } else {
                let mut start = region.1;
                let mut end = region.1 + binsize;
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, start, end, counts[0]).unwrap();
                for (ix, count) in counts.into_iter().skip(1).enumerate() {
                    let bin = (ix + 1) as u32;
                    start = (bin * binsize) + region.1;
                    end = min(start + binsize, region.2);
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, start, end, count).unwrap();
                }
            }
        }

    // binsize = 1
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
                    if collapse {
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, pos, 0 as f32).unwrap();
                    } else {
                        for p in l_start..pos {
                            writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, 0 as f32).unwrap();
                        }
                    }
                    written = true;
                }
                pileup_start = false;
                l_start = pos;
                l_cov = cov;
            } else {
                if pos != l_end + 1 {
                    if collapse {
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, l_end + 1, l_cov as f32).unwrap();
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_end + 1, pos, 0 as f32).unwrap();
                    } else {
                        for p in l_start..l_end + 1 {
                            writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, l_cov as f32).unwrap();
                        }
                        for p in l_end + 1..pos {
                            writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, 0 as f32).unwrap();
                        }
                    }
                    written = true;
                    l_start = pos;
                    l_cov = cov;
                } else if l_cov != cov {
                    if collapse {
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, pos, l_cov as f32).unwrap();
                    } else {
                        for p in l_start..pos {
                            writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, l_cov as f32).unwrap();
                        }
                    }
                    written = true;
                    l_start = pos;
                    l_cov = cov;
                }
            }
            l_end = pos;
            l_cov = cov;
        }
        // if bg is empty, whole region is 0 coverage
        if !written {
            if collapse {
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, region.2, 0 as f32).unwrap();
            } else {
                for p in l_start..region.2 {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, 0 as f32).unwrap();
                }
            }
        } else {
            // Still need to write the last pileup(s)
            if collapse {
                writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_start, l_end + 1, l_cov as f32).unwrap();
                // Make sure that if we didn't reach end of chromosome, we still write 0 cov.
                if l_end + 1 < region.2 {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, l_end + 1, region.2, 0 as f32).unwrap();
                }
            } else {
                for p in l_start..l_end + 1 {
                    writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, l_cov as f32).unwrap();
                }
                if l_end + 1 < region.2 {
                    for p in l_end + 1..region.2 + 1 {
                        writeln!(writer, "{}\t{}\t{}\t{}", region.0, p, p + 1, 0 as f32).unwrap();
                    }
                }
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
                let end = pos + *len as i64 -1;
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

pub struct Alignmentfilters {
    pub minmappingquality: u8,
    pub samflaginclude: u16,
    pub samflagexclude: u16,
    pub minfraglen: u32,
    pub maxfraglen: u32,
}