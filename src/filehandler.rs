use bigtools::utils::misc::Name;
use rust_htslib::bam::{Read, Reader};
use itertools::Itertools;
use std::io::{BufReader, BufWriter, Write};
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use bigtools::{BigWigRead, BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::bed::bedparser::parse_bedgraph;
use std::collections::HashMap;
use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::{TempPath};
use crate::computematrix::{Scalingregions, Revalue, Region, Bin};
use crate::calc::{mean_float, median_float, min_float, max_float, sum_float, std_float};

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

pub fn write_covfile<LI>(lines: LI, ofile: &str, filetype: &str, chromsizes: HashMap<String, u32>)
where
    LI: Iterator<Item = (String, Value)>,
 {
    if filetype == "bedgraph" {
        // write output file, bedgraph
        let mut writer = BufWriter::new(File::create(ofile).unwrap());
        for (chrom, val) in lines {
            writeln!(writer, "{}\t{}\t{}\t{}", chrom, val.start, val.end, val.value).unwrap();
        }
    } else {
        let vals = BedParserStreamingIterator::wrap_infallible_iter(
            lines,
            false
        );
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()
            .expect("Unable to create tokio runtime for bw writing.");
        let writer = BigWigWrite::create_file(ofile, chromsizes).unwrap();
        let _ = writer.write(vals, runtime);
    }
}

pub fn read_bedfiles(bed_files: &Vec<String>, metagene: bool) -> (Vec<Region>, HashMap<String, u32>) {
    // read all bedfiles in a Vec of strings (filepaths)
    // returns a vec of Region
    let mut regionsizes: HashMap<String, u32> = HashMap::new();
    let mut regions: Vec<Region> = Vec::new();

    //let mut regions: Vec<(String, u32, u32, String, String, String)> = Vec::new();
    let mut nonbed12: bool = false;

    for bed in bed_files {
        let entryname = Path::new(bed)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        
        let bedfile = BufReader::new(File::open(bed).unwrap());
        let mut entries: u32 = 0;
        for line in bedfile.lines() {
            let line = line.unwrap();
            let fields: Vec<&str> = line.split('\t').collect();
            // Depending on bedfile, we have either BED3, BED6 or BED12
            // Note that this approach could allow somebody to have a 'mixed' bedfile, why not.
            match fields.len() {
                3 => {
                    if !nonbed12 {
                        nonbed12 = true;
                    }
                    let start = fields[1].parse().unwrap();
                    let end = fields[2].parse().unwrap();
                    regions.push(
                        Region {
                            chrom: fields[0].to_string(), //chrom
                            start: Revalue::U(start), //start
                            end: Revalue::U(end), //end
                            score: ".".to_string(), //score
                            strand: ".".to_string(), //score
                            name: entryname.to_string(), //bedfile_name
                            regionlength: end - start // regionlength
                        }
                    );
                },
                6 => {
                    if !nonbed12 {
                        nonbed12 = true;
                    }
                    let start = fields[1].parse().unwrap();
                    let end = fields[2].parse().unwrap();
                    regions.push(
                        Region {
                            chrom: fields[0].to_string(), //chrom
                            start: Revalue::U(start), //start
                            end: Revalue::U(end), //end
                            score: ".".to_string(), //score
                            strand: ".".to_string(), //score
                            name: entryname.to_string(), //bedfile_name
                            regionlength: end - start // regionlength
                        }
                    );
                },
                12 => {
                    if metagene {                        
                        let start: u32 = fields[1].parse().unwrap();
                        let blocksizes: Vec<u32> = fields[10]
                            .split(',')
                            .filter(|x| !x.is_empty())
                            .map(|x| x.parse().unwrap())
                            .collect();
                        let length: u32 = blocksizes.iter().sum();
                        let blockstarts: Vec<u32> = fields[11]
                            .split(',')
                            .filter(|x| !x.is_empty())
                            .map(|x| x.parse::<u32>().unwrap() + start)
                            .collect();

                        let (starts, ends) = blocksizes
                            .into_iter()
                            .zip(blockstarts.into_iter())
                            .map(|(s, start)| (start, start + s))
                            .into_iter()
                            .unzip();
                        regions.push(
                            Region {
                                chrom: fields[0].to_string(), //chrom
                                start: Revalue::V(starts), //start
                                end: Revalue::V(ends), //end
                                score: fields[4].to_string(), //score
                                strand: fields[5].to_string(), //score
                                name: entryname.to_string(), //bedfile_name
                                regionlength: length // regionlength
                            }
                        );
                    } else {
                        let start = fields[1].parse().unwrap();
                        let end = fields[2].parse().unwrap();
                        regions.push(
                            Region {
                                chrom: fields[0].to_string(), //chrom
                                start: Revalue::U(start), //start
                                end: Revalue::U(end), //end
                                score: fields[4].to_string(), //score
                                strand: fields[5].to_string(), //score
                                name: entryname.to_string(), //bedfile_name
                                regionlength: end - start // regionlength
                            }
                        );
                    }
                },
                _ => panic!("Invalid BED format. BED file doesn't have 3, 6 or 12 fields."),
            }
            entries += 1;
        }
        regionsizes.insert(entryname, entries);
    }
    if metagene && nonbed12 {
        println!("Warning: Metagene analysis is requested, but not all bedfiles and/or bedfile entries are in BED12 format. Proceed at your own risk.");
    }
    return (regions, regionsizes);
}

pub fn chrombounds_from_bw(bwfile: &str) -> HashMap<String, u32> {
    // define chromsizes hashmap
    let mut chromsizes: HashMap<String, u32> = HashMap::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let mut reader = BigWigRead::open(bwf).unwrap();
    for chrom in reader.chroms() {
        chromsizes.insert(chrom.name.clone(), chrom.length);
    }
    chromsizes
}

pub fn bwintervals(
    bwfile: &str,
    regions: &Vec<Region>,
    slopregions: &Vec<Vec<Bin>>,
    scale_regions: &Scalingregions
) -> Vec<Vec<f32>> {
    // For a given bw file, a vector of slopregions (Bin enum))
    // return a vector with for every region a vector of f64.

    // Make sure regions and slopregions are of equal length
    assert_eq!(
        regions.len(), slopregions.len(),
        "Regions from bed file and parsed regions (slopped) do not have equal length. Something went wrong during computation."
    );
    
    // Define return vector, set up bw reader.
    let mut bwvals: Vec<Vec<f32>> = Vec::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let mut reader = BigWigRead::open(bwf).unwrap();

    // Iterate over regions and slopregions synchronously.
    for (sls, region) in slopregions.iter().zip(regions.iter()) {
        let mut bwval: Vec<f32> = Vec::new();
        // get 'min and max' to query
        // at some point this should become an impl but man am I tired.
        let (min, max) = sls
            .iter()
            .flat_map(|bin| match bin {
                Bin::Conbin(a, b) => vec![*a, *b],
                Bin::Catbin(pairs) => pairs.iter().flat_map(|(x, y)| vec![*x, *y]).collect::<Vec<u32>>(),
            })
            .fold((u32::MAX, u32::MIN), |(min, max), x| (min.min(x), max.max(x)));
        let binvals = reader.get_interval(&region.chrom, min as u32, max as u32).unwrap();
        // since binvals (can) be over binsizes, we expand them to bp and push them to a hashmap
        let mut bwhash: HashMap<u32, f32> = HashMap::new();
        for interval in binvals {
            let interval = interval.unwrap();
            let start = interval.start as u32;
            let end = interval.end as u32;
            let val = interval.value as f32;
            bwhash.extend((start..end).map(|bp| (bp, val)));
        }
        // Now we can iterate over the slopped regions, and get the values from the hashmap.
        for bin in sls {
            match bin {
                Bin::Conbin(a, b) => {
                    if a == b && *b == 0 {
                        if scale_regions.missingdata_as_zero {
                            bwval.push(0.0);
                        } else {
                            bwval.push(std::f32::NAN);
                        }
                    } else {
                        // Get values from the hashmap
                        let vals: Vec<&f32> = (*a..*b)
                            .filter_map(|bp| bwhash.get(&bp))
                            .collect();
                        let val = match scale_regions.avgtype.as_str() {
                            "mean" => mean_float(vals),
                            "median" => median_float(vals),
                            "min" => min_float(vals),
                            "max" => max_float(vals),
                            "std" => std_float(vals),
                            "sum" => sum_float(vals),
                            _ => panic!("Unknown avgtype."),
                        };
                        bwval.push(val);
                    }
                },
                Bin::Catbin(pairs) => {
                    let mut vals: Vec<&f32> = Vec::new();

                    for (start, end) in pairs {
                        if start == end && *end == 0 {
                            if scale_regions.missingdata_as_zero {
                                vals.push(&0.0);
                            } else {
                                vals.push(&std::f32::NAN);
                            }
                        } else {
                            // Get values from the hashmap
                            (*start..*end)
                                .filter_map(|bp| bwhash.get(&bp))
                                .for_each(|v| vals.push(v));
                        }
                    }

                    let val = match scale_regions.avgtype.as_str() {
                        "mean" => mean_float(vals),
                        "median" => median_float(vals),
                        "min" => min_float(vals),
                        "max" => max_float(vals),
                        "std" => std_float(vals),
                        "sum" => sum_float(vals),
                        _ => panic!("Unknown avgtype."),
                    };
                    bwval.push(val);
                }
            }
        }
        assert_eq!(bwval.len(), scale_regions.cols_expected / scale_regions.bwfiles);
        bwvals.push(bwval);
    }
    bwvals
}

pub fn header_matrix(scale_regions: &Scalingregions, regionsizes: HashMap<String, u32>) -> String {
    // Create the header for the matrix.
    // This is quite ugly, but this is mainly because we need to accomodate delta bwfiles.
    let mut headstr = String::new();
    headstr.push_str("@{");
    headstr.push_str(
        &format!("\"upstream\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.upstream).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"downstream\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.downstream).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"body\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.regionbodylength).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"bin size\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.binsize).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"ref point\":[\"{}\"],", (0..scale_regions.bwfiles).map(|_| scale_regions.referencepoint.clone()).collect::<Vec<_>>().into_iter().join("\",\""))
    );
    headstr.push_str(
        &format!("\"verbose\":{},", scale_regions.verbose)
    );
    headstr.push_str(
        &format!("\"bin avg type\":\"{}\",", scale_regions.avgtype)
    );
    headstr.push_str(
        &format!("\"missing data as zero\":{},", scale_regions.missingdata_as_zero)
    );
    // Unimplemented arguments, but they need to be present in header for now anyway.
    headstr.push_str(
        "\"min threshold\":null,\"max threshold\":null,\"scale\":1,\"skip zeros\":false,\"nan after end\":false,"
    );
    headstr.push_str(
        &format!("\"proc number\":{},", scale_regions.proc_number)
    );
    // Unimplemented for now...
    headstr.push_str(
        "\"sort regions\":\"keep\",\"sort using\":\"mean\","
    );
    headstr.push_str(
        &format!("\"unscaled 5 prime\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.unscaled5prime).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"unscaled 3 prime\":[{}],", (0..scale_regions.bwfiles).map(|_| scale_regions.unscaled3prime).collect::<Vec<_>>().into_iter().join(","))
    );
    headstr.push_str(
        &format!("\"group_labels\":[\"{}\"],", scale_regions.bedlabels.join("\",\""))
    );
    // Get cumulative sizes of regions
    let mut groupbounds: Vec<u32> = Vec::new();
    groupbounds.push(0);
    let mut cumsum: u32 = 0;
    for bedlabel in scale_regions.bedlabels.iter() {
        cumsum += regionsizes.get(bedlabel).unwrap();
        groupbounds.push(cumsum);
    }
    let groupbounds = format!("{}", groupbounds.iter()
        .map(|&x| x.to_string())
        .collect::<Vec<String>>()
        .join(","));

    headstr.push_str(
        &format!("\"group_boundaries\":[{}],", groupbounds)
    );
    // Sample labels
    headstr.push_str(
        &format!("\"sample_labels\":[\"{}\"],", scale_regions.bwlabels.join("\",\""))
    );
    // Get cumulative sizes for sample boundaries
    let colsize_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let mut samplebounds: Vec<u64> = Vec::new();
    samplebounds.push(0);
    let mut cumsum: u64 = 0;
    for i in 0..scale_regions.bwfiles {
        cumsum += colsize_per_sample as u64;
        samplebounds.push(cumsum);
    }
    let samplebounds = format!("{}", samplebounds.iter()
        .map(|&x| x.to_string())
        .collect::<Vec<String>>()
        .join(","));

    headstr.push_str(
        &format!("\"sample_boundaries\":[{}]", samplebounds)
    );
    headstr.push_str(
        "}\n"
    );
    headstr
}

pub fn write_matrix(
    header: String,
    mat: Vec<Vec<f32>>,
    ofile: &str,
    regions: Vec<Region>,
    scale_regions: &Scalingregions
) {
    // Write out the matrix to a compressed file.
    let omat = File::create(ofile).unwrap();
    let mut encoder = GzEncoder::new(omat, Compression::default());
    encoder.write_all(header.as_bytes()).unwrap();
    // Final check to make sure our regions and mat iter are of equal length.
    assert_eq!(regions.len(), mat.len());
    for (region, row) in regions.into_iter().zip(mat.into_iter()) {
        // Skipping rules.
        // skip_zeros
        if scale_regions.skipzero && row.iter().all(|&x| x == 0.0) {
            continue;
        }
        // min threshold
        if scale_regions.minthresh != 0.0 && row.iter().any(|&x| x < scale_regions.minthresh) {
            continue;
        }
        // max threshold
        if scale_regions.maxthresh != 0.0 && row.iter().any(|&x| x > scale_regions.maxthresh) {
            continue;
        }
        let mut writerow = format!(
            "{}\t{}\t{}\t{}:{}-{}\t{}\t{}\t",
            region.chrom,                // String field
            region.start.to_string(),    // u64 field 1 converted to string
            region.end.to_string(),    // u64 field 2 converted to string
            region.chrom,                // Chrom field of chr:st-end
            region.start.to_string(),    // st field of chr:st-end
            region.end.to_string(),    // end field of chr:st-end
            region.score,                // String field
            region.strand,                // String field
        );
        writerow.push_str(
            &row.iter().map(|x| (scale_regions.scale * x).to_string()).collect::<Vec<String>>().join("\t")
        );
        writerow.push_str("\n");
        encoder.write_all(writerow.as_bytes()).unwrap();
    }
}