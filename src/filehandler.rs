use rust_htslib::bam::{Read, Reader};
use itertools::Itertools;
use std::io::{BufReader, BufWriter, Write};
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use bigtools::{BigWigRead, BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use std::collections::HashMap;
use crate::computematrix::Scalingregions;
use flate2::write::GzEncoder;
use flate2::Compression;

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

pub fn read_bedfiles(bed_files: &Vec<String>) -> (Vec<(String, u64, u64, String, String, String)>, HashMap<String, u64>) {
    // read all bedfiles in a Vec of strings (filepaths)
    // returns a vec of tuples, with last entry being the filepath (without extension and dirs ~= 'smartLabels')
    let mut regionsizes: HashMap<String, u64> = HashMap::new();
    let mut regions: Vec<(String, u64, u64, String, String, String)> = Vec::new();
    for bed in bed_files {
        let entryname = Path::new(bed)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();

        let bedfile = BufReader::new(File::open(bed).unwrap());
        let mut entries: u64 = 0;
        for line in bedfile.lines() {
            let line = line.unwrap();
            let fields: Vec<&str> = line.split('\t').collect();
            // If score and strand are not provided, we set them to .
            let field_3 = fields.get(3).unwrap_or(&".").to_string();
            let field_4 = fields.get(4).unwrap_or(&".").to_string();
            regions.push(
                (
                    fields[0].to_string(), //chrom
                    fields[1].parse().unwrap(), //start
                    fields[2].parse().unwrap(), //end
                    field_3, //score
                    field_4, //strand
                    entryname.to_string() //bedfile_name
                )
            );
            entries += 1;
        }
        regionsizes.insert(entryname, entries);
    }
    return (regions, regionsizes);
}

pub fn chrombounds_from_bw(bwfile: &str) -> HashMap<String, u64> {
    // define chromsizes hashmap
    let mut chromsizes: HashMap<String, u64> = HashMap::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let mut reader = BigWigRead::open(bwf).unwrap();
    for chrom in reader.chroms() {
        chromsizes.insert(chrom.name.clone(), chrom.length as u64);
    }
    chromsizes
}

pub fn bwintervals(
    bwfile: &str,
    regions: &Vec<(String, u64, u64, String, String, String)>,
    slopregions: &Vec<Vec<(u64, u64)>>,
    scale_regions: &Scalingregions
) -> Vec<Vec<f64>> {
    // For a given bw file, a vector of slop regions (where every vec entry is a vec of (start, end) tuples)
    // return a vector with for every region a vector of f64.

    // Make sure regions and slopregions are of equal length
    assert_eq!(regions.len(), slopregions.len());
    
    // Define return vector, set up bw reader.
    let mut bwvals: Vec<Vec<f64>> = Vec::new();
    let bwf = File::open(bwfile).expect("Failed to open bw file.");
    let mut reader = BigWigRead::open(bwf).unwrap();

    // Iterate over regions and slopregions synchronously.
    for (sls, (chrom, _,  _, _, _, _)) in slopregions.iter().zip(regions.iter()) {
        let mut bwval: Vec<f64> = Vec::new();
        // get 'min and max' to query
        let (min, max) = sls
            .iter()
            .flat_map(|(x, y)| vec![*x, *y])
            .filter(|&value| value != 0)
            .minmax()
            .into_option()
            .expect("Vec is empty.");
        let binvals = reader.get_interval(chrom, min as u32, max as u32).unwrap();
        // since binvals (can) be over binsizes, we expand them to bp and push them to a hashmap
        let mut bwhash: HashMap<u64, f64> = HashMap::new();
        for interval in binvals {
            let interval = interval.unwrap();
            let start = interval.start as u64;
            let end = interval.end as u64;
            let val = interval.value as f64;
            bwhash.extend((start..end).map(|bp| (bp, val)));
        }
        for (start, end) in sls {
            if start == end && *end == 0 {
                if scale_regions.missingdata_as_zero {
                    bwval.push(0.0);
                } else {
                    bwval.push(std::f64::NAN);
                }
            } else {
                // Get values from the hashmap
                let mean: f64 = (*start..*end)
                    .filter_map(|bp| bwhash.get(&bp))
                    .copied()
                    .sum::<f64>()
                    / (end - start) as f64;
                bwval.push(mean);
            }
        }
        // Make sure bwval is of expected length.
        assert_eq!(bwval.len(), scale_regions.cols_expected / scale_regions.bwfiles);
        bwvals.push(bwval);
    }
    bwvals
}

pub fn header_matrix(scale_regions: &Scalingregions, regionsizes: HashMap<String, u64>) -> String {
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
    let mut groupbounds: Vec<u64> = Vec::new();
    groupbounds.push(0);
    let mut cumsum: u64 = 0;
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

pub fn write_matrix(header: String, mat: Vec<Vec<f64>>, ofile: &str, regions: Vec<(String, u64, u64, String, String, String)>) {
    println!("Writing out matrix to file.");
    // Write out the matrix to a compressed file.
    let omat = File::create(ofile).unwrap();
    let mut encoder = GzEncoder::new(omat, Compression::default());
    encoder.write_all(header.as_bytes()).unwrap();
    // Final check to make sure our regions and mat iter are of equal length.
    assert_eq!(regions.len(), mat.len());
    for (region, row) in regions.into_iter().zip(mat.into_iter()) {
        let mut writerow = format!(
            "{}\t{}\t{}\t{}\t{}\t",
            region.0,                // String field
            region.1.to_string(),    // u64 field 1 converted to string
            region.2.to_string(),    // u64 field 2 converted to string
            region.3,                // String field
            region.4,                // String field
        );
        writerow.push_str(
            &row.iter().map(|x| x.to_string()).collect::<Vec<String>>().join("\t")
        );
        writerow.push_str("\n");
        encoder.write_all(writerow.as_bytes()).unwrap();
    }
}