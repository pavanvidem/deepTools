pub fn scale_factor(norm_method: &str, mapped: u64, binsize: u64, effective_genome_size) -> f64 {
    let mut scale_factor = 1.0;
    return match norm_method {
        "RPKM" => {
            // RPKM = # reads per tile / total reads (millions) * tile length (kb)
            let mmr = mapped as f64 / 1e6;
            let bs_kb = binsize as f64 / 1000.0;
            scale_factor *= 1.0 / (mmr * bs_kb);
            scale_factor
        }
        "CPM" => {
            // CPM = # reads per tile / total reads (millions)
            let mmr = mapped as f64 / 1e6;
            scale_factor *= 1.0 / mmr;
            scale_factor
        }
        "BPM" => {
            // BPM = bins per million mapped reads
            let bs_kb: f64 = binsize as f64 / 1000.0;
            let tmp_scalefactor = (mapped as f64 / bs_kb) / 1e6;
            scale_factor *= 1.0 / (tmp_scalefactor * bs_kb);
            scale_factor
        }
        "RPGC" => {
            // RPGC = mapped reads * fragment length / effective genome size
            scale_factor
        }
        _ => {
            scale_factor
        }
    }
}