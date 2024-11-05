pub fn scale_factor(norm_method: &str, mapped: u64, binsize: u64) -> f64 {
    let mut scale_factor = 1.0;
    return match norm_method {
        "RPKM" => {
            // RPKM = # reads per tile / total reads (millions) * tile length (kb)
            let mmr = mapped as f64 / 1e6;
            let bs_kb = binsize as f64 / 1000;
            scale_factor *= 1.0 / (mmr * bs_kb)
        }
        _ => {}
    }
    scale_factor
}