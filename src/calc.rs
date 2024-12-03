pub fn median(mut nvec: Vec<u32>) -> f32 {
    if nvec.is_empty() {
        return 0.0;
    } else if nvec.len() == 1 {
        return nvec[0] as f32;
    } else {
        nvec.sort_unstable();
        let len = nvec.len();
        if len % 2 == 1 {
            return nvec[len / 2] as f32;
        } else {
            return (nvec[len / 2] + nvec[len / 2 - 1]) as f32 / 2.0;
        }
    }
}


pub fn calc_ratio(
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