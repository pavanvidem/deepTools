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