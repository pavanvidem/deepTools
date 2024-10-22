use pyo3::prelude::*;
mod bamcoverage;
mod compute;
mod alignmentsieve;

#[pymodule]
fn hp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bamcoverage::r_bamcoverage, m)?)?;
    Ok(())
}
