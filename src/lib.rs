use pyo3::prelude::*;
mod bamcoverage;
mod covcalc;
mod alignmentsieve;
mod bamhandler;
mod normalization;

#[pymodule]
fn hp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bamcoverage::r_bamcoverage, m)?)?;
    Ok(())
}
