use pyo3::prelude::*;
mod bamcoverage;
mod bamcompare;
// mod computematrix;
mod covcalc;
mod alignmentsieve;
mod filehandler;
mod normalization;
mod calc;

#[pymodule]
fn hp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bamcoverage::r_bamcoverage, m)?)?;
    m.add_function(wrap_pyfunction!(bamcompare::r_bamcompare, m)?)?;
    //m.add_function(wrap_pyfunction!(computematrix::r_computematrix, m)?)?;
    Ok(())
}
