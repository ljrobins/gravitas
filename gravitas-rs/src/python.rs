//! Python bindings for gravitas using pyo3

use pyo3::prelude::*;
use pyo3::types::PyModule;
use numpy::{PyArray2, PyArrayMethods, PyReadonlyArray2};

use crate::gravity::earth_acceleration as rust_earth_acceleration;

/// Acceleration due to the Earth's gravity field, as defined by the EGM96 model.
///
/// :param r_ecef_km: Position in ECEF coordinates, in km.
/// :type r_ecef_km: ndarray [nx3]
/// :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
/// :type nmax: int
/// :return: Acceleration in ECI coordinates, in km/s^2.
/// :rtype: ndarray [nx3]
#[pyfunction]
fn earth_acceleration(
    r_ecef_km: PyReadonlyArray2<'_, f64>,
    nmax: usize,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let view = r_ecef_km.as_array();
    let (n, ncols) = view.dim();

    if ncols != 3 {
        return Err(pyo3::PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input array must have shape (n, 3)"
        ));
    }

    // Convert input to Rust format
    let positions: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let row = view.row(i);
            [row[0], row[1], row[2]]
        })
        .collect();

    // Compute acceleration
    let result = rust_earth_acceleration(&positions, nmax);

    // Convert to numpy array - create 2D array row by row
    let gil = unsafe { pyo3::Python::assume_gil_acquired() };
    let mut output = PyArray2::<f64>::zeros_bound(gil, (n, 3), false);
    for i in 0..n {
        let r = &result[i];
        unsafe {
            *output.uget_mut([i, 0]) = r[0];
            *output.uget_mut([i, 1]) = r[1];
            *output.uget_mut([i, 2]) = r[2];
        }
    }
    Ok(output)
}

/// Moon acceleration (placeholder - not yet implemented)
#[pyfunction]
fn moon_acceleration(_r_mcmf_km: PyReadonlyArray2<'_, f64>, _nmax: usize) -> PyResult<Bound<'_, PyArray2<f64>>> {
    Err(pyo3::PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
        "Moon acceleration not yet implemented"
    ))
}

/// A Python module implemented in Rust
#[pymodule]
fn gravitas(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(earth_acceleration, m)?)?;
    m.add_function(wrap_pyfunction!(moon_acceleration, m)?)?;
    Ok(())
}
