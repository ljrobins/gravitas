//! Gravity field computation using the Pines norm recurrence algorithm
//!
//! Based on the implementation from:
//! https://core.ac.uk/download/pdf/76424485.pdf

use crate::egm96::{EGM96, EGM96_MAX_DEGREE};
use rayon::prelude::*;
use std::f64::consts::SQRT_2;
use std::sync::LazyLock;

/// Earth parameters for EGM96
const EARTH_MU: f64 = 398600.44;   // km^3/s^2
const EARTH_REQ: f64 = 6378.137;   // km

/// EGM96 with pre-computed coefficient cache
struct EGM96WithCache {
    egm: EGM96,
    cache: Vec<(Vec<f64>, Vec<f64>)>,
}

impl EGM96WithCache {
    /// Get pre-computed coefficient arrays for a given nmax
    fn get_cnm_snm(&self, nmax: usize) -> (&[f64], &[f64]) {
        let (cnm, snm) = &self.cache[nmax];
        (cnm.as_slice(), snm.as_slice())
    }
}

/// Thread-local cache for work arrays to avoid repeated allocation
struct WorkArrays {
    anm: Vec<f64>,
    rm: Vec<f64>,
    im: Vec<f64>,
}

impl WorkArrays {
    fn new() -> Self {
        Self {
            anm: Vec::new(),
            rm: Vec::new(),
            im: Vec::new(),
        }
    }

    /// Get work arrays sized for the given nmax, reusing if possible
    fn get(&mut self, nmax: usize) -> (&mut [f64], &mut [f64], &mut [f64]) {
        let nmax_plus_2 = nmax + 2;
        let anm_sz = nm2i(nmax_plus_2 + 1, nmax_plus_2 + 1);
        let rm_sz = nmax + 2;

        // Resize vectors if needed
        if self.anm.len() < anm_sz {
            self.anm.resize(anm_sz, 0.0);
        }
        if self.rm.len() < rm_sz {
            self.rm.resize(rm_sz, 0.0);
        }
        if self.im.len() < rm_sz {
            self.im.resize(rm_sz, 0.0);
        }

        (&mut self.anm[..anm_sz], &mut self.rm[..rm_sz], &mut self.im[..rm_sz])
    }
}

thread_local! {
    static WORK_ARRAYS: std::cell::RefCell<WorkArrays> = std::cell::RefCell::new(WorkArrays::new());
}

/// Cached EGM96 coefficients with pre-computed arrays for all nmax values
static EGM: LazyLock<EGM96WithCache> = LazyLock::new(|| {
    let egm = EGM96::new();
    let cache = egm.precompute_all_nmax();
    EGM96WithCache { egm, cache }
});

/// Compute gravity acceleration for a single position vector
///
/// # Arguments
/// * `r` - Position vector in ECEF coordinates (km)
/// * `cnm` - Cosine coefficients (flat array)
/// * `snm` - Sine coefficients (flat array)
/// * `nmax` - Maximum degree of spherical harmonic expansion
///
/// # Returns
/// Acceleration vector (km/s^2)
pub fn compute_acceleration_single(
    r: &[f64; 3],
    cnm: &[f64],
    snm: &[f64],
    nmax: usize,
) -> [f64; 3] {
    let rmag = r[0].hypot(r[1].hypot(r[2]));

    if rmag == 0.0 {
        return [0.0, 0.0, 0.0];
    }

    let inv_rmag = 1.0 / rmag;
    let stu = [r[0] * inv_rmag, r[1] * inv_rmag, r[2] * inv_rmag];
    let nmax_plus_2 = nmax + 2;

    // Use thread-local work arrays to avoid allocation
    // We use a closure to avoid holding the borrow across the entire function
    let result = WORK_ARRAYS.with(|wa| {
        let mut wa = wa.borrow_mut();
        let (anm, rm, im) = wa.get(nmax);
        
        // Only zero the parts we actually use (not the entire array)
        let anm_sz = nm2i(nmax_plus_2 + 1, nmax_plus_2 + 1);
        let rm_sz = nmax + 2;
        
        // Zero anm array (only the part we use)
        for i in 0..anm_sz {
            anm[i] = 0.0;
        }
        // Zero rm and im arrays (only the part we use)
        for i in 0..rm_sz {
            rm[i] = 0.0;
            im[i] = 0.0;
        }
        
        anm[0] = SQRT_2;

        // Diagonal and off-diagonal recursions for normalization
        for m in 0..=nmax_plus_2 {
            if m != 0 {
                anm[nm2i(m, m)] =
                    (1.0 + 1.0 / (2.0 * m as f64)).sqrt() * anm[nm2i(m - 1, m - 1)];
            }
            if m != nmax_plus_2 {
                anm[nm2i(m + 1, m)] = (2.0 * m as f64 + 3.0).sqrt() * stu[2] * anm[nm2i(m, m)];
            }
            if m < nmax + 1 {
                for n in (m + 2)..=nmax_plus_2 {
                    let alpha = ((2.0 * n as f64 + 1.0) * (2.0 * n as f64 - 1.0)
                        / ((n as f64 - m as f64) * (n as f64 + m as f64))).sqrt();
                    let beta = ((2.0 * n as f64 + 1.0)
                        * (n as f64 - m as f64 - 1.0)
                        * (n as f64 + m as f64 - 1.0)
                        / ((2.0 * n as f64 - 3.0)
                            * (n as f64 + m as f64)
                            * (n as f64 - m as f64))).sqrt();
                    anm[nm2i(n, m)] = alpha * stu[2] * anm[nm2i(n - 1, m)]
                        - beta * anm[nm2i(n - 2, m)];
                }
            }
        }

        for n in 0..=nmax_plus_2 {
            anm[nm2i(n, 0)] *= 0.5_f64.sqrt();
        }

        rm[0] = 0.0;
        rm[1] = 1.0;
        im[0] = 0.0;
        im[1] = 0.0;

        for m in 2..nmax + 1 {
            rm[m] = stu[0] * rm[m - 1] - stu[1] * im[m - 1];
            im[m] = stu[0] * im[m - 1] + stu[1] * rm[m - 1];
        }

        let rho = EARTH_MU / (EARTH_REQ * rmag);
        let rhop = EARTH_REQ / rmag;

        let mut g1: f64 = 0.0;
        let mut g2: f64 = 0.0;
        let mut g3: f64 = 0.0;
        let mut g4: f64 = 0.0;

        for n in 0..=nmax {
            let mut g1t: f64 = 0.0;
            let mut g2t: f64 = 0.0;
            let mut g3t: f64 = 0.0;
            let mut g4t: f64 = 0.0;
            let mut sm = 0.5;

            for m in 0..=n {
                let anmp1 = if n == m { 0.0 } else { anm[nm2i(n, m + 1)] };
                let dnm = cnm[nm2i(n, m)] * rm[m + 1] + snm[nm2i(n, m)] * im[m + 1];
                let enm = cnm[nm2i(n, m)] * rm[m] + snm[nm2i(n, m)] * im[m];
                let fnm = snm[nm2i(n, m)] * rm[m] - cnm[nm2i(n, m)] * im[m];
                let alpha = (sm * (n as f64 - m as f64) * (n as f64 + m as f64 + 1.0)).sqrt();

                g1t += anm[nm2i(n, m)] * m as f64 * enm;
                g2t += anm[nm2i(n, m)] * m as f64 * fnm;
                g3t += alpha * anmp1 * dnm;
                g4t += ((n + m + 1) as f64 * anm[nm2i(n, m)]
                    + alpha * stu[2] * anmp1) * dnm;

                if m == 0 { sm = 1.0; }
            }

            let rho_scaled = rho * rhop;
            g1 += rho_scaled * g1t;
            g2 += rho_scaled * g2t;
            g3 += rho_scaled * g3t;
            g4 += rho_scaled * g4t;
        }

        [g1 - g4 * stu[0], g2 - g4 * stu[1], g3 - g4 * stu[2]]
    });
    
    result
}

/// Threshold below which threading overhead outweighs benefits
/// Below this nmax, serial computation is faster than parallel
const THREADING_THRESHOLD: usize = 16;

/// Compute gravity acceleration for multiple position vectors
///
/// # Arguments
/// * `positions` - Array of position vectors, each of shape [x, y, z] in km
/// * `nmax` - Maximum degree of spherical harmonic expansion
///
/// # Returns
/// Array of acceleration vectors (km/s^2), same length as positions
///
/// # Panics
/// Panics if nmax > EGM96_MAX_DEGREE
pub fn compute_acceleration_batch(positions: &[[f64; 3]], nmax: usize) -> Vec<[f64; 3]> {
    if nmax > EGM96_MAX_DEGREE {
        panic!(
            "nmax ({}) exceeds maximum degree of EGM96 ({})",
            nmax, EGM96_MAX_DEGREE
        );
    }

    // Special case: nmax=0 is just the monopole term (mu/r^2)
    // Skip coefficient loading and parallel iteration for maximum performance
    if nmax == 0 {
        return positions
            .iter()
            .map(|r| {
                let rmag = r[0].hypot(r[1].hypot(r[2]));
                if rmag == 0.0 {
                    return [0.0, 0.0, 0.0];
                }
                let inv_r3 = 1.0 / (rmag * rmag * rmag);
                [
                    -EARTH_MU * r[0] * inv_r3,
                    -EARTH_MU * r[1] * inv_r3,
                    -EARTH_MU * r[2] * inv_r3,
                ]
            })
            .collect();
    }

    // Use pre-computed cached EGM coefficients (zero allocation)
    let (cnm, snm) = EGM.get_cnm_snm(nmax);

    // Use threading only when nmax is high enough to justify the overhead
    // Below THREADING_THRESHOLD, serial computation is faster
    if nmax < THREADING_THRESHOLD {
        positions
            .iter()
            .map(|r| compute_acceleration_single(r, &cnm, &snm, nmax))
            .collect()
    } else {
        positions
            .par_iter()
            .map(|r| compute_acceleration_single(r, &cnm, &snm, nmax))
            .collect()
    }
}

/// Compute Earth's gravity acceleration using EGM96
///
/// # Arguments
/// * `positions` - Array of position vectors in ECEF coordinates (km), shape [N, 3]
/// * `nmax` - Maximum degree of spherical harmonic expansion (0 <= nmax <= 360)
///
/// # Returns
/// Acceleration vectors in ECEF coordinates (km/s^2), shape [N, 3]
///
/// # Panics
/// Panics if nmax < 0 or nmax > 360
pub fn earth_acceleration(positions: &[[f64; 3]], nmax: usize) -> Vec<[f64; 3]> {
    compute_acceleration_batch(positions, nmax)
}

/// Convert (n, m) to flat array index: n*(n+1)/2 + m
#[inline]
pub const fn nm2i(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nm2i() {
        assert_eq!(nm2i(0, 0), 0);
        assert_eq!(nm2i(1, 0), 1);
        assert_eq!(nm2i(1, 1), 2);
        assert_eq!(nm2i(2, 0), 3);
        assert_eq!(nm2i(2, 1), 4);
        assert_eq!(nm2i(2, 2), 5);
    }

    #[test]
    fn test_single_point_nmax0() {
        let positions = [[6778.0, 0.0, 0.0]];
        let result = earth_acceleration(&positions, 0);
        // Compare with known C++ output: [-0.0086763, 0., 0.]
        assert!((result[0][0] - (-0.0086763)).abs() < 1e-8);
        assert!((result[0][1] - 0.0).abs() < 1e-15);
        assert!((result[0][2] - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_single_point_nmax10() {
        let positions = [[6778.0, 0.0, 0.0]];
        let result = earth_acceleration(&positions, 10);
        // Verify result is physically reasonable (monopole + small perturbations)
        let expected_monopole = -EARTH_MU / (6778.0 * 6778.0);
        assert!((result[0][0] - expected_monopole).abs() < 1e-4);
    }

    #[test]
    fn test_single_point_nmax16() {
        let positions = [[7000.0, 0.0, 0.0]];
        let result = earth_acceleration(&positions, 16);
        // Verify result is physically reasonable
        let expected_monopole = -EARTH_MU / (7000.0 * 7000.0);
        assert!((result[0][0] - expected_monopole).abs() < 1e-4);
    }

    #[test]
    fn test_surface_north_pole() {
        let positions = [[0.0, 0.0, 6378.0]];
        let result = earth_acceleration(&positions, 0);
        // nmax=0 should give pure gravitational acceleration at pole
        assert!((result[0][0] - 0.0).abs() < 1e-15);
        assert!((result[0][1] - 0.0).abs() < 1e-15);
        assert!((result[0][2] - (-0.00979871)).abs() < 1e-8);
    }

    #[test]
    fn test_geostationary() {
        let positions = [[42164.0, 0.0, 0.0]];
        let result = earth_acceleration(&positions, 10);
        // Verify result is physically reasonable
        let expected_monopole = -EARTH_MU / (42164.0 * 42164.0);
        assert!((result[0][0] - expected_monopole).abs() < 1e-5);
    }

    #[test]
    fn test_multi_point() {
        let positions = [
            [6778.0, 0.0, 0.0],
            [0.0, 6778.0, 0.0],
            [0.0, 0.0, 6778.0],
        ];
        let result = earth_acceleration(&positions, 5);
        assert_eq!(result.len(), 3);
        // All should have similar magnitude (same radius)
        let mag0 = result[0][0].hypot(result[0][1].hypot(result[0][2]));
        let mag1 = result[1][0].hypot(result[1][1].hypot(result[1][2]));
        let mag2 = result[2][0].hypot(result[2][1].hypot(result[2][2]));
        assert!((mag0 - mag1).abs() < 1e-4);
        assert!((mag1 - mag2).abs() < 1e-4);
    }

    #[test]
    fn test_large_batch_consistency() {
        // Create 1000 identical points
        let positions: Vec<[f64; 3]> = vec![[7000.0, 0.0, 0.0]; 1000];
        let result = earth_acceleration(&positions, 16);
        assert_eq!(result.len(), 1000);

        // All results should be identical
        for i in 1..result.len() {
            assert_eq!(result[i], result[0]);
        }
    }

    #[test]
    fn test_nmax_boundaries() {
        let positions = [[6778.0, 0.0, 0.0]];

        // nmax=0 should work
        let _ = earth_acceleration(&positions, 0);

        // nmax=360 should work
        let _ = earth_acceleration(&positions, 360);

        // nmax=361 should panic
        let result = std::panic::catch_unwind(|| earth_acceleration(&positions, 361));
        assert!(result.is_err(), "nmax=361 should panic");
    }

    #[test]
    fn test_various_positions() {
        // Test various positions and compare with known C++ behavior
        let test_cases = [
            // (position, nmax, expected_x, tolerance)
            ([6378.0, 0.0, 0.0], 0, -0.00979871, 1e-8),
            ([4500.0, 4500.0, 4500.0], 0, -0.00378818, 1e-7),
            ([100000.0, 0.0, 0.0], 0, -3.9860044e-5, 1e-10),
        ];

        for (pos, nmax, expected_x, tol) in test_cases {
            let result = earth_acceleration(&[pos], nmax);
            assert!(
                (result[0][0] - expected_x).abs() < tol,
                "Position {:?}, nmax={}: got x={}, expected {} (tol={})",
                pos,
                nmax,
                result[0][0],
                expected_x,
                tol
            );
        }
    }

    #[test]
    fn test_symmetric_positions() {
        // Test that opposite positions give opposite accelerations
        let positions_pos = [[6778.0, 0.0, 0.0]];
        let positions_neg = [[-6778.0, 0.0, 0.0]];

        for nmax in [0, 10, 100, 360] {
            let result_pos = earth_acceleration(&positions_pos, nmax);
            let result_neg = earth_acceleration(&positions_neg, nmax);

            for i in 0..3 {
                assert!(
                    (result_pos[0][i] + result_neg[0][i]).abs() < 1e-6,
                    "Anti-symmetry failed at nmax={}, component={}: {} vs {}",
                    nmax,
                    i,
                    result_pos[0][i],
                    result_neg[0][i]
                );
            }
        }
    }

    #[test]
    fn test_very_far_position() {
        let positions = [[100000.0, 0.0, 0.0]];
        let result = earth_acceleration(&positions, 360);
        // At large distances, acceleration should approach mu/r^2
        let r = 100000.0;
        let expected_mag = EARTH_MU / (r * r);
        let actual_mag = result[0][0].abs();
        assert!((actual_mag - expected_mag).abs() < 1e-6);
    }
}
