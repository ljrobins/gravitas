//! EGM96 gravity model coefficients
//!
//! Loads and stores the EGM96 spherical harmonic coefficients.
//! Pre-computes and caches coefficient arrays for all nmax values.

/// Maximum degree of EGM96
pub const EGM96_MAX_DEGREE: usize = 360;

/// Number of coefficients in EGM96
pub const EGM96_NCOEFF: usize = 65340;

/// EGM96 gravity model
#[derive(Clone)]
pub struct EGM96 {
    /// n values for each coefficient
    pub n: Vec<i32>,
    /// m values for each coefficient
    pub m: Vec<i32>,
    /// Cosine coefficients
    pub c: Vec<f64>,
    /// Sine coefficients
    pub s: Vec<f64>,
}

impl EGM96 {
    /// Load EGM96 coefficients from embedded data
    pub fn new() -> Self {
        Self {
            n: include_str!("n_EGM96.txt")
                .lines()
                .map(|line| line.trim().parse().expect("Failed to parse n_EGM96"))
                .collect(),
            m: include_str!("m_EGM96.txt")
                .lines()
                .map(|line| line.trim().parse().expect("Failed to parse m_EGM96"))
                .collect(),
            c: include_str!("c_EGM96.txt")
                .lines()
                .map(|line| line.trim().parse().expect("Failed to parse c_EGM96"))
                .collect(),
            s: include_str!("s_EGM96.txt")
                .lines()
                .map(|line| line.trim().parse().expect("Failed to parse s_EGM96"))
                .collect(),
        }
    }

    /// Get the maximum degree present in the model
    pub fn max_degree(&self) -> usize {
        self.n.iter().map(|&n| n as usize).max().unwrap_or(0)
    }

    /// Pre-compute and cache coefficient arrays for all nmax values
    ///
    /// Returns a vector where index i contains the (cnm, snm) arrays for nmax=i
    pub fn precompute_all_nmax(&self) -> Vec<(Vec<f64>, Vec<f64>)> {
        let mut cache = Vec::with_capacity(EGM96_MAX_DEGREE + 1);
        
        // Pre-allocate the full coefficient arrays
        let full_sz = nm2i(EGM96_MAX_DEGREE + 2, EGM96_MAX_DEGREE + 2);
        let mut full_cnm = vec![0.0f64; full_sz];
        let mut full_snm = vec![0.0f64; full_sz];
        
        // Load all coefficients into the full arrays
        let num = self.n.len();
        for i in 0..num {
            let n = self.n[i] as usize;
            let m = self.m[i] as usize;
            let ind = nm2i(n, m);
            full_cnm[ind] = self.c[i];
            full_snm[ind] = self.s[i];
        }
        
        // Set the monopole term
        full_snm[0] = 0.0;
        full_cnm[0] = 1.0;
        
        // For each nmax, create a slice of the full arrays
        for nmax in 0..=EGM96_MAX_DEGREE {
            let sz = nm2i(nmax + 2, nmax + 2);
            let cnm = full_cnm[..sz].to_vec();
            let snm = full_snm[..sz].to_vec();
            cache.push((cnm, snm));
        }
        
        cache
    }

    /// Read coefficients up to a given nmax, storing in flat arrays
    ///
    /// The coefficients are stored in a flat array indexed by nm2i(n, m) = n*(n+1)/2 + m
    ///
    /// # Arguments
    /// * `nmax` - Maximum degree to include
    ///
    /// # Returns
    /// (cnm, snm) - Flat arrays of cosine and sine coefficients
    ///
    /// Note: Only coefficients up to the given nmax are loaded. The loop breaks
    /// as soon as m == nmax is encountered, ensuring minimal coefficient loading.
    pub fn read_cnm_snm(&self, nmax: usize) -> (Vec<f64>, Vec<f64>) {
        let sz = nm2i(nmax + 2, nmax + 2);
        let mut cnm = vec![0.0f64; sz];
        let mut snm = vec![0.0f64; sz];

        // Only iterate through coefficients up to the point where m == nmax
        // This ensures we only load the coefficients we need for the given nmax
        let num = self.n.len();
        for i in 0..num {
            let n = self.n[i] as usize;
            let m = self.m[i] as usize;
            let ind = nm2i(n, m);
            cnm[ind] = self.c[i];
            snm[ind] = self.s[i];
            // Stop when we've read all coefficients for this nmax
            if m == nmax {
                break;
            }
        }

        // Set the monopole term
        snm[0] = 0.0;
        cnm[0] = 1.0;

        (cnm, snm)
    }
}

/// Convert (n, m) to flat array index: n*(n+1)/2 + m
#[inline]
pub const fn nm2i(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}
