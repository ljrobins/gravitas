#!/usr/bin/octave -qf
# EGM96 Gravity Model - Octave Implementation
# Computes Earth gravity acceleration up to degree and order 3
# Uses the same Pines normalization as the C++ implementation

# nm2i: maps (n,m) to 0-based linear index, same as C++
function idx = nm2i(n, m)
    idx = n * (n + 1) / 2 + m;
endfunction

function [cnm, snm] = load_egm96_coeffs(nmax)
    % Load EGM96 cosine and sine coefficients up to nmax
    % Returns arrays indexed by nm2i(n,m) + 1 (Octave 1-based)
    
    sz = nm2i(nmax + 2, nmax + 2) + 1;
    cnm = zeros(1, sz);
    snm = zeros(1, sz);

    # Boundary conditions
    cnm(nm2i(0, 0) + 1) = 1.0;
    snm(nm2i(0, 0) + 1) = 0.0;

    # n=1 coefficients (all zero in EGM96)
    cnm(nm2i(1, 0) + 1) = 0.0;
    snm(nm2i(1, 0) + 1) = 0.0;
    cnm(nm2i(1, 1) + 1) = 0.0;
    snm(nm2i(1, 1) + 1) = 0.0;

    # n=2 coefficients
    cnm(nm2i(2, 0) + 1) = -4.8416537174e-4;
    snm(nm2i(2, 0) + 1) = 0.0;
    cnm(nm2i(2, 1) + 1) = -1.8698763596e-10;
    snm(nm2i(2, 1) + 1) = 1.1952801203e-9;
    cnm(nm2i(2, 2) + 1) = 2.4391435240e-6;
    snm(nm2i(2, 2) + 1) = -1.4001668365e-6;

    # n=3 coefficients
    cnm(nm2i(3, 0) + 1) = 9.5725417379e-7;
    snm(nm2i(3, 0) + 1) = 0.0;
    cnm(nm2i(3, 1) + 1) = 2.0299888218e-6;
    snm(nm2i(3, 1) + 1) = 2.4851315872e-7;
    cnm(nm2i(3, 2) + 1) = 9.0462776860e-7;
    snm(nm2i(3, 2) + 1) = -6.1902594420e-7;
    cnm(nm2i(3, 3) + 1) = 7.2107265706e-7;
    snm(nm2i(3, 3) + 1) = 1.4143562696e-6;
endfunction

function g = pinesnorm(rf, cnm, snm, nmax, mu, req)
    % Implements the Pines normalization spherical harmonic gravity model
    % Same algorithm as in libgrav.hpp::pinesnorm()
    %
    % All internal arrays use 0-based indexing (stored in 1-based Octave arrays)
    
    rmag = norm(rf);
    stu = rf / rmag;
    # stu = [stu0, stu1, stu2] = [x/r, y/r, z/r]

    # Auxiliary coefficients array (Pines normalization)
    anm_sz = nm2i(nmax + 2, nmax + 2) + 1;
    anm = zeros(1, anm_sz);
    anm(1) = sqrt(2.0);  # anm[nm2i(0,0)] = sqrt(2)

    # Compute auxiliary coefficients (Pines normalization)
    for m = 0:nmax + 2
        if m != 0
            # Diagonal recursion: anm[m,m] = sqrt(1 + 1/(2m)) * anm[m-1,m-1]
            anm(nm2i(m, m) + 1) = sqrt(1.0 + 1.0 / (2.0 * m)) * anm(nm2i(m - 1, m - 1) + 1);
        end
        if m != nmax + 2
            # First off-diagonal recursion: anm[m+1,m] = sqrt(2m+3) * stu[2] * anm[m,m]
            anm(nm2i(m + 1, m) + 1) = sqrt(2 * m + 3) * stu(3) * anm(nm2i(m, m) + 1);
        end
        if m < nmax + 1
            for n = m + 2:nmax + 2
                alpha_num = (2 * n + 1) * (2 * n - 1);
                alpha_den = (n - m) * (n + m);
                alpha = sqrt(alpha_num / alpha_den);
                beta_num = (2 * n + 1) * (n - m - 1) * (n + m - 1);
                beta_den = (2 * n - 3) * (n + m) * (n - m);
                beta = sqrt(beta_num / beta_den);
                anm(nm2i(n, m) + 1) = alpha * stu(3) * anm(nm2i(n - 1, m) + 1) - beta * anm(nm2i(n - 2, m) + 1);
            end
        end
    end

    # Normalize first column (m=0) by sqrt(0.5)
    for n = 0:nmax + 2
        anm(nm2i(n, 0) + 1) *= sqrt(0.50);
    end

    # Compute trigonometric terms rm[m] and im[m] for m = 0..nmax+1
    # rm[m] and im[m] are 0-indexed in C++, stored at Octave index m+1
    max_m = nmax + 1;
    rm = zeros(1, max_m + 1);  # indices 0..max_m, stored at 1..max_m+1
    im = zeros(1, max_m + 1);
    rm(1) = 0.0;   # rm[0]
    rm(2) = 1.0;   # rm[1]
    im(1) = 0.0;   # im[0]
    im(2) = 0.0;   # im[1]

    for m = 2:max_m
        # rm[m] = stu[0]*rm[m-1] - stu[1]*im[m-1]
        rm(m + 1) = stu(1) * rm(m) - stu(2) * im(m);
        # im[m] = stu[0]*im[m-1] + stu[1]*rm[m-1]
        im(m + 1) = stu(1) * im(m) + stu(2) * rm(m);
    end

    # Scale factors
    rho = mu / (req * rmag);
    rhop = req / rmag;

    # Initialize gravity components
    g1 = 0.0;
    g2 = 0.0;
    g3 = 0.0;
    g4 = 0.0;

    # Sum over degree n
    for n = 0:nmax
        g1t = 0.0;
        g2t = 0.0;
        g3t = 0.0;
        g4t = 0.0;
        sm = 0.5;

        for m = 0:n
            # Get anm[n, m+1] for off-diagonal terms
            if n == m
                anmp1 = 0.0;
            else
                anmp1 = anm(nm2i(n, m + 1) + 1);
            end

            # Compute dnm, enm, fnm
            # rm[m] stored at Octave index m+1, rm[m+1] at m+2
            dnm = cnm(nm2i(n, m) + 1) * rm(m + 2) + snm(nm2i(n, m) + 1) * im(m + 2);
            enm = cnm(nm2i(n, m) + 1) * rm(m + 1) + snm(nm2i(n, m) + 1) * im(m + 1);
            fnm = snm(nm2i(n, m) + 1) * rm(m + 1) - cnm(nm2i(n, m) + 1) * im(m + 1);

            # Compute alpha
            alpha = sqrt(sm * (n - m) * (n + m + 1));

            # Accumulate gravity components
            g1t += anm(nm2i(n, m) + 1) * m * enm;
            g2t += anm(nm2i(n, m) + 1) * m * fnm;
            g3t += alpha * anmp1 * dnm;
            g4t += ((n + m + 1) * anm(nm2i(n, m) + 1) + alpha * stu(3) * anmp1) * dnm;

            if m == 0
                sm = 1.0;
            end
        end

        rho *= rhop;
        g1 += rho * g1t;
        g2 += rho * g2t;
        g3 += rho * g3t;
        g4 += rho * g4t;
    end

    # Compute final gravity vector
    g = zeros(1, 3);
    g(1) = g1 - g4 * stu(1);
    g(2) = g2 - g4 * stu(2);
    g(3) = g3 - g4 * stu(3);
endfunction

function g = earth_accel_octave(r_ecef, nmax)
    % Earth gravity acceleration using EGM96 model
    
    # Earth parameters (same as C++ implementation)
    mu = 398600.44;       # km^3/s^2
    req = 6378.137;       # km

    # Load coefficients
    [cnm, snm] = load_egm96_coeffs(nmax);

    # Compute gravity
    g = pinesnorm(r_ecef, cnm, snm, nmax, mu, req);
endfunction

function print_coefficients(nmax)
    % Print all coefficients used
    [cnm, snm] = load_egm96_coeffs(nmax);
    
    printf("\nCoefficients used (degree/order <= %d):\n", nmax);
    printf("%-4s %-4s %-20s %-20s\n", "n", "m", "cnm", "snm");
    printf("%s\n", repmat("-", 60));

    for n = 0:nmax
        for m = 0:n
            idx = nm2i(n, m);
            printf("%-4d %-4d %20.15e %20.15e\n", n, m, cnm(idx + 1), snm(idx + 1));
        endfor
    endfor
endfunction

# ===== Main =====
printf("EGM96 Gravity Model - Octave Implementation (degree/order <= 3)\n");
printf("==============================================================\n\n");

# Test points (ECEF coordinates in km)
test_points = [
    7000.0,   0.0,   0.0;
    7000.0, 7000.0,   0.0;
    6378.137,  0.0,   0.0;
      0.0,    0.0, 6378.137;
   26578.137,  0.0,   0.0
];

printf("Test Points:\n");
printf("%-10s %-10s %-10s | %-18s %-18s %-18s | %s\n", ...
    "X (km)", "Y (km)", "Z (km)", "gx (km/s^2)", "gy (km/s^2)", "gz (km/s^2)", "|g| (km/s^2)");
printf("%s\n", repmat("-", 95));

for i = 1:size(test_points, 1)
    r = test_points(i, :);
    g = earth_accel_octave(r, 3);
    gmag = norm(g);
    printf("%10.3f %10.3f %10.3f | %18.15e %18.15e %18.15e | %15.10f\n", ...
        r(1), r(2), r(3), g(1), g(2), g(3), gmag);
endfor

# Print coefficients
print_coefficients(3);

# ===== Validation against C++ reference values =====
printf("\nValidation against C++ implementation (gravitas.earth_acceleration with nmax=3):\n");
printf("============================================================================\n\n");

# Reference values from C++ implementation
ref_values = [
    -8.145722114672491e-03,  1.551220034299247e-08,  5.138614519466628e-09;
    -2.877958971145420e-03, -2.877992964222179e-03, -7.577989268589646e-09;
    -9.814266313695816e-03,  2.988041014616246e-08,  8.182834033573595e-09;
     1.288974930657953e-07,  1.582598136083982e-08, -9.766561041909736e-03;
    -5.643254384419217e-04, -1.101493608479816e-10,  6.494652192280322e-12
];

printf("%-10s %-10s %-10s | %-18s %-18s %-18s | %s\n", ...
    "X (km)", "Y (km)", "Z (km)", "gx (km/s^2)", "gy (km/s^2)", "gz (km/s^2)", "MATCH");
printf("%s\n", repmat("-", 1, 90));

all_pass = true;
for i = 1:size(test_points, 1)
    r = test_points(i, :);
    g = earth_accel_octave(r, 3);
    ref = ref_values(i, :);
    diff = abs(g - ref);
    max_diff = max(diff);
    if max_diff < 1e-15
        match = "YES";
    else
        match = "NO";
        all_pass = false;
    end
    printf("%10.3f %10.3f %10.3f | %18.15e %18.15e %18.15e | %s (diff: %e)\n", ...
        r(1), r(2), r(3), g(1), g(2), g(3), match, max_diff);
endfor

printf("\n");
if all_pass
    printf("ALL TESTS PASSED: Octave results match C++ implementation to machine precision!\n");
else
    printf("SOME TESTS FAILED: Differences detected.\n");
endif
