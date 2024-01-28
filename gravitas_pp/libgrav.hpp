#include <GRGM360.hpp>
#include <EGM96.hpp>
#include <MRO120F.hpp>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Python.h>

using namespace std;

void print_test()
{
    cout << "20" << endl;
}


double req;
double mu;
int model_index;
int body_index;

enum BODY {
    EARTH,
    MOON,
    MARS,
};

enum MODEL {
    EGM96,
    GRGM360, // https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shbdr/gggrx_1200a_shb_l180.lbl
    MRO120F, // https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_120f_sha.lbl
};

void set_indices(char* model_name, int *model_index, int *body_index) {
    
    if (strcmp(model_name, "EGM96") == 0) {
        *model_index = EGM96;
        *body_index = EARTH;
        return;
    }
    if (strcmp(model_name, "GRGM360") == 0) {
        *model_index = GRGM360;
        *body_index = MOON;
        return;
    }
    if (strcmp(model_name, "MRO120F") == 0) {
        *model_index = MRO120F;
        *body_index = MARS;
        return;
    }
    // raise a python error
    PyErr_SetString(PyExc_ValueError, "Invalid model name");
    return;
}


void set_body_params(int body_index, double *mu, double *req) {
    if(body_index == EARTH) {
        *mu = 398600.44;
        *req = 6378.137;
    }
    if(body_index == MOON) {
        *mu = 4902.8001224453001;
        *req = 1738.0;
    }
    if(body_index == MARS) {
        *mu = 42828.3748574;
        *req = 3396.0;
    }
}

int nm2i(int n, int m) {
    return n * (n+1) / 2 + m;
}

void read_cnm_snm(int nmax, int model_index, Eigen::VectorXd &cnm, Eigen::VectorXd &snm) {
    // printf("Starting coefficients read!\n");

    // set pointers
    const Eigen::Map<const Eigen::VectorXi> *n = nullptr;
    const Eigen::Map<const Eigen::VectorXi> *m = nullptr;
    const Eigen::Map<const Eigen::VectorXd> *c = nullptr;
    const Eigen::Map<const Eigen::VectorXd> *s = nullptr;
    int num = 0;

    switch (model_index) {
        case EGM96:
            n = &n_EGM96_eigen; m = &m_EGM96_eigen; c = &c_EGM96_eigen; s = &s_EGM96_eigen;
            num = n_EGM96_eigen.size();
            break;
        case GRGM360:
            break;
    }

    // if nmax is greater than the maximum degree of the model, raise an error
    if(nmax > (*n)[num-1]) {
        PyErr_SetString(PyExc_ValueError, "nmax is greater than the maximum degree of the model");
        return;
    }

    for(int i = 0; i < num; i++) {
        int ind = nm2i((*n)[i], (*m)[i]);
        cnm[ind] = (*c)[i];
        snm[ind] = (*s)[i];
        // printf("n=%d, m=%d, c=%.2e, s=%.2e, cnm=%.2e, snm=%.2e\n", (*n)[i], (*m)[i], (*c)[i], (*s)[i], cnm[ind], snm[ind]);
        if((*m)[i] == nmax) {
            break;
        }
    }

    snm[0] = 0.0;
    cnm[0] = 1.0;
    // printf("Finished coefficients read!\n");
    return;
}

Eigen::VectorXd pinesnorm(Eigen::VectorXd rf, Eigen::VectorXd cnm, Eigen::VectorXd snm, int nmax, double mu, double req) {
    // printf("Starting pinesnorm!\n");   
    // Based on pinesnorm() from: https://core.ac.uk/download/pdf/76424485.pdf

    double rmag = rf.norm();
    Eigen::VectorXd stu = rf.normalized();
    int anm_sz = nm2i(nmax+3, nmax+3);
    Eigen::VectorXd anm(anm_sz);
    anm[0] = sqrt(2.0);

    for(int m = 0; m <= nmax+2; m++) {
        if(m != 0) { // DIAGONAL RECURSION
            anm[nm2i(m,m)] = sqrt(1.0+1.0/(2.0*m))*anm[nm2i(m-1,m-1)];
            // printf("ANM: %d %d %.2e\n", m, m, anm[nm2i(m,m)]);
        }
        if(m != nmax+2) { // FIRST OFF-DIAGONAL RECURSION 
            anm[nm2i(m+1,m)] = sqrt(2*m+3)*stu[2]*anm[nm2i(m,m)];
        }
        if(m < nmax+1) {
            for(int n = m+2; n <= nmax+2; n++) {
                double alpha_num = (2*n+1)*(2*n-1);
                double alpha_den = (n-m)*(n+m);
                double alpha = sqrt(alpha_num/alpha_den);
                double beta_num = (2*n+1)*(n-m-1)*(n+m-1);
                double beta_den = (2*n-3)*(n+m)*(n-m);
                double beta = sqrt(beta_num/beta_den);
                anm[nm2i(n,m)] = alpha*stu[2]*anm[nm2i(n-1,m)] - beta*anm[nm2i(n-2,m)];
            }
        }
    }
    
    for(int n = 0; n <= nmax+2; n++) {
        anm[nm2i(n,0)] *= sqrt(0.50);
    }
     
    Eigen::VectorXd rm(nmax+2);
    Eigen::VectorXd im(nmax+2);
    rm[0] = 0.00; rm[1] = 1.00; 
    im[0] = 0.00; im[1] = 0.00; 
    for(int m = 2; m < nmax+2; m++) {
        rm[m] = stu[0]*rm[m-1] - stu[1]*im[m-1]; 
        im[m] = stu[0]*im[m-1] + stu[1]*rm[m-1];
    }
    double rho  = mu/(req*rmag);
    double rhop = req/rmag;
    double g1 = 0.00; double g2 = 0.00; double g3 = 0.00; double g4 = 0.00;
    for(int n = 0; n <= nmax; n++) {
        double g1t = 0.0; double g2t = 0.0; double g3t = 0.0; double g4t = 0.0;
        double sm = 0.5;
        for(int m = 0; m <= n; m++) {
            double anmp1;
            if(n == m) {
                anmp1 = 0.0;
            }
            else {
                anmp1 = anm[nm2i(n,m+1)];
            }

            double dnm = cnm[nm2i(n,m)]*rm[m+1] + snm[nm2i(n,m)]*im[m+1];
            double enm = cnm[nm2i(n,m)]*rm[m] + snm[nm2i(n,m)]*im[m];
            double fnm = snm[nm2i(n,m)]*rm[m] - cnm[nm2i(n,m)]*im[m];
            double alpha  = sqrt(sm*(n-m)*(n+m+1));
            g1t += anm[nm2i(n,m)]*m*enm;
            g2t += anm[nm2i(n,m)]*m*fnm;
            g3t += alpha*anmp1*dnm;
            g4t += ((n+m+1)*anm[nm2i(n,m)]+alpha*stu[2]*anmp1)*dnm;
            // printf("ANM:   %d %d %.2e %.2e\n", n, m, anm[nm2i(n,m)], anmp1);
            // printf("DEF:   %d %d %.2e %.2e %.2e\n", n, m, dnm, enm, fnm);
            // printf("G1-4t: %d %d %.2e %.2e %.2e %.2e\n", n, m, g1t, g2t, g3t, g4t);
            // printf("CS:    %d %d %.2e %.2e\n", n, m, cnm[nm2i(n,m)], snm[nm2i(n,m)]);
            if(m == 0) sm = 1.0;
        }
        rho *= rhop;
        g1 += rho*g1t; 
        g2 += rho*g2t; 
        g3 += rho*g3t; 
        g4 += rho*g4t;
        // printf("n=%d, g1 = %.2e, g2 = %.2e, g3 = %.2e, g4 = %.2e\n", 
        // n, g1, g2, g3, g4);
    }

    Eigen::VectorXd result(3);
    result[0] = g1-g4*stu[0];
    result[1] = g2-g4*stu[1];
    result[2] = g3-g4*stu[2];
    return result;
}

Eigen::MatrixXd acceleration(Eigen::MatrixXd r_ecef, int nmax, char* model_name) {

    // if nmax < 0, raise an error
    if(nmax < 0) {
        PyErr_SetString(PyExc_ValueError, "nmax must be greater than or equal to 0");
        return Eigen::MatrixXd::Zero(1,1);
    }

    int npts = r_ecef.rows();

    set_indices(model_name, &model_index, &body_index);
    set_body_params(body_index, &mu, &req);
    int sz = nm2i(nmax+2, nmax+2);
    Eigen::VectorXd cnm(sz);
    Eigen::VectorXd snm(sz);
    read_cnm_snm(nmax, model_index, cnm, snm);

    Eigen::MatrixXd accel_vector(npts, 3);

    for(int i = 0; i < npts; i++) {
        Eigen::VectorXd result = pinesnorm(r_ecef.row(i), cnm, snm, nmax, mu, req);
        accel_vector.row(i) = result;
    }

    return accel_vector;
}
