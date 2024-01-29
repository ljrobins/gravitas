#include <Python.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>  // If you're using Eigen types

namespace py = pybind11;
using namespace std;

#include "libgrav.hpp"

PYBIND11_MODULE(gravitas, m) {
    m.def("earth_acceleration", &earth_acceleration, R"mydelimiter(
        Acceleration due to the Earth's gravity field, as defined by the EGM96 model.

        :param r_ecef_km: Position in ECEF coordinates, in km.
        :type r_ecef_km: ndarray [nx3]
        :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
        :type nmax: int
        :return: Acceleration in ECI coordinates, in km/s^2.
        :rtype: float
        )mydelimiter");
    m.def("moon_acceleration", &moon_acceleration, R"mydelimiter(
        Acceleration due to the Moon's gravity field, as defined by the GRGM360 model.

        :param r_mcmf_km: Position in MCMF coordinates, in km.
        :type r_mcmf_km: ndarray [nx3]
        :param nmax: Maximum degree of the spherical harmonic expansion, 0 <= nmax <= 360.
        :type nmax: int
        :return: Acceleration in MCI coordinates, in km/s^2.
        :rtype: float
        )mydelimiter");
}