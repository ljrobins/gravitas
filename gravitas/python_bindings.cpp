#include <Python.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>  // If you're using Eigen types

namespace py = pybind11;
using namespace std;

#include "libgrav.hpp"

PYBIND11_MODULE(gravitas, m) {
    m.doc() = "gravitas module";  // Optional module docstring
    m.def("earth_acceleration", &earth_acceleration, "A function that returns the acceleration of a body due to a gravity model");
    m.def("moon_acceleration", &moon_acceleration, "A function that returns the acceleration of a body due to a gravity model");
}