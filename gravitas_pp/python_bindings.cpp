#include <Python.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>  // If you're using Eigen types

namespace py = pybind11;
using namespace std;

#include "libgrav.hpp"

PYBIND11_MODULE(gravitas_cpp, m) {
    m.doc() = "gravitas module";  // Optional module docstring
    m.def("print_test", &print_test, "A function that prints the first value of n_coef_EGM96");
    m.def("acceleration", &acceleration, "A function that returns the acceleration of a body due to a gravity model");
}