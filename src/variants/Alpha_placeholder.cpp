// // #include "../cpp/include/petls_headers/petls.hpp"

// #include <pybind11/stl.h>

// #include <pybind11/pybind11.h>
// #include "Complex.hpp"
// // #include "Simplex.hpp"
// #include <pybind11/eigen.h> //https://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#Using-the-C++-eigen-library-to-calculate-matrix-inverse-and-determinant

// namespace py = pybind11;

// // Give dummy class and constructor for when alpha complex not installed
// void init_Alpha(py::module &m) {
//  py::class_<petls::Simplex>(m, "Alpha")
//      .def(py::init<>());
// }