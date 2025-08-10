// #include "../cpp/include/petls_headers/petls.hpp"

#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
// #include "Complex.hpp"
#include "dFlag.hpp"

#include <pybind11/eigen.h> //https://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#Using-the-C++-eigen-library-to-calculate-matrix-inverse-and-determinant

namespace py = pybind11;

void define_dflag(py::module &m){
     py::class_<petls::dFlag, petls::Complex>(m,"dFlag")
     .def(py::init<const char* , int>()); // Alpha(const char* filename, int max_dim); // .OFF file 
     
}

void init_dFlag(py::module &m) {
     define_dflag(m);
     // define_dflag<petls::schur<int>>(m,"eigensolver","schur");
     // define_dflag<petls::schur<int>>(m,"bdcsvd","schur");
}