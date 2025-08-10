#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
#include "Rips.hpp"
#include <pybind11/eigen.h> //https://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#Using-the-C++-eigen-library-to-calculate-matrix-inverse-and-determinant

namespace py = pybind11;

void define_rips(py::module &m){
    py::class_<petls::Rips, petls::Complex>(m,"Rips")
     .def(py::init<const char* , int>()) // Rips(const char* filename, int max_dim);
     .def(py::init<const char* , int, float>()) // Rips(const char* filename, int max_dim, float threshold);
     .def_static("from_distances",
          [](DenseMatrix_PL distances, int max_dim) {
               return petls::Rips::from_distances(distances, max_dim);
          })
     .def_static("from_distances",
          [](DenseMatrix_PL distances, int max_dim, float threshold) {
               return petls::Rips::from_distances(distances, max_dim, threshold);
          })
     .def(py::init<std::vector<std::vector<double>>, int>()) // Rips(std::vector<std::vector<double>> points, int max_dim);
     .def(py::init<std::vector<std::vector<double>>, int, float>()); // Rips(std::vector<std::vector<double>> points, int max_dim, float threshold);
}


void init_Rips(py::module &m) {
     define_rips(m);
}