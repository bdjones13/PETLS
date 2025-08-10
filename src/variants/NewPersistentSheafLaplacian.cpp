#include <pybind11/stl.h>

#include <pybind11/pybind11.h>
#include "newPSL.hpp"
#include <pybind11/eigen.h> //https://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#Using-the-C++-eigen-library-to-calculate-matrix-inverse-and-determinant
#include <string>
#include <pybind11/functional.h>
namespace py = pybind11;
typedef float spectra_type;

void define_psl(py::module &m){
    using storage = float;
    // using pl = petls::NewPersistentSheafLaplacian;
    // Each python class will be Complex_EigsAlgorithm_UpAlgorithm
    // e.g. Complex_selfadjoint_schur
    py::class_<petls::NewPersistentSheafLaplacian>(m,"NewPersistentSheafLaplacian")                                                  
    .def(py::init< std::vector<Eigen::SparseMatrix<storage>>, std::vector<std::vector<double>> >()) // NewPersistentSheafLaplacian(std::vector<SparseMatrixInt> boundaries,
                                                                                    //    std::vector<std::vector<filtration_type>> filtrations);
    
    .def("set_eigs_algorithm_func",py::overload_cast<std::string>(&petls::NewPersistentSheafLaplacian::set_eigs_algorithm_func), "Set eigenvalue algorithm by name")
    .def("set_eigs_algorithm_func",py::overload_cast<std::function<spectra_vec(DenseMatrix_PL&)>>(&petls::NewPersistentSheafLaplacian::set_eigs_algorithm_func), "Set eigenvalue algorithm to the function pointer")
    .def("set_boundaries_filtrations",
        py::overload_cast<std::vector<Eigen::SparseMatrix<storage>>, std::vector<std::vector<double>> >( &petls::NewPersistentSheafLaplacian::set_boundaries_filtrations), "Set the boundary matrices and filtrations of the Persistent Laplacian")
    .def("set_verbose",py::overload_cast<bool>(&petls::NewPersistentSheafLaplacian::set_verbose))
    .def("set_flipped",py::overload_cast<bool>(&petls::NewPersistentSheafLaplacian::set_flipped))
    .def("get_L", py::overload_cast<int, double, double>(&petls::NewPersistentSheafLaplacian::get_L))
    .def("get_up", py::overload_cast<int, double, double>(&petls::NewPersistentSheafLaplacian::get_up))
    .def("get_down", py::overload_cast<int, double>(&petls::NewPersistentSheafLaplacian::get_down))
    .def("nonzero_spectra",py::overload_cast<int,double,double,Eigen::SparseMatrix<spectra_type>,bool>(&petls::NewPersistentSheafLaplacian::nonzero_spectra))
    .def("spectra", py::overload_cast< std::vector<std::tuple<int,double,double>> > (&petls::NewPersistentSheafLaplacian::spectra))
    .def("spectra", py::overload_cast<int, double, double>                          (&petls::NewPersistentSheafLaplacian::spectra))
    .def("spectra", py::overload_cast<>                                             (&petls::NewPersistentSheafLaplacian::spectra))
    .def("eigenpairs", py::overload_cast< std::vector<std::tuple<int,double,double>> > (&petls::NewPersistentSheafLaplacian::eigenpairs))
    .def("eigenpairs", py::overload_cast<int, double, double>                          (&petls::NewPersistentSheafLaplacian::eigenpairs))
    .def("eigenpairs", py::overload_cast<>                                             (&petls::NewPersistentSheafLaplacian::eigenpairs))
    .def("eigenvalues_summarize", py::overload_cast<std::vector<spectra_type>>( &petls::NewPersistentSheafLaplacian::eigenvalues_summarize))
    .def("store_L", py::overload_cast<int, double, double, std::string>(&petls::NewPersistentSheafLaplacian::store_L))
    .def("print_boundaries", py::overload_cast<>( &petls::NewPersistentSheafLaplacian::print_boundaries))
    .def("store_spectra",py::overload_cast<std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>>,std::string>(&petls::NewPersistentSheafLaplacian::store_spectra))
    .def("store_spectra_summary",py::overload_cast<std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>>,std::string>(&petls::NewPersistentSheafLaplacian::store_spectra_summary))
    .def("filtration_list_to_spectra_request",py::overload_cast<std::vector<filtration_type>, std::vector<int>>(&petls::NewPersistentSheafLaplacian::filtration_list_to_spectra_request))
    .def("get_all_filtrations",py::overload_cast<>(&petls::NewPersistentSheafLaplacian::get_all_filtrations));
}

void init_PersistentSheafLaplacian(py::module &m) {
    // Call the above function multiple times, one for each pair of eigensolver and up algorithm combination
    define_psl(m);
    // define_pl<int, petls::schur<int>>(m,"schur","");
    // define_pl<int, petls::schur<int>>(m,"schur","");

    // define_pl<float>(m,"_float");
    // define_pl<float, petls::schur<float>>(m,"schur","_float");
    // define_pl<float, petls::schur<float>>(m,"schur","_float");
}