
#ifndef HELPERS
#define HELPERS

#include <vector>
#include <iostream>
#include <string>
// #include "pch.hpp"
#include "Complex.hpp"
#include "petls.hpp"



void print_std_vec(std::vector<float> vec){
    std::cout << "[";
    for (int i = 0; i < (int) vec.size(); i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << "]" << std::endl;
}

void print_std_vec(std::vector<double> vec){
    std::cout << "[";
    for (int i = 0; i < (int) vec.size(); i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << "]" << std::endl;
}

void print_all_spectra(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra){
    for (int i = 0; i < (int) spectra.size(); i++ ){
        std::cout << "dim = " << std::get<0>(spectra[i]) << ", a = " << std::get<1>(spectra[i]) << ", b = " << std::get<2>(spectra[i]);
        std::cout << ", eigs = ";
        print_std_vec(std::get<3>(spectra[i]));
        std::cout << std::endl;
    }
}

petls::Complex filtered_triangle_complex_alg(){
    //Boundary matrices corresponding to directed flag laplacian of
    // dim 0:
    // 0 1 2
    // dim 1:
    // 0 1 3
    // 1 2 4
    // 0 2 5 
    SparseMatrixInt d1(3,3);
    SparseMatrixInt d2(3,1);
    d1.coeffRef(0,0) = -1;
    d1.coeffRef(1,0) = 1;
    d1.coeffRef(1,1) = -1;
    d1.coeffRef(2,1) = 1;
    d1.coeffRef(0,2) = -1;
    d1.coeffRef(2,2) = 1;

    d2.coeffRef(0,0) = 1;
    d2.coeffRef(1,0) = 1;
    d2.coeffRef(2,0) = -1;

    std::vector<filtration_type> c0_filtrations = {0, 1, 2};
    std::vector<filtration_type> c1_filtrations = {3, 4, 5};
    std::vector<filtration_type> c2_filtrations = {5};

    std::vector<SparseMatrixInt> boundaries;
    boundaries.push_back(d1);
    boundaries.push_back(d2);
    std::vector<std::vector<filtration_type>> filtrations;
    filtrations.push_back(c0_filtrations);
    filtrations.push_back(c1_filtrations);
    filtrations.push_back(c2_filtrations);
    
    petls::Complex complex(boundaries,filtrations);
    return complex;
}

bool test_sample(petls::Complex &complex, std::vector<spectra_type> reference_spectra, int dim, double a, double b){
    spectra_vec ref_spectra_eigen = Eigen::Map<spectra_vec, Eigen::Unaligned>(reference_spectra.data(),reference_spectra.size());
    std::vector<spectra_type> spectra_std_vec = complex.spectra(dim,a,b);
    spectra_vec spectra = Eigen::Map<spectra_vec, Eigen::Unaligned>(spectra_std_vec.data(),spectra_std_vec.size());
    float tol = 1e-1;

    if (spectra.size() != ref_spectra_eigen.size()){
        std::cout << "Fail compare spectra, different size vectors. a = " << a << ", b = " << b << ", dim = " << dim <<  std::endl;
        std::cout << "spectra.size() = " << spectra.size() << ", ref_spectra_eigen.size() = " << ref_spectra_eigen.size() << std::endl;
        std::cout << "Expected:";
        petls::print_vector_precise(ref_spectra_eigen);
        std::cout << "Actual:";
        petls::print_vector_precise(spectra);
        return false;
    }

    if(ref_spectra_eigen.isApprox(spectra, tol)){
        return true;
    } else{
        std::string message = "Fail spectra a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", dim=" + std::to_string(dim);
        std::cout << message << std::endl;
        std::cout << "Expected:";
        petls::print_vector_precise(ref_spectra_eigen);
        std::cout << "Actual:";
        petls::print_vector_precise(spectra);
        return false;
    }
}

bool test_sample(petls::Complex &complex, std::vector<spectra_type> reference_spectra, DenseMatrix_PL reference_eigenvectors, int dim, double a, double b){
    spectra_vec ref_spectra_eigen = Eigen::Map<spectra_vec, Eigen::Unaligned>(reference_spectra.data(),reference_spectra.size());

    std::pair<std::vector<spectra_type>,DenseMatrix_PL> eigenpairs = complex.eigenpairs(dim, a, b);
    std::vector<spectra_type> eigs = eigenpairs.first;
    DenseMatrix_PL eigvs = eigenpairs.second;

    spectra_vec spectra = Eigen::Map<spectra_vec, Eigen::Unaligned>(eigs.data(),eigs.size());
    spectra_type tol = 1e-1;
    if (spectra.size() != ref_spectra_eigen.size()){
        std::cout << "Fail compare spectra, different size vectors. a = " << a << ", b = " << b << ", dim = " << dim <<  std::endl;
        std::cout << "spectra.size() = " << spectra.size() << ", ref_spectra_eigen.size() = " << ref_spectra_eigen.size() << std::endl;
        return false;
    }
    // fail if eigenvalues wrong
    if(!ref_spectra_eigen.isApprox(spectra, tol)){
        std::string message = "Fail spectra basic triangle a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", dim=" + std::to_string(dim);
        std::cout << message << std::endl;
        std::cout << "Expected:";
        // petls::print_vector_precise(ref_spectra_eigen);
        std::cout << "Actual:";
        // petls::print_vector_precise(spectra);
        return false;
    }
    // fail if eigenvectors wrong
    if (!eigvs.isApprox(reference_eigenvectors,tol)){
        std::string message = "Fail eigenvectors basic triangle a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", dim=" + std::to_string(dim);
        std::cout << message << std::endl;
        std::cout << "Expected:";
        std::cout << reference_eigenvectors << std::endl;
        std::cout << "Actual:";
        std::cout << eigvs << std::endl;
        return false;
    }
    return true;
}


// template<class up_alg>
// bool test_sample_vec(petls::Complex<up_alg> &complex, std::vector<float> reference_spectra, int dim, double a, double b){
//     Eigen::VectorXf ref_spectra_eigen = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(reference_spectra.data(),reference_spectra.size());
//     std::vector<float> spectra_std_vec = complex.spectra(dim,a,b);
//     Eigen::VectorXf spectra = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(spectra_std_vec.data(),spectra_std_vec.size());
//     float tol = 1e-1;
//     if(ref_spectra_eigen.isApprox(spectra, tol)){
//         return true;
//     } else{
//         std::string message = "Fail spectra basic triangle a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", dim=" + std::to_string(dim);
//         std::cout << message << std::endl;
//         std::cout << "Expected:";
//         // petls::print_vector_precise(ref_spectra_eigen);
//         std::cout << "Actual:";
//         // petls::print_vector_precise(spectra);
//         return false;
//     }
// }

petls::Complex filtered_triangle_complex(){
    //Boundary matrices corresponding to directed flag laplacian of
    // dim 0:
    // 0 1 2
    // dim 1:
    // 0 1 3
    // 1 2 4
    // 0 2 5 
    SparseMatrixInt d1(3,3);
    SparseMatrixInt d2(3,1);
    d1.coeffRef(0,0) = -1;
    d1.coeffRef(1,0) = 1;
    d1.coeffRef(1,1) = -1;
    d1.coeffRef(2,1) = 1;
    d1.coeffRef(0,2) = -1;
    d1.coeffRef(2,2) = 1;

    d2.coeffRef(0,0) = 1;
    d2.coeffRef(1,0) = 1;
    d2.coeffRef(2,0) = -1;

    std::vector<filtration_type> c0_filtrations = {0, 1, 2};
    std::vector<filtration_type> c1_filtrations = {3, 4, 5};
    std::vector<filtration_type> c2_filtrations = {5};

    std::vector<SparseMatrixInt> boundaries;
    boundaries.push_back(d1);
    boundaries.push_back(d2);
    std::vector<std::vector<filtration_type>> filtrations;
    filtrations.push_back(c0_filtrations);
    filtrations.push_back(c1_filtrations);
    filtrations.push_back(c2_filtrations);
    
    petls::Complex complex(boundaries,filtrations);
    return complex;
}

bool compare_spectra(std::vector<float> observed_spectra, std::vector<float> reference_spectra, float tol){
    // convert to eigen to use isApprox
    Eigen::VectorXf ref_spectra_eigen = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(reference_spectra.data(),reference_spectra.size());
    Eigen::VectorXf observed_spectra_eigen = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(observed_spectra.data(),observed_spectra.size());
    return observed_spectra_eigen.isApprox(ref_spectra_eigen, tol);
}
#endif