#include "Complex.hpp"
#include "typedefs.hpp"
#include "petls.hpp"
#include "dFlag.hpp"
#include <iostream>
#include "../helpers.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/SparseExtra>

void test_rand(int size){
    Eigen::MatrixXf m = Eigen::MatrixXf::Random(size,size);
    Eigen::MatrixXf a = m.transpose()*m;
    auto start_spectra = std::chrono::high_resolution_clock::now();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(a, Eigen::EigenvaluesOnly);
    es.eigenvalues();
    auto end_spectra = std::chrono::high_resolution_clock::now();
    auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
    std::cout << "duration rand eigs: " << duration_spectra.count() << std::endl;
    assert(1 == 0);
    return;
}