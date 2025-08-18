// #include "Complex.hpp"
// #include "typedefs.hpp"
// #include "petls.hpp"
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <chrono>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>


using SparseMatrixInt = Eigen::SparseMatrix<int>;
using SparseMatrixDouble = Eigen::SparseMatrix<double>;
using DenseMatrix_PL = Eigen::MatrixXf;




Eigen::VectorXf DenseSymSpectra(DenseMatrix_PL &L,Spectra::SortRule sort_rule){
        Spectra::DenseSymMatProd<float> op(L);

        // Construct eigen solver object, requesting the largest three eigenvalues
        int nev = 10; // Number of eigenvalues to compute
        int ncv = 20; // Dimension of the Krylov subspace, should be larger than nev

        nev = std::min(nev, (int) L.rows()-1); // algorithm requires  1 < nev <= M.rows()-1
        ncv = std::min(std::max(ncv, 2*nev), (int) L.rows());// Try to make ncv >= 2*nev, but not larger than the size of the matrix
        if (nev < 1){
            Eigen::VectorXf eigs(1);
            eigs(0) = L(0,0);
            return eigs;
        }

        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<float>> eigensolver(op, nev, ncv);

        // Initialize and compute
        eigensolver.init();
        int nconv = eigensolver.compute(sort_rule);
        if (eigensolver.info() != Spectra::CompInfo::Successful) {
            std::cout << "Eigenvalue computation failed" << std::endl;
            std::cout << "nev="<<nev << ", ncv="<<ncv << ", nconv="<<nconv << ", rows=" <<L.rows() << std::endl;
            // Fall back to Eigen if Spectra fails
            Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L, Eigen::EigenvaluesOnly);
        
            Eigen::VectorXf eigs = es.eigenvalues();
            
            return eigs;
            // throw std::runtime_error("Eigenvalue computation failed using Spectra. This can happen if the matrix has a large null-space (high betti-number).");
        }
        Eigen::VectorXf eigs = eigensolver.eigenvalues();
        return eigs;
}

Eigen::VectorXf InverseSpectra(DenseMatrix_PL &L){
    Spectra::DenseSymShiftSolve<float> op(L);

    int nev = 10; // Number of eigenvalues to compute
    int ncv = 20; // Dimension of the Krylov subspace, should be larger than nev

    nev = std::min(nev, (int) L.rows()-1); // algorithm requires  1 < nev <= M.rows()-1
    ncv = std::min(std::max(ncv, 2*nev), (int) L.rows());// Try to make ncv >= 2*nev, but not larger than the size of the matrix
    if (nev < 1){
        Eigen::VectorXf eigs(1);
        eigs(0) = L(0,0);
        return eigs;
    }
    Spectra::SymEigsShiftSolver<Spectra::DenseSymShiftSolve<float>> eigensolver(op, nev, ncv, 0.01);
    eigensolver.init();
    int nconv = eigensolver.compute(Spectra::SortRule::LargestMagn);
    if (eigensolver.info() != Spectra::CompInfo::Successful){
        std::cout << "Eigenvalue computation failed" << std::endl;
        std::cout << "nev="<<nev << ", ncv="<<ncv << ", nconv="<<nconv << ", rows=" <<L.rows() << std::endl;
        // Fall back to Eigen if Spectra fails
        Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L, Eigen::EigenvaluesOnly);
    
        Eigen::VectorXf eigs = es.eigenvalues();
        
        return eigs;
        // throw std::runtime_error("Eigenvalue computation failed using Spectra. This can happen if the matrix has a large null-space (high betti-number).");
    }
    Eigen::VectorXf eigs = eigensolver.eigenvalues();
    return eigs;
}


// std::vector<double> get_filtrations(double delta){
//     double start = 0.0;
//     double end = 5.0;
//     double current = start;
//     std::vector<double> filtrations;
//     while (current <= end){
//         filtrations.push_back(current);
//         current += delta;
//     }
//     return filtrations;
// }

int main(int argc, char *argv[]){
    if (argc < 3){
        std::cout << "Usage: " << argv[0] << " <algorithm> <replicate>" << std::endl;
    }
    // double delta = 0.1;
    int replicate = std::stoi(argv[2]);
    std::string algorithm = argv[1];

    std::srand(replicate);

    // std::vector<double> filtrations = get_filtrations(delta);
    // std::vector<int> dims = {0, 1, 2};
    std::vector<int> sizes = {10, 20};//, 50, 100, 200, 500, 1000, 2000, 5000, 10000};
    int duration;
    Eigen::VectorXf eigenvalues;
    // output file
    std::ofstream out("./random_profiles/" + algorithm + "/r" + std::to_string(replicate) + ".csv");
    out << "simplices," << algorithm << std::endl;
    for (int size : sizes){
       
        Eigen::VectorXf eigenvalues = Eigen::VectorXf::Random(size);
        
        std::chrono::_V2::steady_clock::time_point start_t = std::chrono::steady_clock::now();

        if (algorithm == "eigen.selfadjoint"){
            Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es(L,Eigen::EigenvaluesOnly);
            eigenvalues = es.eigenvalues();
        } 
        else if (algorithm == "eigen.bdcsvd"){
            Eigen::BDCSVD<DenseMatrix_PL> bdcsvd(L);
            eigenvalues = bdcsvd.singularValues();
        }
        else if (algorithm == "eigen.eigensolver"){
            Eigen::EigenSolver<DenseMatrix_PL> es(L);
            Eigen::VectorXcf eigenvalues_complex = es.eigenvalues();
            eigenvalues = eigenvalues_complex.real();
            std::sort(eigenvalues.begin(), eigenvalues.end()); // Eigen::Eigensolver returns unsorted
        }
        else if (algorithm == "spectra.dense.largest"){
            eigenvalues = DenseSymSpectra(L, Spectra::SortRule::LargestMagn);
        }
        else if (algorithm == "spectra.dense.smallest"){
            eigenvalues = DenseSymSpectra(L, Spectra::SortRule::SmallestMagn);
        }
        else if (algorithm == "spectra.inverse.smallest"){
            eigenvalues = InverseSpectra(L);
        }
        auto end_t = std::chrono::steady_clock::now();
        duration = (std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t)).count();
        out << dim << "," << filtration << "," << filtration + delta << "," 
            << simplices << "," << duration << std::endl;
        std::cout << "Duration: " << duration << " ms" << std::endl;
        
    }
    out.close();
    // int timing;
    // std::vector<int> samples = {402, 950};//, 2182, 6081};//, 22577, 30000};
    // // std::vector<int> samples = {22577};
    // for (int i = 0; i < samples.size(); i++){
    //     DenseMatrix_PL L = Eigen::MatrixXf(read_eigen_sparse("data/eigenvalues/matrix_" + std::to_string(samples[i]) + "rows.mtx").cast<float>());
    //     std::cout << "Loaded matrix L with samples[" << i << "] = " << samples[i] << std::endl;
    //     std::cout << "L.rows() = " << L.rows() << std::endl;
        
    //     timing = eigs_bdcsvd(L);
    //     std::cout << "timing bdcsvd = " << timing << std::endl;

    //     // timing = eigs_standard(L);
    //     // std::cout << "timing standard = " << timing << std::endl;

    //     timing = eigs_selfadjoint(L);
    //     std::cout << "timing selfadjoint = " << timing << std::endl;

    //     DenseMatrix_PL Schur;
    //     timing =  schur_inv(L,Schur);
    //     std::cout << "timing setup schur inverse = " << timing << std::endl;   
    //     std::cout << "Schur.rows() = " << Schur.rows() << std::endl;

    //     timing =  schur_colpiv(L,Schur);
    //     std::cout << "timing setup schur colpiv = " << timing << std::endl;   

    //     timing =  schur_partialpiv(L,Schur);
    //     std::cout << "timing setup schur partial piv = " << timing << std::endl;   

    //     timing =  schur_householder(L,Schur);
    //     std::cout << "timing setup schur householder = " << timing << std::endl;   

    //     timing =  schur_bdcsvd(L,Schur);
    //     std::cout << "timing setup schur bdcsvd = " << timing << std::endl;   


    //     timing = eigs_standard(Schur);
    //     std::cout << "timing schur standard = " << timing << std::endl; 

    //     timing = eigs_bdcsvd(Schur);
    //     std::cout << "timing schur bdcsvd = " << timing << std::endl;

    // }
    // assert(1 == 0);
    
}