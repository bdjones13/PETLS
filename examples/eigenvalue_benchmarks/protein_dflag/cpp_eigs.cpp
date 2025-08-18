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


bool DenseSymSpectra(DenseMatrix_PL &L,Spectra::SortRule sort_rule){
        Spectra::DenseSymMatProd<float> op(L);

        // Construct eigen solver object, requesting the largest three eigenvalues
        int nev = 10; // Number of eigenvalues to compute
        int ncv = 20; // Dimension of the Krylov subspace, should be larger than nev

        nev = std::min(nev, (int) L.rows()-1); // algorithm requires  1 < nev <= M.rows()-1
        ncv = std::min(std::max(ncv, 2*nev), (int) L.rows());// Try to make ncv >= 2*nev, but not larger than the size of the matrix
        if (nev < 1){
            Eigen::VectorXf eigs(1);
            eigs(0) = L(0,0);
            return true;
            // return eigs;
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
            return false;//failed
            // return eigs;
            // throw std::runtime_error("Eigenvalue computation failed using Spectra. This can happen if the matrix has a large null-space (high betti-number).");
        }
        Eigen::VectorXf eigs = eigensolver.eigenvalues();
        std::cout << "Eigs=" << eigs << std::endl;
        return true;//succeeded
}

bool InverseSpectra(DenseMatrix_PL &L){
    Spectra::DenseSymShiftSolve<float> op(L);

    int nev = 10; // Number of eigenvalues to compute
    int ncv = 20; // Dimension of the Krylov subspace, should be larger than nev

    nev = std::min(nev, (int) L.rows()-1); // algorithm requires  1 < nev <= M.rows()-1
    ncv = std::min(std::max(ncv, 2*nev), (int) L.rows());// Try to make ncv >= 2*nev, but not larger than the size of the matrix
    if (nev < 1){
        Eigen::VectorXf eigs(1);
        eigs(0) = L(0,0);
        return true;//succeeded
        // return eigs;
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
        
        return false;//failed
        // return eigs;
        // throw std::runtime_error("Eigenvalue computation failed using Spectra. This can happen if the matrix has a large null-space (high betti-number).");
    }
    Eigen::VectorXf eigs = eigensolver.eigenvalues();
    // std::cout << "Eigs=" << eigs << std::endl;
    // return eigs;
    return true;//succeeded
}


std::vector<double> get_filtrations(double delta){
    double start = 0.0;
    double end = 10.0;
    double current = start;
    std::vector<double> filtrations;
    while (current <= end){
        filtrations.push_back(current);
        current += delta;
    }
    return filtrations;
}

int main(int argc, char *argv[]){
    if (argc < 3){
        std::cout << "Usage: " << argv[0] << " <algorithm> <replicate>" << std::endl;
    }
    double delta = 0.1;
    int replicate = std::stoi(argv[2]);
    std::string algorithm = argv[1];
    std::vector<double> filtrations = get_filtrations(delta);
    // std::vector<int> dims = {0, 1, 2};
    std::vector<int> dims = {1, 2};
    int duration;
    Eigen::VectorXf eigenvalues;
    // output file
    std::ofstream out("./scratch_profiles/" + algorithm + "/r" + std::to_string(replicate) + ".csv");
    out << "dimension,a,b,simplices," << algorithm << "," << algorithm << "_failed" << std::endl;
    for (double filtration : filtrations){
        for (int dim : dims){
            
            bool success = true;
            
            std::cout << "Filtration: " << filtration << ", Dimension: " << dim;
            
            // open matrix file and read matrix
            std::stringstream stream;
            stream << "./protein_L/dim" << dim << "_a" << std::fixed << std::setprecision(2)
                << filtration << "_b" << std::fixed << std::setprecision(2) << filtration+delta << ".mkt";
            std::string filename = stream.str();
            SparseMatrixDouble Ld;
            bool s = Eigen::loadMarket(Ld, filename);
            if (!s){
                std::cout << "Failed to load matrix from file: " << filename << std::endl;
                out << dim << "," << filtration << "," << filtration + delta << "," << 0 << "," << 0 << std::endl;
                continue;
            }
            Eigen::MatrixXf L_sparse = Ld.cast<float>();
            DenseMatrix_PL L = Eigen::MatrixXf(L_sparse);
            // empty case
            int simplices = L_sparse.rows();
            if (simplices == 0){
                std::cout << "No simplices found in file: " << filename << std::endl;
                out << dim << "," << filtration << "," << filtration + delta << "," << 0 << "," << 0 << ",False" << std::endl;
                continue;
            }
            
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
                success = DenseSymSpectra(L, Spectra::SortRule::LargestMagn);
            }
            // else if (algorithm == "spectra.dense.smallest"){
            //     eigenvalues = DenseSymSpectra(L, Spectra::SortRule::SmallestMagn);
            // }
            else if (algorithm == "spectra.inverse.smallest"){
                success = InverseSpectra(L);
            }
            auto end_t = std::chrono::steady_clock::now();
            duration = (std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t)).count();
            out << dim << "," << filtration << "," << filtration + delta << "," 
                << simplices << "," << duration;
                
            if (success){
                out << ",False" << std::endl; // the sheet will have "True" whenever it failed, because we want to count the failures
            } else{
                out << ",True" << std::endl;
                success = true;
            }
            std::cout << "Duration: " << duration << " ms" << std::endl;
        }
    }
    out.close();
    std::cout << "Finished processing all filtrations and dimensions." << std::endl;
    return 0;
}