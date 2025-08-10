#include "../../include/petls_headers/eigenvalues/eigs_algorithms.hpp"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

#include "../../include/petls_headers/petls.hpp"

#include <cmath> // std::abs

namespace petls{

    bool matrix_is_diagonal(DenseMatrix_PL &M){
        float tol = 1e-4;
        for (int j = 0; j < M.cols(); j++){
            for (int i = j+1; i < M.rows(); i++){
                if (std::abs(M(i,j)) > tol){
                    // std::cout << "matrix is not diagonal, M(" << i << ", " << j << ")=" << M(i,j) << std::endl;
                    return false; 
                }
            }
        }
        return true;
    }


    spectra_vec SelfAdjointEigen(DenseMatrix_PL& L) {
        // std::cout << "eigs_algorithms.cpp line: " << __LINE__ << std::endl;
        if (matrix_is_diagonal(L)){
            spectra_vec eigs = L.diagonal();
            round_zeros(eigs, 1e-3);
            std::sort(eigs.begin(), eigs.end()); // Diagonal is unsorted
            return eigs;
        }
        Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L, Eigen::EigenvaluesOnly);
        
        spectra_vec eigs = es.eigenvalues();
        
        round_zeros(eigs, 1e-3);
        
        return eigs;
    }

    
    spectra_vec SelfAdjointEigenSparse(SparseMatrix_PL& L) {
        
        // std::cout << "eigs_algorithms.cpp line: " << __LINE__ << std::endl;
        Eigen::SelfAdjointEigenSolver<SparseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<SparseMatrix_PL>(L, Eigen::EigenvaluesOnly);
        
        // std::cout << "eigs_algorithms.cpp line: " << __LINE__ << std::endl;
        spectra_vec eigs = es.eigenvalues();
        
        // std::cout << "eigs_algorithms.cpp line: " << __LINE__ << std::endl;
        round_zeros(eigs, 1e-3);
        
        // std::cout << "eigs_algorithms.cpp line: " << __LINE__ << std::endl;
        return eigs;
    }

    std::pair<spectra_vec,DenseMatrix_spectra_PL> SelfAdjointEigenpairsEigen(DenseMatrix_PL &L){
        Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L);
        spectra_vec eigs = es.eigenvalues();
        DenseMatrix_spectra_PL eigvs = es.eigenvectors();
        round_zeros(eigs, 1e-3);
        return std::pair<spectra_vec,DenseMatrix_spectra_PL>(eigs, eigvs);
    }

    spectra_vec EigensolverEigen(DenseMatrix_PL &L){
        Eigen::EigenSolver<DenseMatrix_PL> es = Eigen::EigenSolver<DenseMatrix_PL>(L, false); // false for compte eigenvalues only
        spectra_vec eigs = es.eigenvalues().real();
        round_zeros(eigs, 1e-3);
        std::sort(eigs.begin(), eigs.end()); // Eigen::Eigensolver returns unsorted
        return eigs;
    }

    spectra_vec BDCSVDEigen(DenseMatrix_PL &L){
        Eigen::BDCSVD<DenseMatrix_PL> bdcsvd(L.rows(), L.cols());
        bdcsvd.compute(L);
        spectra_vec eigs = bdcsvd.singularValues();
        round_zeros(eigs, 1e-3);
        std::sort(eigs.begin(), eigs.end());
      return eigs;
  }

}
