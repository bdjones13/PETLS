#include "../../include/petls_headers/eigenvalues/spectra_algorithms.hpp"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

#include "../../include/petls_headers/petls.hpp"
#include <iostream>
#include <cmath> // std::abs

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/GenEigsSolver.h>

namespace petls{

    spectra_vec DenseSymSpectra(DenseMatrix_PL &L){
         if (petls::matrix_is_diagonal(L)){
            spectra_vec eigs = L.diagonal();
            round_zeros(eigs, 1e-3);
            std::sort(eigs.begin(), eigs.end()); // Diagonal is unsorted
            return eigs;
        }
        Spectra::DenseSymMatProd<coefficient_type> op(L);

        // Construct eigen solver object, requesting the largest three eigenvalues
        int nev = 10; // Number of eigenvalues to compute
        int ncv = 20; // Dimension of the Krylov subspace, should be larger than nev

        nev = std::min(nev, (int) L.rows()-1); // algorithm requires  1 < nev <= M.rows()-1
        ncv = std::min(std::max(ncv, 2*nev), (int) L.rows());// Try to make ncv >= 2*nev, but not larger than the size of the matrix


        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<coefficient_type>> eigensolver(op, nev, ncv);

        // Initialize and compute
        eigensolver.init();
        int nconv = eigensolver.compute(Spectra::SortRule::SmallestMagn);
        if (eigensolver.info() != Spectra::CompInfo::Successful) {
            std::cout << "Eigenvalue computation failed" << std::endl;
            std::cout << "nev="<<nev << ", ncv="<<ncv << ", nconv="<<nconv << ", rows=" <<L.rows() << std::endl;
            // Fall back to Eigen if Spectra fails
            Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L, Eigen::EigenvaluesOnly);
        
            spectra_vec eigs = es.eigenvalues();
            
            round_zeros(eigs, 1e-3);
            
            return eigs;
            // throw std::runtime_error("Eigenvalue computation failed using Spectra. This can happen if the matrix has a large null-space (high betti-number).");
        }
        spectra_vec eigs = eigensolver.eigenvalues();
        round_zeros(eigs, 1e-3);
        return eigs;
    }

  

    // std::pair<spectra_vec,DenseMatrix_spectra_PL> SelfAdjointEigenpairsEigen(DenseMatrix_PL &L){
    //     Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L);
    //     spectra_vec eigs = es.eigenvalues();
    //     DenseMatrix_spectra_PL eigvs = es.eigenvectors();
    //     round_zeros(eigs, 1e-3);
    //     return std::pair<spectra_vec,DenseMatrix_spectra_PL>(eigs, eigvs);
    // }



}
