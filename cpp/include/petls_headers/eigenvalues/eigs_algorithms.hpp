#ifndef eigs_algs_H
#define eigs_algs_H

#include "../typedefs.hpp"


namespace petls{

    spectra_vec SelfAdjointEigen(DenseMatrix_PL& L) ;
    spectra_vec SelfAdjointEigenSparse(SparseMatrix_PL& L);
    std::pair<spectra_vec,DenseMatrix_spectra_PL> SelfAdjointEigenpairsEigen(DenseMatrix_PL &L);
    spectra_vec EigensolverEigen(DenseMatrix_PL &L);
    spectra_vec BDCSVDEigen(DenseMatrix_PL &L);
}

#endif