#ifndef typedefs_h
#define typedefs_h

#include "Eigen/SparseCore" //also includes core





typedef float coefficient_type;
typedef float spectra_type;

typedef Eigen::SparseMatrix<coefficient_type, Eigen::ColMajor> SparseMatrix_PL;
typedef Eigen::Matrix<coefficient_type, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix_PL;
typedef Eigen::Matrix<spectra_type, Eigen::Dynamic, Eigen::Dynamic> SparseMatrix_spectra_PL;
typedef Eigen::Matrix<spectra_type, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix_spectra_PL;
typedef Eigen::Matrix<spectra_type, Eigen::Dynamic, 1> spectra_vec;

typedef Eigen::SparseMatrix<int, Eigen::ColMajor> SparseMatrixInt;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrixDouble;
typedef Eigen::SparseMatrix<float, Eigen::ColMajor> SparseMatrixFloat;


typedef double filtration_type; // TODO: single precision
typedef int index_type; // TODO: change from int to uint



//TODO: remove PRUNE_CONSTANT or determine a justified value
#define PRUNE_CONSTANT 0.0000001


#endif