#include "helpers.hpp"

#include "Complex.hpp"
#include "typedefs.hpp"
#include "petls.hpp"
#include "dFlag.hpp"
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <vector>
#include <string>

// #include <catch2/catch_test_macros.hpp>
// #include <catch2/catch_all.hpp>

SparseMatrixInt read_eigen_sparse(std::string filename){
    SparseMatrixDouble m;
    bool s = Eigen::loadMarket(m, filename);
    return m.cast<int>();
    }

int eigs_bdcsvd(DenseMatrix_PL &L){
    Eigen::VectorXf eigenvalues;
    auto start_eigs = std::chrono::high_resolution_clock::now();

    Eigen::BDCSVD<DenseMatrix_PL> bdcsvd(L);
    eigenvalues = bdcsvd.singularValues();
    
    auto end_eigs = std::chrono::high_resolution_clock::now();
    auto duration_eigs = std::chrono::duration_cast<std::chrono::milliseconds>(end_eigs - start_eigs);
    // petls::print_vector_precise(eigenvalues);
    return (int) duration_eigs.count();
}


int eigs_standard(DenseMatrix_PL &L){
    Eigen::VectorXf eigenvalues;

    auto start_eigs = std::chrono::high_resolution_clock::now();   

    Eigen::EigenSolver<DenseMatrix_PL> es(L); // Matrix is no longer self adjoint, must use standard Eigensolver
    eigenvalues = es.eigenvalues().real();

    auto end_eigs = std::chrono::high_resolution_clock::now();
    auto duration_eigs = std::chrono::duration_cast<std::chrono::milliseconds>(end_eigs - start_eigs);
    // petls::print_vector_precise(eigenvalues);
    return (int) duration_eigs.count();
}

int eigs_selfadjoint(DenseMatrix_PL &L){
    Eigen::VectorXf eigenvalues;
    auto start_eigs = std::chrono::high_resolution_clock::now();   

    Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es(L,Eigen::EigenvaluesOnly);
    eigenvalues = es.eigenvalues();

    auto end_eigs = std::chrono::high_resolution_clock::now();
    auto duration_eigs = std::chrono::duration_cast<std::chrono::milliseconds>(end_eigs - start_eigs);
    // petls::print_vector_precise(eigenvalues);
    return (int) duration_eigs.count();
}

int schur_inv(DenseMatrix_PL &L, DenseMatrix_PL &Schur){
     
    int m = L.rows();
    int k;
    DenseMatrix_PL change_of_basis(m,m);

    // compute null space inefficiently for testing
    // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807

    Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
    cod.compute(L);
    DenseMatrix_PL V = cod.matrixZ().transpose();
    DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    DenseMatrix_PL P = cod.colsPermutation();
    DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns

    auto start_setup = std::chrono::high_resolution_clock::now();  
    int n = PH_basis_dense.cols();
    k = m - n;
    change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;

    // Random matrix will be linearly independent with probability 1
    DenseMatrix_PL nonharmonic_basis = Eigen::MatrixXf::Random(m,k);
    change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;
    
    DenseMatrix_PL temp = change_of_basis.inverse()*L*change_of_basis; // TODO: optimize this
    Schur = temp.topLeftCorner(k,k);
    auto end_setup = std::chrono::high_resolution_clock::now();
    auto duration_setup = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
    return duration_setup.count();
}


int schur_colpiv(DenseMatrix_PL &L, DenseMatrix_PL &Schur){
     
    int m = L.rows();
    int k;
    DenseMatrix_PL change_of_basis(m,m);

    // compute null space inefficiently for testing
    // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807

    Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
    cod.compute(L);
    DenseMatrix_PL V = cod.matrixZ().transpose();
    DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    DenseMatrix_PL P = cod.colsPermutation();
    DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns

    auto start_setup = std::chrono::high_resolution_clock::now();  
    int n = PH_basis_dense.cols();
    k = m - n;
    change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;

    // Random matrix will be linearly independent with probability 1
    DenseMatrix_PL nonharmonic_basis = Eigen::MatrixXf::Random(m,k);
    change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;
    
    Eigen::ColPivHouseholderQR<Eigen::MatrixXf> solver(change_of_basis);
    DenseMatrix_PL temp = solver.solve(L)*change_of_basis;
    Schur = temp.topLeftCorner(k,k);
    auto end_setup = std::chrono::high_resolution_clock::now();
    auto duration_setup = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
    return duration_setup.count();
}



int schur_partialpiv(DenseMatrix_PL &L, DenseMatrix_PL &Schur){
     
    int m = L.rows();
    int k;
    DenseMatrix_PL change_of_basis(m,m);

    // compute null space inefficiently for testing
    // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807

    Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
    cod.compute(L);
    DenseMatrix_PL V = cod.matrixZ().transpose();
    DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    DenseMatrix_PL P = cod.colsPermutation();
    DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns

    auto start_setup = std::chrono::high_resolution_clock::now();  
    int n = PH_basis_dense.cols();
    k = m - n;
    change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;

    // Random matrix will be linearly independent with probability 1
    DenseMatrix_PL nonharmonic_basis = Eigen::MatrixXf::Random(m,k);
    change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;
    
    Eigen::PartialPivLU<Eigen::MatrixXf> solver(change_of_basis);
    DenseMatrix_PL temp = solver.solve(L)*change_of_basis;
    Schur = temp.topLeftCorner(k,k);
    auto end_setup = std::chrono::high_resolution_clock::now();
    auto duration_setup = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
    return duration_setup.count();
}

int schur_householder(DenseMatrix_PL &L, DenseMatrix_PL &Schur){
     
    int m = L.rows();
    int k;
    DenseMatrix_PL change_of_basis(m,m);

    // compute null space inefficiently for testing
    // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807

    Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
    cod.compute(L);
    DenseMatrix_PL V = cod.matrixZ().transpose();
    DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    DenseMatrix_PL P = cod.colsPermutation();
    DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns

    auto start_setup = std::chrono::high_resolution_clock::now();  
    int n = PH_basis_dense.cols();
    k = m - n;
    change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;

    // Random matrix will be linearly independent with probability 1
    DenseMatrix_PL nonharmonic_basis = Eigen::MatrixXf::Random(m,k);
    change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;
    
    Eigen::HouseholderQR<Eigen::MatrixXf> solver(change_of_basis);
    DenseMatrix_PL temp = solver.solve(L)*change_of_basis;
    Schur = temp.topLeftCorner(k,k);
    auto end_setup = std::chrono::high_resolution_clock::now();
    auto duration_setup = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
    return duration_setup.count();
}

int schur_bdcsvd(DenseMatrix_PL &L, DenseMatrix_PL &Schur){
     
    int m = L.rows();
    int k;
    DenseMatrix_PL change_of_basis(m,m);

    // compute null space inefficiently for testing
    // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807

    Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
    cod.compute(L);
    DenseMatrix_PL V = cod.matrixZ().transpose();
    DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
    DenseMatrix_PL P = cod.colsPermutation();
    DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns

    auto start_setup = std::chrono::high_resolution_clock::now();  
    int n = PH_basis_dense.cols();
    k = m - n;
    change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;

    // Random matrix will be linearly independent with probability 1
    DenseMatrix_PL nonharmonic_basis = Eigen::MatrixXf::Random(m,k);
    change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;
    
    Eigen::BDCSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> solver(change_of_basis);
    DenseMatrix_PL temp = solver.solve(L)*change_of_basis;
    Schur = temp.topLeftCorner(k,k);
    auto end_setup = std::chrono::high_resolution_clock::now();
    auto duration_setup = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
    return duration_setup.count();
}

int main(){
    int timing;
    std::vector<int> samples = {402, 950};//, 2182, 6081};//, 22577, 30000};
    // std::vector<int> samples = {22577};
    for (int i = 0; i < samples.size(); i++){
        DenseMatrix_PL L = Eigen::MatrixXf(read_eigen_sparse("data/eigenvalues/matrix_" + std::to_string(samples[i]) + "rows.mtx").cast<float>());
        std::cout << "Loaded matrix L with samples[" << i << "] = " << samples[i] << std::endl;
        std::cout << "L.rows() = " << L.rows() << std::endl;
        
        timing = eigs_bdcsvd(L);
        std::cout << "timing bdcsvd = " << timing << std::endl;

        // timing = eigs_standard(L);
        // std::cout << "timing standard = " << timing << std::endl;

        timing = eigs_selfadjoint(L);
        std::cout << "timing selfadjoint = " << timing << std::endl;

        DenseMatrix_PL Schur;
        timing =  schur_inv(L,Schur);
        std::cout << "timing setup schur inverse = " << timing << std::endl;   
        std::cout << "Schur.rows() = " << Schur.rows() << std::endl;

        timing =  schur_colpiv(L,Schur);
        std::cout << "timing setup schur colpiv = " << timing << std::endl;   

        timing =  schur_partialpiv(L,Schur);
        std::cout << "timing setup schur partial piv = " << timing << std::endl;   

        timing =  schur_householder(L,Schur);
        std::cout << "timing setup schur householder = " << timing << std::endl;   

        timing =  schur_bdcsvd(L,Schur);
        std::cout << "timing setup schur bdcsvd = " << timing << std::endl;   


        timing = eigs_standard(Schur);
        std::cout << "timing schur standard = " << timing << std::endl; 

        timing = eigs_bdcsvd(Schur);
        std::cout << "timing schur bdcsvd = " << timing << std::endl;

    }
    assert(1 == 0);
    
}