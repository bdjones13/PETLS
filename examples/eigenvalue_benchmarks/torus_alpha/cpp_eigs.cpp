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



// SparseMatrixInt read_eigen_sparse(std::string filename){
//     SparseMatrixDouble m;
//     bool s = Eigen::loadMarket(m, filename);
//     return m.cast<int>();
//     }

// int eigs_bdcsvd(DenseMatrix_PL &L){
//     Eigen::VectorXf eigenvalues;
//     auto start_eigs = std::chrono::high_resolution_clock::now();

//     Eigen::BDCSVD<DenseMatrix_PL> bdcsvd(L);
//     eigenvalues = bdcsvd.singularValues();
    
//     auto end_eigs = std::chrono::high_resolution_clock::now();
//     auto duration_eigs = std::chrono::duration_cast<std::chrono::milliseconds>(end_eigs - start_eigs);
//     // petls::print_vector_precise(eigenvalues);
//     return (int) duration_eigs.count();
// }


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

// int eigs_selfadjoint(DenseMatrix_PL &L){
//     Eigen::VectorXf eigenvalues;
//     auto start_eigs = std::chrono::high_resolution_clock::now();   

//     Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es(L,Eigen::EigenvaluesOnly);
//     eigenvalues = es.eigenvalues();

//     auto end_eigs = std::chrono::high_resolution_clock::now();
//     auto duration_eigs = std::chrono::duration_cast<std::chrono::milliseconds>(end_eigs - start_eigs);
//     // petls::print_vector_precise(eigenvalues);
//     return (int) duration_eigs.count();
// }

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


std::vector<double> get_filtrations(double delta){
    double start = 0.0;
    double end = 5.0;
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
    std::vector<int> dims = {0, 1, 2};
    int duration;
    Eigen::VectorXf eigenvalues;
    // output file
    std::ofstream out("./profiles/" + algorithm + "/r" + std::to_string(replicate) + ".csv");
    out << "dimension,a,b,simplices," << algorithm << std::endl;
    for (double filtration : filtrations){
        for (int dim : dims){
            

            
            std::cout << "Filtration: " << filtration << ", Dimension: " << dim;
            
            // open matrix file and read matrix
            std::stringstream stream;
            stream << "./torus_L/dim" << dim << "_a" << std::fixed << std::setprecision(2)
                << filtration << "_b" << std::fixed << std::setprecision(2) << filtration+delta << "_r"
                << replicate << ".mkt";
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
                out << dim << "," << filtration << "," << filtration + delta << "," << 0 << "," << 0 << std::endl;
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