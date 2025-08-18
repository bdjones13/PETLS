#include "filtered_triangle.cpp"
#include "petls.hpp"
#include "up_algorithms.hpp"
#include "Eigen/Eigenvalues"

// Copy of the schur_algorithm function for testing purposes
inline void Test_Schur(petls::FilteredBoundaryMatrix<int>* fbm, filtration_type a, filtration_type b, DenseMatrix_PL &L_up){
        int a_row_index = fbm->index_of_filtration(false,a);
        int b_row_index = fbm->index_of_filtration(false,b);
        int b_col_index = fbm->index_of_filtration(true,b);

        if (a_row_index == b_row_index){
            Eigen::SparseMatrix<int> B_pers(b_row_index, b_col_index);
            fbm->submatrix_at_filtration(b, B_pers);
            L_up = (B_pers*B_pers.transpose()).template cast<coefficient_type>();
            return;
        } else if (a_row_index == -1){ // no rows, return 0x0 empty matrix
            L_up.setZero(0,0);
            return;

        } else if (b_col_index == -1){ // no columns, return nxn 0-matrix for n = number of cells in C_{up_dim}^a
            L_up.setZero(a_row_index+1,a_row_index+1);
            return;

        }
        // NOTE: doing the self adjoint lower triangular trick here is slower
        // possibly because of the conversions between dense and sparse

        Eigen::SparseMatrix<int> B_pers_int(b_row_index, b_col_index);
        fbm->submatrix_at_filtration(b, B_pers_int);
        
    
        Eigen::SparseMatrix<int> L_up_b_int(B_pers_int.rows(), B_pers_int.rows());
        // L_up_b_int = B_pers_int*B_pers_int.transpose();
        L_up_b_int.template selfadjointView<Eigen::Lower>().rankUpdate(B_pers_int);
        SparseMatrix_PL L_up_b(L_up_b_int.rows(), L_up_b_int.cols());
        L_up_b = L_up_b_int.template cast<coefficient_type>();

        int a_rows = a_row_index + 1;
        int b_rows = b_row_index + 1;

        
        
        // split L_up_b into 4 along those borders
        // A B
        // C D
        // and compute A - B*D^{-1}*C
        // by the structure of L, we know B = C^T
        // and we want A - B * (things)
        // Also compute D^{-1}*C by solving the sparse lower-triangular linear system

        SparseMatrix_PL A = L_up_b.topLeftCorner(a_rows, a_rows);
        SparseMatrix_PL C = L_up_b.bottomLeftCorner(b_rows-a_rows, a_rows);
        SparseMatrix_PL B = C.transpose();
        SparseMatrix_PL D = L_up_b.bottomRightCorner(b_rows-a_rows, b_rows-a_rows);
        
        L_up = A - B * DenseMatrix_PL(D).ldlt().solve(DenseMatrix_PL(C));
        
        // The computed matrix is only correct in the lower triangular portion.
        // Eigen::SelfAdjointEigensolver is okay with this, but to avoid confusion and issues with other eigensolvers,
        // We symmetrize the matrix:
        L_up = DenseMatrix_PL(L_up.selfadjointView<Eigen::Lower>());
        
        return;
    }



int test_up_algs(){
    // return -1;
    //should be spectra(1,4,5) = [1,3]
    petls::Complex mycomplex = filtered_triangle_complex_alg(); // default is schur
    petls::Complex mycomplex2 = filtered_triangle_complex_alg();
    petls::Complex mycomplex3 = filtered_triangle_complex_alg();
    petls::Complex mycomplex4 = filtered_triangle_complex_alg();


    mycomplex2.set_eigs_algorithm_func(petls::schur_algorithm<int>);
    mycomplex3.set_eigs_algorithm_func("schur");
    mycomplex4.set_eigs_algorithm_func(Test_Schur);

    bool passed_all_tests = true;
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex2);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex3);
    if (passed_all_tests)
        return 0;
    else
        return -1;
}

int main(){
    return test_up_algs();
}