#include "filtered_triangle.cpp"
#include "petls.hpp"
#include "eigs_algorithms.hpp"
#include "Eigen/Eigenvalues"


inline spectra_vec Test_SelfAdjointEigen(DenseMatrix_PL& L) {
    Eigen::SelfAdjointEigenSolver<DenseMatrix_PL> es = Eigen::SelfAdjointEigenSolver<DenseMatrix_PL>(L, Eigen::EigenvaluesOnly);
    spectra_vec eigs = es.eigenvalues();
    petls::round_zeros(eigs, 1e-3);
    return eigs;
}

int test_filtered_triangle_eigs_algs(){
    // return -1;
    //should be spectra(1,4,5) = [1,3]
    petls::Complex mycomplex = filtered_triangle_complex_alg();
    petls::Complex mycomplex2 = filtered_triangle_complex_alg();
    petls::Complex mycomplex3 = filtered_triangle_complex_alg();
    petls::Complex mycomplex4 = filtered_triangle_complex_alg();
    petls::Complex mycomplex5 = filtered_triangle_complex_alg();
    petls::Complex mycomplex6 = filtered_triangle_complex_alg();
    petls::Complex mycomplex7 = filtered_triangle_complex_alg();
    
    mycomplex2.set_eigs_algorithm_func(petls::EigensolverEigen);
    mycomplex3.set_eigs_algorithm_func(petls::BDCSVDEigen);
    mycomplex4.set_eigs_algorithm_func("selfadjoint");
    mycomplex5.set_eigs_algorithm_func("eigensolver");
    mycomplex6.set_eigs_algorithm_func("bdcsvd");
    mycomplex7.set_eigs_algorithm_func(Test_SelfAdjointEigen);


    bool passed_all_tests = true;
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex2);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex3);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex4);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex5);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex6);
    passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex7);
    if (passed_all_tests)
        return 0;
    else
        return -1;
}

int main(){
    return test_filtered_triangle_eigs_algs();
}