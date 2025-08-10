// #include "pch.hpp"
// #include "Complex.hpp"
// #include "typedefs.hpp"
// #include "petls.hpp"
#include "../helpers.hpp"
#include <vector>
// #include <utility>
#include <cmath> //sqrt

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <iostream>


void print_all_eigenpairs(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL>> eigenpairs){
    
    for (int i = 0; i < (int)eigenpairs.size(); i++){
        std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL> cur = eigenpairs[i];
        std::cout << "\n dim = " << std::get<0>(cur) << ", a = " << std::get<1>(cur) << ", b = " << std::get<2>(cur);
        std::cout << ", eigs = [";
        std::vector<spectra_type> eigs = std::get<3>(cur);
        for(int j = 0; j < (int) eigs.size(); j++){
            std::cout << eigs[j] << ", ";
        }
        std::cout << "], eigenvectors = " << std::get<4>(cur) << std::endl;
    }
}

template<class up_alg>
bool test_filtered_triangle_complex_vec(petls::Complex<int, up_alg> &complex){
    float x = std::sqrt(2)/2;
    DenseMatrix_PL refv {{x,x },{x,-x}};
    // bool passed_all_tests = test_sample(complex, {0,2}, refv, 0, 1.0, 3.0);
    // std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL>> eigenpairs = complex.eigenpairs();
    // print_all_eigenpairs(eigenpairs);

    return test_sample(complex, {0,2}, refv, 0, 1.0, 3.0);
}



int test_filtered_triangle_eigenvectors(){
    petls::Complex<int> mycomplex = filtered_triangle_complex_alg<petls::schur<int>>();
    bool passed_all_tests = true;
    passed_all_tests = passed_all_tests & test_filtered_triangle_complex_vec(mycomplex);
    if (passed_all_tests)
        return 0;
    else
        return -1;
}


int main(){
    return -1;// fail
    // return test_filtered_triangle_eigenvectors();
}