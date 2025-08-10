#include "../helpers.hpp"
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <iostream>


bool test_filtered_triangle_complex(petls::Complex &complex){

    bool passed_all_tests = true;

    // dim 0
    passed_all_tests = passed_all_tests && test_sample(complex, {0}, 0, 0.0, 0.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0}, 0, 0.0, 1.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0, 0}, 0, 1.0, 2.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0, 0, 2}, 0, 2.0, 3.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0, 1, 3}, 0, 3.0, 4.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0, 3, 3}, 0, 4.0, 5.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {0, 3, 3}, 0, 5.0, 6.0);

    // dim 1
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 1, 0.0, 0.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 1, 0.0, 1.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 1, 1.0, 2.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 1, 2.0, 3.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {2}, 1, 3.0, 4.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {1, 3}, 1, 4.0, 5.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {3, 3, 3}, 1, 5.0, 6.0);

    // dim 2
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 0.0, 0.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 0.0, 1.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 1.0, 2.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 2.0, 3.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 3.0, 4.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {}, 2, 4.0, 5.0);
    passed_all_tests = passed_all_tests && test_sample(complex, {3}, 2, 5.0, 6.0);

    return passed_all_tests;
}


bool test_filtered_triangle_all(petls::Complex complex){
    std::cout << "call all spectra" << std::endl;
    std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> all_spectra = complex.spectra();
    std::cout << "got all spectra. printing" << std::endl;
    for (int i = 0; i < (int) all_spectra.size(); i++){
        std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>> requested_spectra = all_spectra[i];
        std::cout << "dim=" << std::get<0>(requested_spectra) << ", a=" << std::get<1>(requested_spectra) << ", b=" << std::get<2>(requested_spectra) << ", spectra=";
        // std::cout << "got " << i << "th specta vector" << std::endl;
        print_std_vec(std::get<3>(requested_spectra));
    }
    std::cout << "printed all spectra" << std::endl;
    return true; //TODO turn into a test of success
}


// bool test_filtered_triangle_eigs_algs(){
//     //should be spectra(1,4,5) = [1,3]
//     petls::Complex<petls::selfadjoint> mycomplex = filtered_triangle_complex_alg<petls::selfadjoint,petls::schur<>>();
//     petls::Complex<petls::bdcsvd> mycomplex2 = filtered_triangle_complex_alg<petls::bdcsvd,petls::schur<>>();
//     petls::Complex<petls::eigensolver> mycomplex3 = filtered_triangle_complex_alg<petls::eigensolver, petls::schur<>>();
    
//     bool passed_all_tests = true;
//     passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex);
//     passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex2);
//     passed_all_tests = passed_all_tests && test_filtered_triangle_complex(mycomplex3);
//     passed_all_tests = false;//force failure
//     if (passed_all_tests)
//         return 0;
//     else
//         return -1;
// }