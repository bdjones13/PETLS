#include "helpers.hpp"

// #include "Complex.hpp"
// #include <unsupported/Eigen/SparseExtra>
// #include <Eigen/Sparse>
// #include <Eigen/Dense>

#include <vector>

// TEST_CASE("reduction larger"){
//     int max_dim = 3;
    

//     petls::dFlag my_dflag("data/flag/1a99_8.flag", max_dim);
//     // my_dflag.print_boundaries();
//     // my_PAL.spectra();
//     std::cout << "Persistent Directed Flag Laplacian Spectra" << std::endl;
//     my_dflag.set_verbose(true);
//     SparseMatrixFloat dummy;
//     auto start_spectra = std::chrono::high_resolution_clock::now();
//     auto spectra = my_dflag.nonzero_spectra(2,7.996,7.996,dummy,true);
//     // my_dflag.store_spectra(spectra,"myprefix");
//     // my_dflag.store_spectra_summary(spectra,"myprefix");
//     std::cout << "Finished Persistent Directed Flag Laplacian Spectra" << std::endl;
//     auto end_spectra = std::chrono::high_resolution_clock::now();
//     auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
//     std::cout << "duration directed flag  spectra: " << duration_spectra.count() << std::endl;
//     REQUIRE(0 == 1);
// }

// TEST_CASE("non-reducation large"){
//     int max_dim = 3;
    

//     petls::dFlag my_dflag("data/flag/1a99_8.flag", max_dim);
//     // my_dflag.print_boundaries();
//     // my_PAL.spectra();
//     std::cout << "Persistent Directed Flag Laplacian Spectra" << std::endl;
//     my_dflag.set_verbose(true);
//     SparseMatrixFloat dummy;
//     auto start_spectra = std::chrono::high_resolution_clock::now();
//     auto spectra = my_dflag.spectra(2,7.996,7.996);
//     // my_dflag.store_spectra(spectra,"myprefix");
//     // my_dflag.store_spectra_summary(spectra,"myprefix");
//     std::cout << "Finished Persistent Directed Flag Laplacian Spectra" << std::endl;
//     auto end_spectra = std::chrono::high_resolution_clock::now();
//     auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
//     std::cout << "duration directed flag  spectra: " << duration_spectra.count() << std::endl;
//     REQUIRE(0 == 1);
// }

int test_reduction(){
    //Boundary matrices corresponding to directed flag laplacian of
    // dim 0:
    // 0 1 2
    // dim 1:
    // 0 1 3
    // 1 2 4
    // 0 2 5


    /* Visually:
    
           0
          / \
       3 /   \ 5  
        v     v
        1 ---> 2
           4
    */
    petls::Complex pl = filtered_triangle_complex();
   
    // basis of PH 2,3 would be [1,1,0] and [0,0,1] 
    // because the connected components alive at 2 and persist until 3 can be represented by v_0,v_1 and v_2
    SparseMatrixFloat ph2_3(3,2);
    ph2_3.coeffRef(0,0) = 1;
    ph2_3.coeffRef(1,0) = 1;
    ph2_3.coeffRef(2,1) = 1;
    bool use_dummy_harmonic_basis = false;
    std::vector<spectra_type> result2_3 = pl.nonzero_spectra(0,2,3, ph2_3, use_dummy_harmonic_basis); // should be 0,0,2 -> 2

    float tol = 1e-1;
    bool passed = true;
    passed = passed && compare_spectra(result2_3,{2},tol);

    // basis of PH  3,4 would be [1, 1, 1] 
    // because the connected component alive at 3 that persists until 4 can be represented by v_0, v_1, v_2
    SparseMatrixFloat ph3_4(3,1);
    ph3_4.coeffRef(0,0) = 1;
    ph3_4.coeffRef(1,0) = 1;
    ph3_4.coeffRef(2,0) = 1;
    std::vector<spectra_type> result3_4 = pl.nonzero_spectra(0,3,4,ph3_4, use_dummy_harmonic_basis); // should be 0,1,3 -> 1,3

    
    passed = passed && compare_spectra(result3_4,{1,3},tol);
    if (passed){
        return 0;
    } else{
        return -1;
    }
}

int main(){
    return test_reduction();
}
