//#include "../helpers.hpp"

#include "Complex.hpp"
#include "typedefs.hpp"
#include "petls.hpp"
#include <vector>
#include <string>

void test_save(){
//Boundary matrices corresponding to directed flag laplacian of
    // dim 0:
    // 0 0 0 1.0
    // dim 1:
    // 0 1 0
    // 1 2 0
    // 0 2 0
    // 2 3 1.0 
    SparseMatrixInt d1(4,4);
    SparseMatrixInt d2(4,1);
    d1.coeffRef(0,0) = -1;
    d1.coeffRef(1,0) = 1;
    d1.coeffRef(1,1) = -1;
    d1.coeffRef(2,1) = 1;
    d1.coeffRef(0,2) = -1;
    d1.coeffRef(2,2) = 1;
    d1.coeffRef(2,3) = -1;
    d1.coeffRef(3,3) = 1;

    d2.coeffRef(0,0) = 1;
    d2.coeffRef(1,0) = 1;
    d2.coeffRef(2,0) = -1;

    std::vector<filtration_type> c0_filtrations = {0.0, 0.0, 0.0, 1.0};
    std::vector<filtration_type> c1_filtrations = {0.0, 0.0, 0.0, 1.0};
    std::vector<filtration_type> c2_filtrations = {1.0};

    std::vector<SparseMatrixInt> boundaries;
    boundaries.push_back(d1);
    boundaries.push_back(d2);
    std::vector<std::vector<filtration_type>> filtrations;
    filtrations.push_back(c0_filtrations);
    filtrations.push_back(c1_filtrations);
    filtrations.push_back(c2_filtrations);
    
    petls::Complex<int, petls::schur<>> pl(boundaries,filtrations);
    pl.store_L(0,0,1,"saved_matrix.mtx");
}