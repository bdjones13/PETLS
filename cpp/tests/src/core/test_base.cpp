#include "Complex.hpp"
#include "../helpers.hpp"
#include <fstream>
using PTL = petls::Complex;

int test_specific(){
    PTL pl = filtered_triangle_complex();
    float tol = 1e-3;
    bool passed = true;
    passed = passed && compare_spectra(pl.spectra(0,5,6),{0,3,3},tol);
    if (passed){
        return 0;
    } else{
        return -1;
    }
}

bool test_L(){
    PTL pl = filtered_triangle_complex();
    Eigen::MatrixXf L = pl.get_L(1,5,5);
    // petls::print_full_matrix_precise(L);
    std::vector<float> diag = {3,3,3};
    Eigen::Vector3f ref_diag(diag.data());
    Eigen::MatrixXf ref = ref_diag.asDiagonal();
    // petls::print_full_matrix_precise(ref);
    return L.isApprox(ref); 
}

bool test_up(){
    PTL pl = filtered_triangle_complex();
    Eigen::MatrixXf up = pl.get_up(0,1,3);
    Eigen::Matrix2f ref;
    ref << 1, -1, -1, 1;
    return up.isApprox(ref);
}

bool test_down(){
    PTL pl = filtered_triangle_complex();
    Eigen::MatrixXf down = pl.get_down(1,5);
    Eigen::Matrix3f ref;
    ref << 2, -1, 1, 
            -1, 2, 1,
            1, 1, 2;
    return down.isApprox(ref);
}

bool test_sum_up_down(){
    PTL pl = filtered_triangle_complex();
    Eigen::MatrixXf down = pl.get_down(1,5);
    Eigen::MatrixXf up = pl.get_up(1,5,5);
    Eigen::MatrixXf L = pl.get_L(1,5,5);
    return L.isApprox(up+down.cast<float>());
}

bool test_profile(){
    petls::Complex complex = filtered_triangle_complex();
    complex.spectra(0,1,2);
    complex.spectra(0,2,3);
    complex.spectra(1,2,3);
    std::string filename = "test_profile.csv";
    std::ifstream file;
    file.open(filename);

    if(file.is_open()){
        std::cout << "test_profile.csv already exists. Deleting before proceeding." << std::endl;
        file.close();
        std::remove(filename.c_str());
    }

    complex.time_to_csv(filename);
    file.open(filename);
    bool status = file.good();
    file.close();
    return status;
}

int main(){
    bool passed = true;
    passed = passed && test_L();
    passed = passed && test_up();
    passed = passed && test_down();
    passed = passed && test_sum_up_down();
    passed = passed && test_profile();
    int status = passed ? 0 : -1;
    std::cout << "status: " << status << std::endl;
    return passed ? 0 : -1;
    // return test_specific();
}