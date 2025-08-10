#include "dFlag.hpp"
#include <iostream>

#include "../core/filtered_triangle.cpp"
#include <chrono>



int test_flag_transpose_flip(){
    const char* filename = "data/flag/triangle.flag";
    int max_dim = 3;
    petls::dFlag my_dflag(filename, max_dim);

    std::cout << "Persistent Directed Flag Laplacian Spectra" << std::endl;
    my_dflag.set_verbose(true);
    my_dflag.set_flipped(true);
    bool passed = test_filtered_triangle_complex(my_dflag);
    if (passed){
        return 0;
    } else{
        return -1;
    }
}

int main(){
    return test_flag_transpose_flip();
}