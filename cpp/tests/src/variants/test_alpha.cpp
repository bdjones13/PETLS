#include "Complex.hpp"
#include "typedefs.hpp"
#include "petls.hpp"
#include "Alpha.hpp"
#include <iostream>
//#include "../helpers.hpp"

bool test_alpha_off(const char* filename){
    int max_dim = 3;
    petls::Alpha my_PAL(filename, max_dim);
    
    // my_PAL.spectra();
    std::cout << "Persistent Alpha Laplacian Spectra" << std::endl;
    petls::print_spectra(my_PAL.spectra());
    return true;
}

bool test_alpha_points(){
    int max_dim = 3;
    std::vector<std::tuple<double,double,double>> points;
    points.push_back(std::make_tuple(0.0,0.0,0.0));
    points.push_back(std::make_tuple(0.0,0.0,1.0));
    points.push_back(std::make_tuple(0.0,1.0,0.0));
    points.push_back(std::make_tuple(1.0,0.0,0.0));
    petls::Alpha my_PAL(points, max_dim);
    my_PAL.set_verbose(true);
    std::cout << "Persistent Alpha Laplacian Spectra From points" << std::endl;
    petls::print_spectra(my_PAL.spectra());
    return true;
}

int main(){
    bool passed = true;
    passed = passed && test_alpha_off("data/alpha/input");
    passed = passed && test_alpha_points();
}