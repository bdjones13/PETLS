#include "Complex.hpp"
#include "Rips.hpp"
#include <vector>
#include <iostream>

int main() {
    std::vector<std::vector<double>> points = {
        {0.0, 0.0},
        {3.0, 0.0},
        {0.0, 4.0},
        {3.0, 4.0}
    };

    petls::Rips rips(points, 3);
    std::cout << "Constructed Rips complex" << std::endl;
    
    std::vector<float> eigs = rips.spectra(1, 3, 4); // result is [2.0, 2.0]
    
    // print
    std::cout << "Eigenvalues of L_{1}^{3,4} are:";
        std::cout << "[";
        if (eigs.size())
        for (int j = 0; j < (int) eigs.size() - 1; j++){
            std::cout << eigs[j] << ", ";
        }
        if (eigs.size() > 0) // if no eigs, don't try to print 'last'
            std::cout << eigs[eigs.size()-1];
        std::cout << "]" << std::endl;

    // all pairs of spectra (a)
    std::vector<std::tuple<int, double, double, std::vector<float>>> s = rips.spectra();
    for (int i = 0; i < (int) s.size(); i++){
        std::tuple<int, double, double, std::vector<float>> spectrum = s[i];
        int dim = std::get<0>(spectrum);
        double a = std::get<1>(spectrum);
        double b = std::get<2>(spectrum);
        std::vector<float> eigs = std::get<3>(spectrum);
        
        std::cout << "Eigenvalues of L_{dim}^{a,b} for dim = " << dim << ", a = " << a;
        std::cout << ", b = " << b << " are:";
        std::cout << "[";
        if (eigs.size())
        for (int j = 0; j < (int) eigs.size() - 1; j++){
            std::cout << eigs[j] << ", ";
        }
        if (eigs.size() > 0) // if no eigs, don't try to print 'last'
            std::cout << eigs[eigs.size()-1];
        std::cout << "]" << std::endl;
    }
    std::cout << "Finished spectra" << std::endl;
    rips.store_spectra(s, "myprefix"); // check myprefix_spectra_d.txt for d=0, 1, 2, 3
}