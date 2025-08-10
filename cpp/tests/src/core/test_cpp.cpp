
#include <tuple>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <string>
#include "filtered_triangle.cpp"
// #include "test_eigenvectors.cpp"
#include <chrono>

typedef float spectra_type;

SparseMatrixInt read_eigen_sparse(std::string filename){
        SparseMatrixDouble m;
        bool s = Eigen::loadMarket(m, filename);
        if (!s){
            std::cout << "failed to load Eigen::loadMarket, filename: " << filename << std::endl;
        }
        return m.cast<int>();
    }

std::vector<filtration_type> read_filtration(std::string filename){
        std::ifstream instream(filename);
        std::string str_contents;
        std::string token;
        std::string delimeter = ",";
        size_t pos = 0;
        if (instream.is_open()){
            instream >> str_contents;
        }

        std::vector<filtration_type> filtrations;
        while ((pos = str_contents.find(delimeter)) != std::string::npos){
            token = str_contents.substr(0,pos);
            filtrations.push_back(std::stod(token));
            // std::cout << token << std::endl;
            str_contents.erase(0, pos + delimeter.length());
        }
        // std::cout << str_contents << std::endl;
        return filtrations;
    }
// void test_MWW(){
//     // std::cout << "the first thing" << std::endl;
//     SparseMatrixInt d1 = read_eigen_sparse("../data/MWW_example/matrix_1");
//     SparseMatrixInt d2 = read_eigen_sparse("../data/MWW_example/matrix_2");
//     std::vector<filtration_type> c0 = read_filtration("../data/MWW_example/_filtrations_0.txt");
//     std::vector<filtration_type> c1 = read_filtration("../data/MWW_example/_filtrations_1.txt");
//     std::vector<filtration_type> c2 = read_filtration("../data/MWW_example/_filtrations_2.txt");
//     std::vector<SparseMatrixInt> boundaries;
//     boundaries.push_back(d1);
//     boundaries.push_back(d2);
//     std::vector<std::vector<filtration_type>> filtrations;
//     filtrations.push_back(c0);
//     filtrations.push_back(c1);
//     filtrations.push_back(c2);
//     std::cout << "assemble" << std::endl;
//     petls::Complex complex(boundaries,filtrations);
//     // complex.print_boundaries();

//     // petls::print_vector_precise(c1);
//     // petls::print_vector_precise(c2);
//     auto start_spectra = std::chrono::high_resolution_clock::now();

//     // std::vector<std::tuple<int, filtration_type, filtration_type, Eigen::VectorXd>> all_spectra = complex.spectra();
//     // std::cout << "got all spectra. printing" << std::endl;
//     // for (int i = 0; i < all_spectra.size(); i++){
//     //     std::tuple<int, filtration_type, filtration_type, Eigen::VectorXd> requested_spectra = all_spectra[i];
//     //     std::cout << "dim=" << std::get<0>(requested_spectra) << ", a=" << std::get<1>(requested_spectra) << ", b=" << std::get<2>(requested_spectra) << ", spectra=";
//     //     // std::cout << "got " << i << "th specta vector" << std::endl;
//     //     petls::print_vector_precise(std::get<3>(requested_spectra));
//     // }
//     // std::cout << "printed all spectra" << std::endl;


//     // petls::print_vector_precise(complex.spectra(0, 0.9, 1.6));
    
//     print_std_vec(complex.spectra(1, 1.0, 1.5));
    
//     // petls::print_vector_precise(complex.spectra(2, 1.0, 1.6));
//     std::cout << "got spectra" << std::endl;
//     auto end_spectra = std::chrono::high_resolution_clock::now();
//     auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
//     std::cout << "duration spectra MWW_example: " << duration_spectra.count() << std::endl;

// }

// void test_1a99_cutoff15(){
    
//     std::string prefix = "../data/1a99_cutoff15/";
//     SparseMatrixInt d1 = read_eigen_sparse(prefix + std::string("matrix_0"));
//     SparseMatrixInt d2 = read_eigen_sparse(prefix + std::string("matrix_1"));

//     std::vector<filtration_type> c0 = read_filtration(prefix + "_filtrations_0.txt");
//     std::vector<filtration_type> c1 = read_filtration(prefix + "_filtrations_1.txt");
//     std::vector<filtration_type> c2 = read_filtration(prefix + "_filtrations_2.txt");
//     std::vector<SparseMatrixInt> boundaries;
//     boundaries.push_back(d1);
//     boundaries.push_back(d2);
//     std::vector<std::vector<filtration_type>> filtrations;
//     filtrations.push_back(c0);
//     filtrations.push_back(c1);
//     filtrations.push_back(c2);

//     petls::Complex complex(boundaries,filtrations);
//     // complex.print_boundaries();
//     // std::cout << "d1= " << d1.rows() << " x " << d1.cols() << ", d2 = " << d2.rows() << " x " << d2.cols() << std::endl;
//     // std::cout << "c0 =" << c0.size() << ", c1 =" << c1.size() << ", c2 =" << c2.size();
//     // std::cout << " c2=" << c2.rows() << "x" << c2.cols() << std::endl; 

//     auto start_spectra = std::chrono::high_resolution_clock::now();

//     // std::vector<std::tuple<int, filtration_type, filtration_type, Eigen::VectorXd>> all_spectra = complex.spectra();
//     // std::cout << "got all spectra. printing" << std::endl;
//     // for (int i = 0; i < all_spectra.size(); i++){
//     //     std::tuple<int, filtration_type, filtration_type, Eigen::VectorXd> requested_spectra = all_spectra[i];
//     //     std::cout << "dim=" << std::get<0>(requested_spectra) << ", a=" << std::get<1>(requested_spectra) << ", b=" << std::get<2>(requested_spectra) << ", spectra=";
//     //     // std::cout << "got " << i << "th specta vector" << std::endl;
//     //     petls::print_vector_precise(std::get<3>(requested_spectra));
//     // }
//     // std::cout << "printed all spectra" << std::endl;

//     complex.spectra();
//     // petls::print_vector_precise(complex.spectra(0, 0.9, 1.6));
    
//     // petls::print_vector_precise(complex.spectra(1, 4.0, 4.6));
    
//     // petls::print_vector_precise(complex.spectra(2, 1.0, 1.6));
//     std::cout << "got all spectra " + prefix << std::endl;
//     auto end_spectra = std::chrono::high_resolution_clock::now();
//     auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
//     std::cout << "duration spectra " << prefix << ": " << duration_spectra.count() << std::endl;

// }

// void test_1a99_cutoff8(){
    
//     std::string prefix = "../data/1a99_cutoff8/";
//     SparseMatrixInt d1 = read_eigen_sparse(prefix + std::string("matrix_0"));
//     SparseMatrixInt d2 = read_eigen_sparse(prefix + std::string("matrix_1"));

//     std::vector<filtration_type> c0 = read_filtration(prefix + "_filtrations_0.txt");
//     std::vector<filtration_type> c1 = read_filtration(prefix + "_filtrations_1.txt");
//     std::vector<filtration_type> c2 = read_filtration(prefix + "_filtrations_2.txt");
//     std::vector<SparseMatrixInt> boundaries;
//     boundaries.push_back(d1);
//     boundaries.push_back(d2);
//     std::vector<std::vector<filtration_type>> filtrations;
//     filtrations.push_back(c0);
//     filtrations.push_back(c1);
//     filtrations.push_back(c2);

//     petls::Complex complex(boundaries,filtrations);

//     // std::cout << "d1= " << d1.rows() << " x " << d1.cols() << ", d2 = " << d2.rows() << " x " << d2.cols() << std::endl;
//     // std::cout << "c0 =" << c0.size() << ", c1 =" << c1.size() << ", c2 =" << c2.size();


//     auto start_spectra = std::chrono::high_resolution_clock::now();

 
//     std::vector<std::vector<std::vector<float>>> my_spectra;
//     for (int i = 0; i <= 3; i++){
//         my_spectra.push_back(std::vector<std::vector<float>>());
//     }
//     std::vector<std::tuple<int, double, double, std::vector<float>>> computed_results = complex.spectra();
//     for (int i = 0; i < (int) computed_results.size(); i++){
//             int dim = std::get<0>(computed_results[i]);
//             // double a = std::get<1>(computed_results[i]); //not used
//             // double b = std::get<2>(computed_results[i]); //not used
//             std::vector<float> current_spectra = std::get<3>(computed_results[i]);
//             my_spectra[dim].push_back(current_spectra);
//         }

//     std::cout << "got all spectra " + prefix << std::endl;
//     auto end_spectra = std::chrono::high_resolution_clock::now();
//     auto duration_spectra = std::chrono::duration_cast<std::chrono::milliseconds>(end_spectra - start_spectra);
//     std::cout << "duration spectra " << prefix << ": " << duration_spectra.count() << std::endl;

// }


int main(int argc, char **argv) {
//     petls::Complex defaulted;
//     return 0;
    // std::cout << "size of float " << sizeof(float) << std::endl;

    if (argc > 1){
        std::cout << argv[1] << std::endl;
        if (std::string(argv[1]) == "8"){
            // test_1a99_cutoff8();
        }
        else if (std::string(argv[1]) == "MWW"){
            std::cout << "test_MWW()" << std::endl;
            // test_MWW();
        }
        // // else if (std::string(argv[1]) == "alpha"){
        // //     std::cout << "testing alpha" << std::endl;
        // //     test_alpha_off(argv[2]);
        // // }
        // // else if (std::string(argv[1]) == "alpha_p"){
        // //     std::cout << "testing alpha points" << std::endl;
        // //     test_alpha_points();
        // // }
        // else if (std::string(argv[1]) == "save"){
        //     std::cout << "testing store_L" << std::endl;
        //     test_save();
        // }
        // else if (std::string(argv[1]) == "eigs"){
        //     std::cout << "testing eigenvalue algorithms" << std::endl;
        //     return test_filtered_triangle_eigs_algs();
        // }

    }
    return 0;
}