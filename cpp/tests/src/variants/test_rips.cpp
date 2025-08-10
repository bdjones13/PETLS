// #include "pch.hpp"
// #include "typedefs.hpp"
#include "Rips.hpp"
#include <iostream>
#include "helpers.hpp"



petls::Rips rips_points(){
    int max_dim = 3;
    std::vector<double> p1 = {0, 0};
    std::vector<double> p2 = {0, 3};
    std::vector<double> p3 = {4, 0};
    std::vector<double> p4 = {4, 3};
    std::vector<std::vector<double>> points = {p1, p2, p3, p4};    
    petls::Rips my_rips(points, max_dim);
    return my_rips;
}

petls::Rips rips_points(float threshold){
    int max_dim = 3;
    std::vector<double> p1 = {0, 0};
    std::vector<double> p2 = {0, 3};
    std::vector<double> p3 = {4, 0};
    std::vector<double> p4 = {4, 3};
    std::vector<std::vector<double>> points = {p1, p2, p3, p4};    
    petls::Rips my_rips(points, max_dim, threshold);
    return my_rips;
}

petls::Rips rips_matrix(){
    int max_dim = 3;
    DenseMatrix_PL dists {
        {0, 0, 0, 0},
        {3, 0, 0, 0},
        {4, 5, 0, 0},
        {5, 4, 3, 0}
    };
    return petls::Rips(dists, max_dim);
}

petls::Rips rips_matrix(float threshold){
    int max_dim = 3;
    DenseMatrix_PL dists {
        {0, 0, 0, 0},
        {3, 0, 0, 0},
        {4, 5, 0, 0},
        {5, 4, 3, 0}
    };
    return petls::Rips(dists, max_dim, threshold);
}

petls::Rips rips_matrix_from_dists(){
    // std::cout << "from_distances" << std::endl;
    int max_dim = 3;
    DenseMatrix_PL dists {
        {0, 0, 0, 0},
        {3, 0, 0, 0},
        {4, 5, 0, 0},
        {5, 4, 3, 0}
    };
    return petls::Rips::from_distances(dists, max_dim);
}

bool test_rips(petls::Rips &my_rips){
    // Expects a rips complex that was constructed with points: (0,0), (0,3), (4,0), and (4,3) (4x3 rectangle)
    bool passed_all_tests = true;
    //dim 0
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,0,2,2}, 0, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 0, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,4,4,4}, 0, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,4,4,4}, 0, 5, 5);
    //dim 1
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 1, 0, 3);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {2,2}, 1, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {2,2,4,4}, 1, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {4,4,4,4,4,4}, 1, 5, 5);

    //dim 2
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {4,4,4,4}, 2, 5, 5);
    //dim 3
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {4}, 3, 5, 5);
    return passed_all_tests;
}


bool test_rips_threshold(petls::Rips &my_rips){
    // Expects a rips complex that was constructed with points: (0,0), (0,3), (4,0), and (4,3) (4x3 rectangle)
    // and threshold set to 4.5
    bool passed_all_tests = true;
    //dim 0
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,0,2,2}, 0, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 0, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 0, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 0, 5, 5);
    //dim 1
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 1, 0, 3);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {2,2}, 1, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 1, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {0,2,2,4}, 1, 5, 5);

    //dim 2
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 2, 5, 5);
    //dim 3
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 0, 3);   
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 3, 4);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 4, 5);
    passed_all_tests = passed_all_tests && test_sample(my_rips, {}, 3, 5, 5);
    return passed_all_tests;
}

int main(){
    
    
    bool passed_all_tests = true;
    petls::Rips my_rips_pts = rips_points();
    petls::Rips my_rips_file = petls::Rips("data/rips/rect.lower_distance_matrix",3);
    petls::Rips my_rips_matrix = rips_matrix();
    petls::Rips my_rips_matrix_from_dists = rips_matrix_from_dists();
    passed_all_tests = passed_all_tests && test_rips(my_rips_pts);
    passed_all_tests = passed_all_tests && test_rips(my_rips_file); 
    passed_all_tests = passed_all_tests && test_rips(my_rips_matrix);
    passed_all_tests = passed_all_tests && test_rips(my_rips_matrix_from_dists);
    
    float threshold = 4.5;
    petls::Rips my_rips_pts_threshold = rips_points(threshold);
    petls::Rips my_rips_file_threshold = petls::Rips("data/rips/rect.lower_distance_matrix",3, threshold);
    petls::Rips my_rips_matrix_threshold = rips_matrix(threshold);
    passed_all_tests = passed_all_tests && test_rips_threshold(my_rips_pts_threshold);
    passed_all_tests = passed_all_tests && test_rips_threshold(my_rips_file_threshold); 
    passed_all_tests = passed_all_tests && test_rips_threshold(my_rips_matrix_threshold);
    
    if (passed_all_tests){
        std::cout << "rips passed all tests" << std::endl;
        return 0;
    } else{
        return -1;
    }
}
