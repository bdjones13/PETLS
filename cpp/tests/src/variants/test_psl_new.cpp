// #include "pch.hpp"
// #include "Complex.hpp"
// #include "typedefs.hpp"
// #include "petls.hpp"
#include "newPSL.hpp"
#include "sheaf_simplex_tree.hpp"
#include <iostream>
// #include "../helpers.hpp"


bool test_sample(petls::NewPersistentSheafLaplacian &complex, std::vector<spectra_type> reference_spectra, int dim, double a, double b){
    spectra_vec ref_spectra_eigen = Eigen::Map<spectra_vec, Eigen::Unaligned>(reference_spectra.data(),reference_spectra.size());
    std::vector<spectra_type> spectra_std_vec = complex.spectra(dim,a,b);
    spectra_vec spectra = Eigen::Map<spectra_vec, Eigen::Unaligned>(spectra_std_vec.data(),spectra_std_vec.size());
    float tol = 1e-1;

    if (spectra.size() != ref_spectra_eigen.size()){
        std::cout << "Fail compare spectra, different size vectors. a = " << a << ", b = " << b << ", dim = " << dim <<  std::endl;
        std::cout << "spectra.size() = " << spectra.size() << ", ref_spectra_eigen.size() = " << ref_spectra_eigen.size() << std::endl;
        std::cout << "Expected:";
        petls::print_vector_precise(ref_spectra_eigen);
        std::cout << "Actual:";
        petls::print_vector_precise(spectra);
        return false;
    }

    if(ref_spectra_eigen.isApprox(spectra, tol)){
        return true;
    } else{
        std::string message = "Fail spectra a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", dim=" + std::to_string(dim);
        std::cout << message << std::endl;
        std::cout << "Expected:";
        petls::print_vector_precise(ref_spectra_eigen);
        std::cout << "Actual:";
        petls::print_vector_precise(spectra);
        return false;
    }
}

using Simplex_tree = petls::Simplex_tree; 
float my_restriction(Simplex_tree::Simplex_handle simplex, Simplex_tree::Simplex_handle coface, petls::sheaf_simplex_tree& sst){
    // User will write their own restriction function

    /********************Read this comment block carefully***************************************
    *    If you do not satisfy the conditions below,                                            *
    *                                                                                           *
    *    YOU WILL NOT HAVE A CELLULAR SHEAF                                                     *
    *                                                                                           *
    *    and any conclusions based on having a cellular sheaf will not be valid.                *
    *                                                                                           *
    *    You must seperately verify that your restriction satisfies the conditions.             *
    *********************************************************************************************
    *                                                                                           *
    * The restriction map MUST satisfy the conditions:                                          *
    * (1) Composition: if a is a face of b, denoted a <= b, and b <= c, then                    *
    *                     restriction(b, c) * restriction(a,b) = restriction (a, c)             *
    * (2) Identity: restriction(a, a) = 1                                                       *
    *                                                                                           *
    * Note that the sheaf simplex tree only ever applies restriction from a simplex             *
    *     to a codimension-1 coface (vertex to edge or edge to 2-simplex, etc),                 *
    *     but that does not prevent you from accidentally defining an invalid restriction map   *
    *                                                                                           *
    *********************************************************************************************/

    // This is an example using the restriction from example 2.5 of
    // Xiaoqi Wei and Guo-Wei Wei, "Persistent Sheaf Laplacians," Foundations of Data Science, 2024. 
    // http://dx.doi.org/10.3934/fods.2024033
    // 
    // This sheaf is given by: 
    //      stalks are 1-dimensional R
    //      restriction of vertex i to edge ij is multiplication by charge_j / distance_ij
    //      restriction edge ij to 2-simplex  ijk is multiplication by charge_k / (distance_ik * distance_jk)

    // restriction of vertex to an edge
    if (sst.st.dimension(simplex) == 0){ 
        
        // get the other vertex of the edge
        std::pair<Simplex_tree::Simplex_handle, Simplex_tree::Simplex_handle> endpoints = sst.st.endpoints(coface);
        Simplex_tree::Simplex_key sibling;
        Simplex_tree::Simplex_key key = sst.st.key(simplex); 
        if (sst.st.key(endpoints.first) == key)
            sibling = sst.st.key(endpoints.second);
        else 
            sibling = sst.st.key(endpoints.first);

        // distance
        double distance = sqrt(pow(sst.extra_data[key][0]-sst.extra_data[sibling][0],2)
                                +pow(sst.extra_data[key][1]-sst.extra_data[sibling][1],2)
                                +pow(sst.extra_data[key][2]-sst.extra_data[sibling][2],2));        
        // distance = coface.filtration() // this would work fine instead, 
        //                                   but the above illustrates the use of the extra_data
        
        // charge / distance
        return sst.extra_data[sibling][3]/distance;        
    } else if (sst.st.dimension(simplex) == 1){
        // restriction of edge to a 2-simplex

        float coeff = 1.0;

        // For the simplex ijk, loop over the pairs (edge, missing_vertex): (jk, i), (ik, j), (ij, k)
        Simplex_tree::Boundary_opposite_vertex_simplex_range bovsr = sst.st.boundary_opposite_vertex_simplex_range(coface);
        for (std::pair<Simplex_tree::Simplex_handle, Simplex_tree::Vertex_handle> p : bovsr){
            Simplex_tree::Simplex_handle sibling_sh = std::get<0>(p);
            int sibling_key = sst.st.key(sibling_sh);
            Simplex_tree::Vertex_handle opposite_vh = std::get<1>(p);
            

            // if we want restriction of ij to ijk,
            // when the sibling of ij in ijk is ij itself, multiply by the charge at vertex k 
            if(sibling_key == (int) sst.st.key(simplex)){
                coeff *= sst.extra_data[opposite_vh][3];//charge
            } else{
                // when the sibling is jk or ik, divide by distance_jk or distance_ik, respectively
                coeff /= sst.st.filtration(sibling_sh);
            }
        }
        return coeff;
    }
    // Choose that restriction from 2-simplex or higher is identity function
    return 1.0;
}

petls::sheaf_simplex_tree get_sst(std::vector<std::vector<float>> points, std::vector<float> charges){
    petls::sheaf_simplex_tree sst = petls::rips_sheaf_simplex_tree(points, 6, my_restriction);

    // Add charge info to extra_data
    // You can add extra_data after the sst is constructed
    int counter = 0;
    for (Simplex_tree::Vertex_handle v : sst.st.complex_vertex_range()){
        sst.extra_data[v].push_back(charges[counter]);
        counter++;
    }
    return sst;
}

bool test_sst(std::vector<std::vector<float>> points, std::vector<float> charges){
    petls::sheaf_simplex_tree sst = get_sst(points, charges);
    
  
    std::vector<petls::FilteredBoundaryMatrix<float>> fbms = sst.apply_restriction_function();
    bool passed = true;
    float tol = 1e-3;
    
    // Note: this follows filtration-lexicographic ordering, ex 2.5 in paper does not
                            //v_0           v_1       v_2
    DenseMatrix_PL d1_dense {{-charges[1], charges[0],0}, // v_0 v_1
                            {-charges[2], 0, charges[0]}, // v_0 v_2
                            {0, -charges[2], charges[1]}}; // v_1 v_2
    DenseMatrix_PL d2_dense {{charges[2], -charges[1], charges[0]}};
    SparseMatrix_PL d1 = d1_dense.sparseView();
    SparseMatrix_PL d2 = d2_dense.sparseView();
    std::vector<SparseMatrix_PL> reference_boundaries = {d1, d2};
    for (int i = 0; i <  (int) fbms.size(); i++){
        SparseMatrix_PL d;
        fbms[i].submatrix_at_filtration(2,d);
        passed = passed && d.isApprox(reference_boundaries[i],tol);
    }
    return passed;
}

bool test_psl(std::vector<std::vector<float>> points, std::vector<float> charges){
    petls::sheaf_simplex_tree sst = get_sst(points, charges);
    petls::NewPersistentSheafLaplacian psl = petls::NewPersistentSheafLaplacian(sst);
    float q0 = charges[0];
    float q1 = charges[1];
    float q2 = charges[2];
    float s = q0*q0 + q1*q1 + q2*q2;
    std::vector<DenseMatrix_PL> expected_Laplacians = {
        DenseMatrix_PL({{q1*q1+q2*q2, -q0*q1, -q0*q2}, {-q0*q1, q0*q0+q2*q2, -q1*q2}, {-q0*q2, -q1*q2, q0*q0+q1*q1}}),
        DenseMatrix_PL({{s, 0, 0}, {0, s, 0}, {0, 0, s}}),
        DenseMatrix_PL({{s}})
    };

    bool passed = true;
    float tol = 1e-3;
    for (int i = 0; i < (int) expected_Laplacians.size(); i++){
        DenseMatrix_PL L;
        psl.get_L(i, 1, 1, L);
        passed = passed && L.isApprox(expected_Laplacians[i], tol);
    }
    if (!passed) return passed;
    std::vector<std::vector<spectra_type>> expected_spectra = {{0, 38, 38}, {38,38,38}, {38},{} };
    
    std::vector< std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra = psl.spectra();
    for (int i = 0; i < (int) expected_spectra.size(); i++){
        passed = passed && test_sample(psl, expected_spectra[i], i, 1.0, 1.0);
        // std::vector<spectra_type> s = expected_spectra[i];
        // std::vector<spectra_type> observed = std::get<3>(spectra[i]);
        // spectra_vec s_eigen = spectra_vec(s);
        // spectra_vec observed_eigen = spectra_vec(observed);
        // passed = passed && s_eigen.isApprox(observed_eigen);
        if (!passed){
            std::cout << "i = " << i <<" passed = false!" << std::endl;
        }
    }
    return true;
}

int main(){
    std::vector<std::vector<float>> points = {{0.0, 0.0, 0.0},//equilateral all sides length 1
    {1.0, 0.0, 0.0},
    {0.5, (float)( sqrt(3.0)/2.0), 0.0}};
    std::vector<float> charges = {2.0, 3.0, 5.0};
    bool passed = test_sst(points, charges);
    passed = passed && test_psl(points, charges);
    if (passed){
        return 0;
    } else{
        return -1;
    }
}