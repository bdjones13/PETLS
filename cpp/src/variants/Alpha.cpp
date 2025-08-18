#include "../../include/petls_headers/variants/Alpha.hpp"
#include <vector>
#include <set>
#include <iostream>
#ifdef PETLS_USE_ALPHA_COMPLEX



// Imports and typedefs for Gudhi
#include <gudhi/Simplex_tree.h>
#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Points_3D_off_io.h>

using Alpha_complex_3d = Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
using Point = Alpha_complex_3d::Point_3;
using Vector_of_points = std::vector<Point>;
using Simplex_tree = Gudhi::Simplex_tree<>;


namespace petls{
    void alpha_points(std::vector<Point> &points, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations);
    void alpha_OFF(const char* filename, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations);
            

    Alpha::Alpha(const char* filename, int max_dim) : petls::Complex() { // call default constructor of parent class

        std::vector<std::vector<std::tuple<int,int,int>>> boundaries_triples(max_dim);
        std::vector<std::vector<filtration_type>> filtrations(max_dim+1);
        verbose = false;
        alpha_OFF(filename, max_dim, boundaries_triples, filtrations);
        std::vector<SparseMatrixInt> reindexed_boundaries(boundaries_triples.size());
        reindex_boundaries(boundaries_triples, reindexed_boundaries);
        set_boundaries_filtrations(reindexed_boundaries, filtrations);

    }

    Alpha::Alpha(std::vector<std::tuple<double,double,double>> points, int max_dim) : petls::Complex() { // call default constructor of parent class
        std::vector<std::vector<std::tuple<int,int,int>>> boundaries_triples(max_dim);
        std::vector<std::vector<filtration_type>> filtrations(max_dim+1);
        verbose = false;
        std::vector<Point> points_gudhi(points.size());
        // std::cout << "number of points:" << points.size() << std::endl;
        for (int i = 0; i <  (int) points.size(); i++){
            points_gudhi[i] = Point(std::get<0>(points[i]),
                                    std::get<1>(points[i]),
                                    std::get<2>(points[i])
            );
        }
        alpha_points(points_gudhi, max_dim, boundaries_triples, filtrations);
        // alpha_OFF(filename, max_dim, boundaries_triples, filtrations);
        std::vector<SparseMatrixInt> reindexed_boundaries(boundaries_triples.size());
        reindex_boundaries(boundaries_triples, reindexed_boundaries);
        set_boundaries_filtrations(reindexed_boundaries, filtrations);
    }



    void get_boundaries_and_filtrations(Simplex_tree simplex_tree, int dim_max,  std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        // give each simplex a unique index (assign_key)
        // method from Gudhi rips_persistence_via_boundary_matrix.cpp
        int count = 0;
        for (auto simplex_handle : simplex_tree.filtration_simplex_range()){
            simplex_tree.assign_key(simplex_handle, count++);
        }

        // int max_dim = 3;

        //initialize list of boundaries
        for (int dim = 0; dim < dim_max; dim++){
            boundaries_triples[dim] = std::vector<std::tuple<int,int,int>>();
        }

        // initialize list of filtrations
        for (int dim = 0; dim <= dim_max; dim++){
            filtrations[dim] = std::vector<filtration_type>();
        }

        // iterate over all simplices to store their boundary information
        for (auto f_simplex : simplex_tree.filtration_simplex_range())
        {
            int dim_f = simplex_tree.dimension(f_simplex);
            if (dim_f > dim_max) continue; // skip simplices of dimension greater than max_dim
            // store filtration
            int f_key = simplex_tree.key(f_simplex);
            filtrations[dim_f].push_back(simplex_tree.filtration(f_simplex));

            double f_filtration_value = simplex_tree.filtration(f_simplex);
            
            int sign = 1 - 2 * (dim_f % 2); // from Gudhi Persistent_cohomology.h
            for (auto b_simplex : simplex_tree.boundary_simplex_range(f_simplex))
            {

                boundaries_triples[dim_f-1].push_back(std::make_tuple(simplex_tree.key(b_simplex),
                                                                    f_key,
                                                                    sign));
                sign = -sign;
            }
        }
        
    }

    void alpha_points(std::vector<Point>& points, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        Alpha_complex_3d alpha_complex_from_points(points);
        Gudhi::Simplex_tree<> simplex_tree;
        alpha_complex_from_points.create_complex(simplex_tree);
        
        get_boundaries_and_filtrations(simplex_tree, dim_max, boundaries_triples, filtrations);
        
    }
    

    void alpha_OFF(const char* filename, int dim_max, std::vector<std::vector<std::tuple<int,int,int>>>& boundaries_triples,std::vector<std::vector<filtration_type>>& filtrations){
        // Read points from file
        std::string offInputFile(filename);
        // Read the OFF file (input file name given as parameter) and triangulate points
        Gudhi::Points_3D_off_reader<Point> off_reader(offInputFile);
        // Check the read operation was correct
        if (!off_reader.is_valid())
        {
            std::cerr << "Unable to read file " << filename << std::endl;
            return;
        }
        // Retrieve the triangulation
        std::vector<Point> points = off_reader.get_point_cloud();
        alpha_points(points, dim_max, boundaries_triples, filtrations);
    }
}



#else // PETLS_USE_ALPHA_COMPLEX not defined
namespace petls{
    Alpha::Alpha(const char* filename, int max_dim){
        std::cout << "PETLS_USE_ALPHA_COMPLEX not defined; Alpha complex not being used." << std::endl;
    } 
    Alpha::Alpha(std::vector<std::tuple<double,double,double>> points, int max_dim){
        std::cout << "PETLS_USE_ALPHA_COMPLEX not defined; Alpha complex not being used." << std::endl;
    } 
}
#endif // PETLS_USE_ALPHA_COMPLEX not defined