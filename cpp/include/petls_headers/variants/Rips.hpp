#ifndef PL_R_H
#define PL_R_H

#include "../core/Complex.hpp"
#include "../core/up_algorithms.hpp"
#include "../eigenvalues/eigs_algorithms.hpp"
#include "../petls.hpp"
#include "../../include/Ripser_modified/ripser_modified.h"


namespace petls{
    class Rips : public Complex {
        public:

            // Lower_Distance_Matrix file format
            Rips(const char* filename, int max_dim) : petls::Complex() { // call default constructor of parent class
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                float threshold = std::numeric_limits<value_t>::max();
                rips_distance_matrix_file(filename, (long long int) max_dim, threshold, boundaries_triples, filtrations);
                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);
                this->set_boundaries_filtrations(boundaries, filtrations);
                
            }

            // Lower_Distance_Matrix file format 
            Rips(const char* filename, int max_dim, float threshold) : petls::Complex() { // call default constructor of parent class
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                rips_distance_matrix_file(filename, (long long int) max_dim, threshold, boundaries_triples, filtrations);

                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);
                this->set_boundaries_filtrations(boundaries, filtrations);
                
            }

            // Eigen DenseMatrix with distances (correlations) in lower triangular part
            Rips(DenseMatrix_PL distances, int max_dim) : petls::Complex() { // call default constructor of parent class
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                int num_rows = distances.rows();
                std::vector<value_t> distances_vec;
                distances_vec.reserve(num_rows*(num_rows-1)/2);
                for (int i = 1; i < num_rows; i++){    
                    for (int j = 0; j < i; j++){        
                        distances_vec.push_back(distances(i,j));        
                    }   
                }
                compressed_lower_distance_matrix ripser_mat(std::move(distances_vec));

                basic_rips(ripser_mat,max_dim, std::numeric_limits<value_t>::infinity(), boundaries_triples, filtrations);

                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);

                this->set_boundaries_filtrations(boundaries, filtrations);
            }

            // Eigen DenseMatrix with distances (correlations) in lower triangular part
            Rips(DenseMatrix_PL distances, int max_dim, float threshold) : petls::Complex() { // call default constructor of parent class
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                int num_rows = distances.rows();
                std::vector<value_t> distances_vec;
                distances_vec.reserve(num_rows*(num_rows-1)/2);
                for (int i = 1; i < num_rows; i++){
                    for (int j = 0; j < i; j++){
                        distances_vec.push_back(distances(i,j));
                    }
                }
                compressed_lower_distance_matrix ripser_mat(std::move(distances_vec));

                basic_rips(ripser_mat,max_dim, threshold, boundaries_triples, filtrations);

                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);

                this->set_boundaries_filtrations(boundaries, filtrations);
            }


            // Point cloud format
            Rips(std::vector<std::vector<double>> points, int max_dim){
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                float threshold = std::numeric_limits<value_t>::max();
                rips_point_cloud(points, max_dim, threshold, boundaries_triples, filtrations);

                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);
                
                this->set_boundaries_filtrations(boundaries, filtrations);
                
            }

            Rips(std::vector<std::vector<double>> points, int max_dim, float threshold){
                std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> boundaries_triples;
                std::vector<std::vector<filtration_type>> filtrations;
                rips_point_cloud(points, max_dim, threshold, boundaries_triples, filtrations);

                std::vector<SparseMatrixInt> boundaries(boundaries_triples.size());
                boundaries_eigen(boundaries_triples, boundaries, filtrations);
                this->set_boundaries_filtrations(boundaries, filtrations);
                
            }

            // These static methods are to unambiguously wrap with Pybind11
            // NumPy arrays can get converted to Eigen::Matrix or std::vector<std::vector< >>
            static Rips from_distances(DenseMatrix_PL distances, int max_dim){
                return Rips(distances, max_dim);
            }

            static Rips from_distances(DenseMatrix_PL distances, int max_dim, float threshold){
                return Rips(distances, max_dim, threshold);
            }

            static Rips FromPoints(std::vector<std::vector<double>> points, int max_dim){
                return Rips(points, max_dim);
            }
            static Rips FromPoints(std::vector<std::vector<double>> points, int max_dim, float threshold){
                return Rips(points, max_dim, threshold);
            }


            private:
            void boundaries_eigen(std::vector<std::vector<std::tuple<int64_t,int64_t,int>>> &boundaries_triples, std::vector<SparseMatrixInt> &boundaries, std::vector<std::vector<value_t>>& filtrations){
                for (int i = 0; i < (int) boundaries_triples.size(); i++){
                    
                    int dim = i + 1;
                    std::vector<Eigen::Triplet<int>> tripletList;
                    std::vector<std::tuple<int64_t, int64_t, int>> boundary_triples = boundaries_triples[i]; // is this an unnecessary copy?
                    // convert to eigen triple format and use the mapping
                    for (int j = 0; j < (int) boundary_triples.size(); j++){
                        int row = std::get<0>(boundary_triples[j]);
                        int col = std::get<1>(boundary_triples[j]);
                        int coeff = std::get<2>(boundary_triples[j]);
                        tripletList.push_back(Eigen::Triplet<int>(row, col, coeff));
                    }
                    SparseMatrixInt boundary(filtrations[dim-1].size(), filtrations[dim].size()); 
        
                    boundary.setFromTriplets(tripletList.begin(), tripletList.end());
                    boundaries[i]=boundary;
                }
            }
    };
}

#endif