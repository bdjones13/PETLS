#ifndef SST_H
#define SST_H
// wrap gudhi simplex tree with extra data
#include "FilteredBoundaryMatrix.hpp"


#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <vector>
#include "Eigen/SparseCore"
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace petls{
    using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_default>;
    using Filtration_value = Simplex_tree::Filtration_value;
    using Distance_matrix = std::vector<std::vector<Filtration_value>>;
    class sheaf_simplex_tree{
        public:
            Simplex_tree st;
            std::unordered_map<Simplex_tree::Simplex_key, std::vector<float>> extra_data;
            std::function<float(Simplex_tree::Simplex_handle, Simplex_tree::Simplex_handle, sheaf_simplex_tree&)> restriction;
            sheaf_simplex_tree(Simplex_tree _st,
                            std::unordered_map<Simplex_tree::Simplex_key, std::vector<float>> _extra_data,
                            std::function<float(Simplex_tree::Simplex_handle, Simplex_tree::Simplex_handle, sheaf_simplex_tree&)> _restriction){
                st = Simplex_tree(_st);
                extra_data = _extra_data;
                restriction = _restriction;
            }
            int coface_index(Simplex_tree::Simplex_handle simplex, Simplex_tree::Simplex_handle coface){
                // get iterator for vertices of simplex and coface
                Simplex_tree::Simplex_vertex_range vertex_range_simplex = this->st.simplex_vertex_range(simplex);
                Simplex_tree::Simplex_vertex_range vertex_range_coface = this->st.simplex_vertex_range(coface);
                
                // convert to vector
                std::vector<Simplex_tree::Vertex_handle> vec_simplex(vertex_range_simplex.begin(), vertex_range_simplex.end()); 
                std::vector<Simplex_tree::Vertex_handle> vec_coface(vertex_range_coface.begin(), vertex_range_coface.end()); 
                
                // simplex_vertex_range returns simplices in reverse order
                std::reverse(vec_simplex.begin(), vec_simplex.end());
                std::reverse(vec_coface.begin(), vec_coface.end());
                
                // loop until no match
                for (int i = 0; i < (int) vec_simplex.size(); i++){
                    if (vec_simplex[i] != vec_coface[i]) 
                        return i;
                }
                return vec_coface.size() - 1; // missing vertex is in last position 
            }

            std::vector<FilteredBoundaryMatrix<float>> apply_restriction_function(){
                // std::function<float(Simplex_tree::Simplex_handle, Simplex_tree::Simplex_handle, sheaf_simplex_tree&)> f){
                int complex_dim = this->st.dimension();
                std::vector<size_t> dims = this->st.num_simplices_by_dimension();

                // intermediate sparse matrix storage as (row, col, coefficient)
                // this is later reindexed then given to Eigen to create a SparseMatrixFloat (MatrixXf)
                std::vector<std::vector<std::tuple<int64_t,int64_t,float>>> boundaries_triples;
                std::vector<std::vector<std::pair<int,double>>> filtrations;

                std::vector<SparseMatrixFloat> coboundary_matrices;
                std::vector<std::vector<double>> reindexed_filtrations;
                
                // Wrapper to keep explicit boundary matrices with the filtration of each simplex 
                std::vector<FilteredBoundaryMatrix<float>> coboundaries;

                // Initialize vectors for filtrations and boundary matrices
                for (int dim = 0; dim < (int) complex_dim; dim++){
                    filtrations.push_back({});
                    boundaries_triples.push_back({});
                }
                filtrations.push_back({}); //for top dim
                
                // loop over all simplices

                for (auto sh : st.filtration_simplex_range()){
                    int dim = st.dimension(sh);
                    filtrations[dim].push_back(std::make_pair<int, float>(st.key(sh),st.filtration(sh)));
                    if (dim == complex_dim){ // zero coboundary from top dimension
                        continue;
                    }

                    // loop over all codimension-1 cofaces
                    for (Simplex_tree::Simplex_handle c : st.cofaces_simplex_range(sh, 1)){
                        
                        float sign = (float) std::pow(-1.0, (double) (this->coface_index(sh, c) % 2));
                        
                        // apply the restriction function
                        float coeff = sign*restriction(sh, c, *this);

                        // store as (row, col, coefficient)
                        boundaries_triples[dim].push_back(std::tuple<int64_t,int64_t,float>(st.key(c),st.key(sh),coeff));
                    }
                }
            reindex_boundaries_map(boundaries_triples,coboundary_matrices,filtrations, reindexed_filtrations);

                // convert from reindexed (row, col, coeff) to FilteredBoundaryMatrix (wraps Eigen::MatrixXf)
            for (int dim = 0; dim < (int) complex_dim; dim++){
                SparseMatrixFloat d = coboundary_matrices[dim];
                std::vector<double> filt_domain = reindexed_filtrations[dim];
                std::vector<double> filt_range = reindexed_filtrations[dim+1];
                FilteredBoundaryMatrix<float> fbm(coboundary_matrices[dim],filt_domain, filt_range);
                coboundaries.push_back(fbm);
            }
            return coboundaries;
        }
        private:
            void reindex_boundaries_map(std::vector<std::vector<std::tuple<int64_t,int64_t,float>>> &boundaries_triples,std::vector<SparseMatrixFloat> &reindexed_boundaries, 
                            std::vector<std::vector<std::pair<int,double>>> &filtrations, std::vector<std::vector<double>> &reindexed_filtrations){
                std::vector<std::set<int64_t>> indices_of_actual_simplices_set(boundaries_triples.size()+1);
                std::vector<std::vector<int>> indices_of_actual_simplices_vector(boundaries_triples.size()+1);
                
                std::vector<std::unordered_map<int64_t,int>> index_maps(filtrations.size());
	     
                for (int i = 0; i < (int) filtrations.size(); i++){
                    // int dim = i+1;
                    std::vector<std::pair<int,double>> filtration_pairs_dim = filtrations[i]; // starts with dim 1 simplices 
                
                    // sort by filtration then by simplex number
                    std::sort(filtration_pairs_dim.begin(), filtration_pairs_dim.end(), [](auto &left, auto &right){
                        return (left.first < right.first) || (left.first == right.first && left.second < right.second);
                    });
                
                    // build mappings
                    std::unordered_map<int64_t, int> temp_map(filtration_pairs_dim.size());
                    std::vector<double> filtrations_dim;
                    for (int j = 0; j < (int) filtration_pairs_dim.size(); j++){
                        temp_map[filtration_pairs_dim[j].first] = j;
                        filtrations_dim.push_back(filtration_pairs_dim[j].second);
                    }
                    reindexed_filtrations.push_back(filtrations_dim);
                    index_maps[i] = temp_map;
                }

                for (int i = 0; i < (int) boundaries_triples.size(); i++){
                    int dim = i + 1;
                    std::vector<Eigen::Triplet<float>> tripletList;
                    std::vector<std::tuple<int64_t, int64_t, float>> boundary_triples = boundaries_triples[i];
                    for (int j = 0; j < (int) boundary_triples.size(); j++){
                        int row = index_maps[dim][std::get<0>(boundary_triples[j])];
                        int col = index_maps[dim-1][std::get<1>(boundary_triples[j])];
                        float coeff = std::get<2>(boundary_triples[j]);
                        boundary_triples[j] = std::make_tuple(row, col, coeff);
                        tripletList.push_back(Eigen::Triplet<float>(row, col, coeff));
                    }
                    // boundaries_triples[i] = boundary_triples;
                    SparseMatrixFloat boundary(index_maps[dim].size(), index_maps[dim-1].size()); 
            
                    boundary.setFromTriplets(tripletList.begin(), tripletList.end());
                    reindexed_boundaries.push_back(boundary);
                }
            
            }

    };

    petls::sheaf_simplex_tree rips_sheaf_simplex_tree(std::vector<std::vector<float>> points, Filtration_value max_length,
        std::function<float(Simplex_tree::Simplex_handle, Simplex_tree::Simplex_handle, sheaf_simplex_tree&)> restriction){
        using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

        // This is an example of a function that constructs a sheaf simplex tree
        
        
        // Important: We must wrap the Gudhi Rips_complex rips in our final sheaf_simplex_tree sst 
        //      BEFORE calling rips.create_complex(sst.st, dim_max).
        //      This is because of how Gudhi copies certain information. 
        
        // Declare a rips complex with the given point cloud, 
        //  but do NOT call rips.create_complex yet.
        Rips_complex rips(points,max_length,Gudhi::Euclidean_distance());
        
        // Declare but do NOT give a value to a Simplex_tree
        petls::Simplex_tree st;
        int dim_max = 3;

        // Extra data is essentially a dictionary with keys of type Simplex_key and values of vector<float>
        std::unordered_map<petls::Simplex_tree::Simplex_key, std::vector<float>> extra_data; 
        
        // Construct sst as wrapping the declared simplex tree and extra data 
        petls::sheaf_simplex_tree sst(st, extra_data, restriction);

        // Store the rips complex in the simplex tree that has been wrapped.
        rips.create_complex(sst.st, dim_max);

        // Assuming that Gudhi's rips complex preserves the order of points
        int counter = 0;

        // Assign vertex keys and add (x,y,z) coordinates to extra_data
        for (Simplex_tree::Vertex_handle v : sst.st.complex_vertex_range()){
            Simplex_tree::Simplex_handle as_sh = sst.st.find({v});
            sst.st.assign_key(as_sh, counter);
            sst.extra_data[sst.st.key(as_sh)] = {points[counter][0],
                points[counter][1],
                points[counter][2]};
            counter++;
        }
        // Assign the dimension >= 1 simplex keys
        for (Simplex_tree::Simplex_handle sh : sst.st.filtration_simplex_range()){
            if (sst.st.dimension(sh) == 0)
                continue;
            sst.st.assign_key(sh, counter++);
        }
        return sst;
    }


    

}

#endif