#include "../include/petls_headers/petls.hpp"
#include <string>
#include <vector>
#include <chrono>

#include <iostream>
#include <unordered_map>
#include <fstream>
#include <set>

namespace petls {

    // https://stackoverflow.com/questions/54503795/rounding-numbers-below-a-certain-threshold-to-zero-in-eigen
    void round_zeros(spectra_vec& inout, spectra_type threshold) {
        inout = (threshold < inout.array().abs()).select(inout, 0.0);
    }
    
    void print_full_matrix_precise(SparseMatrixInt m){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << Eigen::MatrixXi(m).format(HeavyFmt) <<std::endl;
    }

    void print_full_matrix_precise(Eigen::MatrixXf m);

    void print_full_matrix_precise(Eigen::MatrixXi m){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << m.format(HeavyFmt) <<std::endl;
    }

    // void print_full_matrix_precise(SparseMatrixInt m);
    void print_full_matrix_precise(SparseMatrix_PL m){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << DenseMatrix_PL(m).format(HeavyFmt) <<std::endl;
    }

    void print_full_matrix_precise(DenseMatrix_PL m){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << m.format(HeavyFmt) <<std::endl;
    }
    void print_vector_precise(Eigen::VectorXd v){
        //TODO: maybe remove newlines?
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << v.format(HeavyFmt) <<std::endl;
    }

    void print_vector_precise(Eigen::VectorXf v){
        //TODO: maybe remove newlines?
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout << v.format(HeavyFmt) <<std::endl;
    }

    void print_spectra(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra){
        for (int i = 0; i < (int) spectra.size(); i++){
            int dim = std::get<0>(spectra[i]);
            filtration_type a = std::get<1>(spectra[i]);
            filtration_type b = std::get<2>(spectra[i]);
            std::vector<spectra_type> spectrum = std::get<3>(spectra[i]);
            std::cout << "spectra (dim, a, b) = (" << dim << ", " << a << ", " << b << ") = [";
            for (int j = 0; j < (int) spectrum.size(); j++){
                std::cout << spectrum[j] << " ";
            }
            std::cout << "]\n";
        }
    }


    void reindex_boundaries(std::vector<std::vector<std::tuple<int,int,int>>> &boundaries_triples,std::vector<SparseMatrixInt> &reindexed_boundaries){
        std::vector<SparseMatrixInt> boundaries;
        std::vector<std::set<int>> indices_of_actual_simplices_set(boundaries_triples.size()+1);
        std::vector<std::vector<int>> indices_of_actual_simplices_vector(boundaries_triples.size()+1);
        for (int i = 0; i < (int) boundaries_triples.size(); i++){
            int dim = i + 1;
            std::vector<Eigen::Triplet<int>> tripletList;
            
            std::vector<std::tuple<int, int, int>> boundary_triples = boundaries_triples[i]; // is this an unnecessary copy?

            tripletList.reserve(boundary_triples.size());
            int max_ind_row = 0;
            int max_ind_col = 0;
            for (int j = 0; j < (int) boundary_triples.size(); j++){
                int row = std::get<0>(boundary_triples[j]);
                int col = std::get<1>(boundary_triples[j]);
                int coeff = std::get<2>(boundary_triples[j]);
                // std::cout << "triplet (" << row << ", " << col << ", " << coeff << ")" << std::endl;
                tripletList.push_back(Eigen::Triplet<int>(row, col, coeff));
                indices_of_actual_simplices_set[dim-1].insert(row);
                indices_of_actual_simplices_set[dim].insert(col);
                if (row > max_ind_row) max_ind_row = row;
                if (col > max_ind_col) max_ind_col = col;
            }

                SparseMatrixInt boundary(max_ind_row+1,max_ind_col+1);

                boundary.setFromTriplets(tripletList.begin(), tripletList.end());

                boundaries.push_back(boundary);

        }
        
        for (int i = 0; i < (int) indices_of_actual_simplices_set.size(); i++){
            std::vector<int> temp_indices;
            temp_indices.assign(indices_of_actual_simplices_set[i].begin(), indices_of_actual_simplices_set[i].end());
            indices_of_actual_simplices_vector[i] = temp_indices;           
        }

        // TODO: add assertions on sizes of these vectors to make sure they are consistent
      
        for (int i = 0; i < (int) boundaries.size(); i++){
          
            Eigen::MatrixXi current = Eigen::MatrixXi(boundaries[i]); // IS THIS SPARSE?
           
            Eigen::MatrixXi removed_cols = current(Eigen::indexing::all,
                indices_of_actual_simplices_vector[i+1]);
            
            Eigen::MatrixXi removed_cols_and_rows = removed_cols(indices_of_actual_simplices_vector[i],
                Eigen::indexing::all);
            
            reindexed_boundaries[i]=removed_cols_and_rows.sparseView();

        }
    }

    void Profile::to_csv(std::string filename){
        std::ofstream out_file(filename);
        out_file << "dim,filtration_a,filtration_b,duration_all,duration_eigs,duration_L,duration_L_up,duration_L_down,duration_sum_up_down,L_rows,betti,lambda" << std::endl;
        for (int i = 0; i < (int) dims.size(); i++){
            out_file << dims[i] << "," << filtration_a[i] << "," << filtration_b[i] << ",";
            out_file << durations_all[i] << "," << durations_eigs[i] << "," << durations_L[i] << "," << durations_L_up[i] << ",";
            out_file << durations_L_down[i] << "," << durations_sum_up_down[i] << "," << L_rows[i] << ",";
            if ((int)bettis.size() == 0) // if no betti numbers/lambdas were computed, output flag=-1
                out_file << -1 << ",";
            else
                out_file << bettis[i] << ",";
            if ((int)lambdas.size() == 0)
                out_file << -1;
            else
                out_file << lambdas[i];
            out_file << std::endl;             
        }
        out_file.close();
    }

}