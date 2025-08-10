#ifndef PLs_H
#define PLs_H

#include "typedefs.hpp"
#include <string>
#include <vector>
#include <chrono>

namespace petls {

    // https://stackoverflow.com/questions/54503795/rounding-numbers-below-a-certain-threshold-to-zero-in-eigen
    void round_zeros(spectra_vec& inout, spectra_type threshold);
    void print_full_matrix_precise(SparseMatrixInt m);
    void print_full_matrix_precise(Eigen::MatrixXf m);

    void print_full_matrix_precise(Eigen::MatrixXi m);
    // void print_full_matrix_precise(SparseMatrixInt m);
    void print_full_matrix_precise(SparseMatrix_PL m);
    void print_full_matrix_precise(DenseMatrix_PL m);
    void print_vector_precise(Eigen::VectorXd v);
    void print_vector_precise(Eigen::VectorXf v);
    void print_spectra(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra);

    void reindex_boundaries(std::vector<std::vector<std::tuple<int,int,int>>> &boundaries_triples,std::vector<SparseMatrixInt> &reindexed_boundaries);
    struct timer{
        public:
            int duration;
            std::chrono::_V2::steady_clock::time_point start_t;
            void start(){
                start_t = std::chrono::steady_clock::now();
            }
            void stop(){
                auto end_t = std::chrono::steady_clock::now();
                duration = (std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t)).count();
            }
    };

    struct Profile{
        public:
            std::vector<int> dims;
            std::vector<filtration_type> filtration_a;
            std::vector<filtration_type> filtration_b;
            timer all;
            timer eigs;
            timer L;
            timer L_up;
            timer L_down;
            timer sum_up_down;
            std::vector<int> durations_all;
            std::vector<int> durations_eigs;
            std::vector<int> durations_L;
            std::vector<int> durations_L_up;
            std::vector<int> durations_L_down;
            std::vector<int> durations_sum_up_down;
            
            std::vector<int> L_rows;
            std::vector<int> bettis;
            std::vector<spectra_type> lambdas;
            // std::vector<double> sparsity_L;
            // std::vector<int> simplices_added;

            void to_csv(std::string filename);
            
            void start_all(){this->all.start();}
            void start_eigs(){eigs.start();}
            void start_L(){L.start();}
            void start_L_up(){L_up.start();}
            void start_L_down(){L_down.start();}
            void start_sum_up_down(){sum_up_down.start();}

            void stop_all(){all.stop(); durations_all.push_back(all.duration);}
            void stop_eigs(){eigs.stop(); durations_eigs.push_back(eigs.duration);}
            void stop_L(){L.stop(); durations_L.push_back(L.duration);}
            void stop_L_up(){L_up.stop(); durations_L_up.push_back(L_up.duration);}
            void stop_L_down(){L_down.stop(); durations_L_down.push_back(L_down.duration);}
            void stop_sum_up_down(){sum_up_down.stop(); durations_sum_up_down.push_back(sum_up_down.duration);}

    };
}

#endif