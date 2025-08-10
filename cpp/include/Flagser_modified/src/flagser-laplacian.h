#ifndef FLAGSER_H
#define FLAGSER_H
#include <fstream>
#include <iostream>

#define ASSEMBLE_REDUCTION_MATRIX
// #define INDICATE_PROGRESS
#define SKIP_APPARENT_PAIRS
#define USE_ARRAY_HASHMAP
#define USE_CELLS_WITHOUT_DIMENSION
#define SORT_COLUMNS_BY_PIVOT


#include "../include/persistence_real.h"
#include "../include/filtration_algorithms.h"
#include "../include/input/flagser.h"
#include "../include/complex/directed_flag_complex_computer_real.h"

// #include "../include/usage/flagser.h"

template <class T> void compute_homology(filtered_directed_graph_t& graph,
										// std::unique_ptr<filtration_algorithm_t> filtration_algorithm, 
										unsigned short min_dim, unsigned short max_dim, 
										std::vector<Eigen::SparseMatrix<int,Eigen::ColMajor>>& boundaries, 
										std::vector<std::vector<double>>& filtrations) {
	// std::unique_ptr<filtration_algorithm_t> filtration_algorithm;
	// filtration_algorithm.reset(get_filtration_computer("max"));
	// filtration_algorithm.reset(get_filtration_computer("max"));

	filtration_algorithm_t *filtration_algorithm;
	filtration_algorithm = get_filtration_computer("max");

	T complex(graph,filtration_algorithm, min_dim, max_dim);
	real_persistence_computer_t<decltype(complex)> persistence_computer(complex);
	persistence_computer.compute_persistent_spectra(min_dim, max_dim, boundaries, filtrations);
}


void basic_flagser(const char* filename, int min_dim, int max_dim,
					std::vector<Eigen::SparseMatrix<int,Eigen::ColMajor>>& boundaries,
					std::vector<std::vector<double>>& filtrations){
	filtered_directed_graph_t graph = read_graph_flagser(filename);
	
	// std::unique_ptr<filtration_algorithm_t> filtration_algorithm;
	// filtration_algorithm.reset(get_filtration_computer("max"));
	
	compute_homology<real_directed_flag_complex_computer::real_directed_flag_complex_computer_t>(graph, min_dim, max_dim, boundaries, filtrations);
}

#endif //FLAGSER_H