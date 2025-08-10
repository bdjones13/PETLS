#ifndef PL_H
#define PL_H

#include <vector>

#include "../typedefs.hpp"
#include "../petls.hpp"
// #include "up_algorithms.hpp"
// #include "../eigenvalues/eigs_algorithms.hpp"
#include "FilteredBoundaryMatrix.hpp"

// #include <unsupported/Eigen/SparseExtra>


// #include <chrono>
// #include <numeric> //std::iota
// #include <map>
// #include <set>

// #include <cassert>

namespace petls{
/**
 * Primary class for computing persistent Laplacians
 * @tparam eigs_Algorithm Algorithm wrapper class to use for computing eigenvalues of the Laplacian in the spectra family of functions. Default is a wrapper for Eigen::SelfAdjointEigenSolver.
 * @tparam up_Algorithm Algorithm wrapper class to use for computing the up-Laplacian. Default is the Schur complement algorithm presented in Memoli, Wan, and Wang 2020.
 * //TODO: \\see up_algorithms and eigs_algorithms 
*/
// template <typename storage = int>
class Complex {
    using storage = int;
    public:
        /********************/
        /* Member variables */
        /********************/

        int top_dim; ///< Top dimension of the complex
        std::vector<FilteredBoundaryMatrix<storage>> filtered_boundaries; ///< Boundary matrix assuming real (or integer) coefficients
        bool verbose; ///< Print progress if spectra() is called
        bool use_flipped; ///< Compute the top-dimensional Laplacian's eigenvalues in spectra function via the eigenvalues of the smaller of B_N B_N^T or B_N^T B_N and possible zero-padding
        Profile profile; ///< Profiler to track time usage of various steps
        std::function<spectra_vec(DenseMatrix_PL&)> eigs_algorithm_func;
        std::function<std::pair<spectra_vec,DenseMatrix_PL>(DenseMatrix_PL&)> eigenpairs_algorithm_func;
        // std::function<void(FilteredBoundaryMatrix<storage>*, filtration_type, filtration_type, DenseMatrix_PL&)> up_algorithm_func;
        std::function<void(FilteredBoundaryMatrix<storage>*, filtration_type, filtration_type, DenseMatrix_PL&)> up_algorithm_func;        
        //may need to be void* instead of FilteredBoundaryMatrix*

        /****************/
        /* Constructors */
        /****************/

        /**
         * Default constructor with no boundary maps or simplices.
         */
        Complex();
        /**
         * @param boundaries a vector of Eigen::SparseMatrix of type int. Boundaries must be sorted in order of dimension.
         * @param filtrations a vector of vector of filtrations. filtrations[dim] is a list of all filtrations of simplices in dimension dim. filtrations must be sorted in order of dimension, each filtrations[i] must be sorted in order of filtration
         *
         * Primary constructor. Important Assumptions:
         * 1) Boundary matrix has real coefficients stored as integers (but not mod 2!)
         * 2) Boundary matrix dimensions agree with filtrations sizes
         * 3) Length(filtrations) = Length(boundaries) + 1
         */
        Complex(std::vector<Eigen::SparseMatrix<storage>> boundaries,
                            std::vector<std::vector<filtration_type>> filtrations);
        /****************** */
        /* Empty destructor */
        /****************** */

        ~Complex(){}

        /***********************/
        /* Getters and setters */
        /***********************/

        void set_eigs_algorithm_func(std::function<spectra_vec(DenseMatrix_PL&)> _eigs_algorithm_func);
        void set_eigs_algorithm_func(std::string name);
        void set_up_func(std::function<void(FilteredBoundaryMatrix<int>*, filtration_type, filtration_type, DenseMatrix_PL&)> _up_algorithm_func);
        /**
         * Set the boundaries and filtrations of a complex, particularly if the default constructor was called.
         * @param boundaries a vector of Eigen::SparseMatrix of type int. Boundaries must be sorted in order of dimension.
         * @param filtrations a vector of vector of filtrations. filtrations[dim] is a list of all filtrations of simplices in dimension dim. filtrations must be sorted in order of dimension, each filtrations[i] must be sorted in order of filtration.
         * Important Assumptions:
         * 1) Boundary matrix has real coefficients stored as integers (but not mod 2!),
         * 2) Boundary matrix dimensions agree with filtrations sizes,
         * 3) Length(filtrations) = Length(boundaries) + 1.
         */
        void set_boundaries_filtrations(std::vector<Eigen::SparseMatrix<storage>> boundaries,
                            std::vector<std::vector<filtration_type>> filtrations);
        
        /**
         * Set verbose
         * @param verbose New setting.
         */
        void set_verbose(bool verbose){this->verbose = verbose;}

        /**
         * Set flipped
         * @param use_flipped New setting
         */
        void set_flipped(bool use_flipped){this->use_flipped = use_flipped;}
        
        /***********************************/
        /* Primary Mathematical Operations */
        /***********************************/


        /**
         * Get the Persistent Laplacian Matrix (by reference).
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value (must be >= a)
         * @param[out] L the matrix where the persistent Laplacian will be stored (by reference).
         */
        void get_L(int dim, filtration_type a, filtration_type b, DenseMatrix_PL &L);
        /**
         * Get the Persistent Laplacian Matrix (by value). 
         * 
         * Warning: this does a potentially expensive copy. The primary reason this function exists is to provide reasonable python
         * binding access to the matrix L itself. If you do not need the persistent Laplacian matrix directly in python (e.g. you want its eigenvalues)
         * it will be more efficient to compute by reference via get_L(int dim, filtration_type a, filtration_type b, DenseMatrix_PL &L).
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value (must be >= a)
         * \return L the persistent Laplacian matrix.
         */
        DenseMatrix_PL get_L(int dim, filtration_type a, filtration_type b);
        /**
         * Get a matrix with the same nonzero eigenvalues as the top-dimensional Persistent Laplacian Matrix (by reference). 
         * @param a start filtration value
         * @param[out] L the matrix (by reference) that has the same eigenvalues as the top-dimensional Persistent Laplacian
         */
        void get_L_top_dim_flipped(filtration_type a, SparseMatrix_PL &L);
        /**
         * Get a matrix with the same nonzero eigenvalues as the top-dimensional Persistent Laplacian Matrix (by value). 
         * 
         * Warning: this does a potentially expensive copy. The primary reason this function exists is to provide reasonable python
         * binding access to the matrix itself. If you do not need this matrix directly in python (e.g. you want its eigenvalues)
         * it will be more efficient to compute by reference via get_L_top_dim_flipped(filtration_type a, SparseMatrix_PL &L).
         * @param a start filtration value
         * \return A matrix with the same nonzero eigenvalues as the top-dimensional persistent Laplacian (but not the persistent Laplacian itself).
         */
        SparseMatrix_PL get_L_top_dim_flipped(filtration_type a);

        /**
         * Get the up persistent Laplacian (by reference). The algorithm used is determined by the template parameter. See up_algorithms.hpp.
         * 
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value
         * @param[out] L_up up persistent Laplacian (by reference)
         */
        void get_up(int dim, filtration_type a, filtration_type b, DenseMatrix_PL &L_up);

        /**
         * Get the up persistent Laplacian (by value). The algorithm used is determined by the template parameter. See up_algorithms.hpp.
         * 
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value
         * \return up persistent Laplacian
         */    
        DenseMatrix_PL get_up(int dim, filtration_type a, filtration_type b);
        /**
         * Get the down persistent Laplacian (by reference).
         * @param dim dimension
         * @param a start filtration value
         * @param[out] L_down the matrix (by reference) of the down presistent Laplacian  
         */
        void get_down(int dim, filtration_type a, Eigen::SparseMatrix<storage> &L_down);
        /**
         * Get the down persistent Laplacian (by value). 
         * 
         * Warning: this does a potentially expensive copy. The primary reason this function exists is to provide reasonable python
         * binding access to the matrix itself. If you do not need this matrix directly in python (e.g. you want its eigenvalues)
         * it will be more efficient to compute by reference via get_down(int dim, filtration_type a, Eigen::SparseMatrix<storage> &L_down).
         * @param dim dimension
         * @param a start filtration value
         * @\return The matrix (by value) of the down presistent Laplacian  
         */
        DenseMatrix_PL get_down(int dim, filtration_type a);
        /**
         * Compute the nonzero eigenvalues of a persistent Laplacian using Schur restriction with the null space, either given by Persistent Homology representatives or by computing the null space.
         * @param dim dimension
         * @param a start filtration level
         * @param b end filtration level
         * @param (optional) PH_basis Basis for the null space of the Laplacian, possibly obtained through persistent homology
         * @param use_dummy_harmonic_basis Compute the null space of the Laplacian here, then do the same projected as if we had been given the null space
         * \return sorted vector of real, positive eigenvalues
         */
        std::vector<spectra_type> nonzero_spectra(int dim, filtration_type a, filtration_type b, SparseMatrix_PL PH_basis, bool use_dummy_harmonic_basis);

        /**********************************************/
        /* Driver functions to get PL and eigenvalues */
        /**********************************************/
        
        /**
         * Get the Persistent Laplacian's eigenvalues at a given dimension and filtration.
         * @param dim dimension
         * @param a start filtration level
         * @param b end filtration level
         * \return sorted vector of real, nonnegative eigenvalues 
         */
        std::vector<spectra_type> spectra(int dim, filtration_type a, filtration_type b);
        
        /**
         * Get all eigenvalues for all combinations of dimension and successive filtration values: a=filtrations[i] and b=filtrations[i+1]. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues.  
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra();
        
        /**
         * This function essentially just calls spectra(dim, a, b) in a loop.
         * @param spectra_quest_list vector of tuples (dim, a, b) to compute the eigenvalues of L_{dim}^{a,b}.
         * \return vector of tuples (dim, a, b, eigenvalues), where eigenvalues it istelf a vector of real, nonnegative eigenvalues.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra(std::vector<std::tuple<int,filtration_type,filtration_type>> spectra_request_list);
        
        /**
         * Get all eigenvalues for all combinations of dimension and filtration values. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues.  
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra_allpairs();
        
        /***************************************************************/
        /* Driver functions to get PL and eigenvalues and eigenvectors */
        /***************************************************************/
        

        /**
         * Get the Persistent Laplacian's eigenvalues and eigenvectors at a given dimension and filtration.
         * @param dim dimension
         * @param a start filtration level
         * @param b end filtration level
         * \return sorted pair of: vector of real, nonnegative eigenvalues and Eigen matrix where column i is the eigenvector for eigenvalue i
         */
        std::pair<std::vector<spectra_type>,DenseMatrix_PL> eigenpairs(int dim, filtration_type a, filtration_type b);
        
        /**
         * Get all eigenvalues for all combinations of dimension and successive filtration values: a=filtrations[i] and b=filtrations[i+1]. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues, eigenvectors) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues and "eigenvectors" is an Eigen::MatrixXf where column i is the eigenvector for eigenvalue i
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL>> eigenpairs();
        /**
         * This function essentially just calls eigenpairs(dim, a, b) in a loop.
         * @param spectra_quest_list vector of tuples (dim, a, b) to compute the eigenvalues of L_{dim}^{a,b}.
         * \return vector of tuples (dim, a, b, eigenvalues, eigenvalues), where eigenvalues it istelf a vector of real, nonnegative eigenvalues and "eigenvectors" is an Eigen::MatrixXf where column i is the eigenvector for eigenvalue i.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>, DenseMatrix_PL>> eigenpairs(std::vector<std::tuple<int,filtration_type,filtration_type>> spectra_request_list);

        /**
         * Utility function to get the Betti number and least nonzero eigenvalue form a vector of eigenvalues.
         * @param eigenvalues eigenvalues
         * \return pair of (Betti number, least nonzero eigenvalue)
         */
        std::pair<int, spectra_type> eigenvalues_summarize(std::vector<spectra_type> eigenvalues);
        /**********************************/
        /* Helpful input/output functions */
        /**********************************/
        

        /**
         * Compute and store a persistent Laplacian matrix in a file in matrix market format.
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value
         * @param filename file to store matrix (typically .mtx extension)
        */
        void store_L(int dim, filtration_type a, filtration_type b, std::string filename);
        /**
         * Print all boundaries and corresponding filtrations.
         */
        void print_boundaries();
        /**
         * Write spectra to files.
         * 
         * for each dimension dim of the complex, files write "{out_prefix}_spectra_{dim}.txt".
         * Each line of the file is a space-separated list of eigenvalues. Lines may be empty. Note the filtration values are not reported.
         * 
         * @param spectra tuples (dim, a, b, eigenvalues)
         * @param out_prefix Eigenvalues will be written to "{out_prefix}_spectra_{dim}.txt"
         *  
         */
        void store_spectra(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra, std::string out_prefix);

        /**
         * Write spectra summary to files "{out_prefix}_spectra_summary.txt"
         * 
         * Each line is a space-separated list of filtrations, bettti numbers, and least nonzero eigenvalues: (filtration a) (filtration b) (betti 0) ... (betti top_dim) (lambda 0) ... (lambda top_dim)
         * @param spectra tuples (dim, a, b, eigenvalues)
         * @param out_prefix Eigenvalues will be written to "{out_prefix}_spectra_summary.txt" 
         */
        void store_spectra_summary(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra, std::string out_prefix);


        /** 
         * Write the profile to a csv file. This is done automatically when spectra() is called, but it can also be done here.BDCSVDEigen
         * @param filename name of the file to write to
         */
        void time_to_csv(std::string filename){profile.to_csv(filename);}

        /**
         * Get tuples (dim, a, b) for all combinations of dimension and successive filtration values: a=filtrations[i], b=filtrations[i+1].
         * @param filtrations vector of filtration values
         * @param dims vector of dimensions
         * \return tuples (dim, a, b) for all combinations of dimension and successive filtration values: a=filtrations[i], b=filtrations[i+1].
         */
        std::vector<std::tuple<int, filtration_type, filtration_type>> filtration_list_to_spectra_request(std::vector<filtration_type> filtrations, std::vector<int> dims);
        
        /**
         * Get tuples (dim, a, b) for all combinations of dimension and filtrations.
         * @param filtrations vector of filtration values
         * @param dims vector of dimensions
         * \return tuples (dim, a, b) for all combinations of dimension and filtration values.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type>> filtration_list_to_spectra_request_allpairs(std::vector<filtration_type> filtrations, std::vector<int> dims);        
        
        
        /**
         * Get all unique filtration values in the complex.
         */
        std::vector<filtration_type> get_all_filtrations();
    protected:
    
        petls::FilteredBoundaryMatrix<storage> dummy_d0();
};


}
#endif