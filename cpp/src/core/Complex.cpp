

#include <vector>

// #include "../../include/petls_headers/typedefs.hpp"
#include "../../include/petls_headers/core/Complex.hpp"
// #include "../../include/petls_headers/petls.hpp"
#include "../../include/petls_headers/core/up_algorithms.hpp"
#include "../../include/petls_headers/eigenvalues/eigs_algorithms.hpp"
// #include "../../include/petls_headers/core/FilteredBoundaryMatrix.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseQR>
#include <unsupported/Eigen/SparseExtra>
// #include <Eigen/Sparse>
#include <Eigen/Dense>

#include <chrono>
#include <numeric> //std::iota
#include <map>
#include <set>

#include <cassert>

namespace petls{

// template <typename storage = int>
    using storage = int;

        /********************/
        /* Member variables */
        /********************/


        /****************/
        /* Constructors */
        /****************/

        /**
         * Default constructor with no boundary maps or simplices.
         */
        Complex::Complex(){
            filtered_boundaries = std::vector<FilteredBoundaryMatrix<storage>>();
            top_dim = 0;
            verbose = false;
            use_flipped=false;
            eigs_algorithm_func = petls::SelfAdjointEigen;
            eigenpairs_algorithm_func = petls::SelfAdjointEigenpairsEigen;
            up_algorithm_func = petls::schur_algorithm<storage>;
            filtered_boundaries.push_back(dummy_d0());
        }

        /**
         * @param boundaries a vector of Eigen::SparseMatrix of type int. Boundaries must be sorted in order of dimension.
         * @param filtrations a vector of vector of filtrations. filtrations[dim] is a list of all filtrations of simplices in dimension dim. filtrations must be sorted in order of dimension, each filtrations[i] must be sorted in order of filtration
         *
         * Primary constructor. Important Assumptions:
         * 1) Boundary matrix has real coefficients stored as integers (but not mod 2!)
         * 2) Boundary matrix dimensions agree with filtrations sizes
         * 3) Length(filtrations) = Length(boundaries) + 1
         */
        Complex::Complex(std::vector<Eigen::SparseMatrix<storage>> boundaries,
                            std::vector<std::vector<filtration_type>> filtrations){
            filtered_boundaries = std::vector<FilteredBoundaryMatrix<storage>>();
            verbose = false;
            use_flipped=false;
            eigs_algorithm_func = petls::SelfAdjointEigen;
            eigenpairs_algorithm_func = petls::SelfAdjointEigenpairsEigen;
            up_algorithm_func = petls::schur_algorithm<storage>;
            this->set_boundaries_filtrations(boundaries, filtrations);
        }


        /***********************/
        /* Getters and setters */
        /***********************/

        void Complex::set_eigs_algorithm_func(std::function<spectra_vec(DenseMatrix_PL&)> _eigs_algorithm_func){
            eigs_algorithm_func = _eigs_algorithm_func;
        }

        void Complex::set_eigs_algorithm_func(std::string name){
            if (name == "selfadjoint"){
                eigs_algorithm_func = petls::SelfAdjointEigen;
            } else if (name == "eigensolver"){
                eigs_algorithm_func = petls::EigensolverEigen;
            } else if (name == "bdcsvd"){
                eigs_algorithm_func = petls::BDCSVDEigen;
            }
        }

        void Complex::set_up_func(std::function<void(FilteredBoundaryMatrix<int>*, filtration_type, filtration_type, DenseMatrix_PL&)> _up_algorithm_func){
            up_algorithm_func = _up_algorithm_func;
        }

        /**
         * Set the boundaries and filtrations of a complex, particularly if the default constructor was called.
         * @param boundaries a vector of Eigen::SparseMatrix of type int. Boundaries must be sorted in order of dimension.
         * @param filtrations a vector of vector of filtrations. filtrations[dim] is a list of all filtrations of simplices in dimension dim. filtrations must be sorted in order of dimension, each filtrations[i] must be sorted in order of filtration.
         * Important Assumptions:
         * 1) Boundary matrix has real coefficients stored as integers (but not mod 2!),
         * 2) Boundary matrix dimensions agree with filtrations sizes,
         * 3) Length(filtrations) = Length(boundaries) + 1.
         */
        void Complex::set_boundaries_filtrations(std::vector<Eigen::SparseMatrix<storage>> boundaries,
                            std::vector<std::vector<filtration_type>> filtrations){
            // Primary setter function called by constructors
            // Input:
            //      boundaries:     vector of Eigen::SparseMatrix of type int.
            //                      boundaries must be sorted in order of dimension
            //      filtrations:    vector of vector of filtrations.
            //                      filtrations[dim] is a list of all filtrations of
            //                      simplices in dimension dim.
            //                      filtrations must be sorted in order of dimension,
            //                      filtrations[i] must be sorted in order of filtration
            // Important Assumptions:
            //      1) Boundary matrix has real coefficients stored as integers (but not mod 2!)
            //      2) Boundary matrix dimensions agree with filtrations sizes
            //      3) Length(filtrations) = Length(boundaries) + 1
            eigs_algorithm_func = petls::SelfAdjointEigen;
            eigenpairs_algorithm_func = petls::SelfAdjointEigenpairsEigen;
            // up_algorithm_func = petls::schur_algorithm;
            filtered_boundaries.clear();
            top_dim = boundaries.size();
            filtered_boundaries.reserve(top_dim+1);
            
            filtered_boundaries.push_back(dummy_d0());

            // Stitch together boundary matrices with filtrations for domain and range
            // Check for consistent sizing
            for (unsigned long int dim = 1; dim <= (unsigned long int) top_dim; dim++){
                if(boundaries[dim-1].cols() != (long int) filtrations[dim].size()){
                    std::cout << "boundaries[" << dim -1 << "].cols()=" << boundaries[dim-1].cols() << " != filtrations[" << dim << "].size()=" << filtrations[dim].size() << std::endl;
                }
                if (boundaries[dim-1].rows() != (long int) filtrations[dim-1].size()){
                    std::cout << "boundaries[" << dim -1 << "].rows()=" << boundaries[dim-1].rows() << " != filtrations[" << dim-1 << "].size()=" << filtrations[dim-1].size() << std::endl;
                }
                filtered_boundaries.push_back(FilteredBoundaryMatrix<storage>(boundaries[dim-1],filtrations[dim],filtrations[dim-1]));
            }
        }
        
        
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
        void Complex::get_L(int dim, filtration_type a, filtration_type b, DenseMatrix_PL &L){ 
            // Get the Persistent Laplacian Matrix in dimension dim from filtration a to filtration b
            // Inputs:   integer dimension, start filtration level, end filtration level
            // Output (by reference to avoid a large copy): L of type Eigen::MatrixXf

            // L_0 = L_up
            if (dim == 0){
                this->profile.start_L_up();
                up_algorithm_func(&filtered_boundaries[dim+1],a,b,L);
                // up_Algorithm up_alg; // see up_algorithms.hpp
                // up_alg(&filtered_boundaries[dim+1],a,b, L);
                
                // Time monitoring
                this->profile.stop_L_up();
                profile.durations_sum_up_down.push_back(0);
                profile.durations_L_down.push_back(0);
                
                return;
            } else if (dim == top_dim){
                // L_{top_dim} = L_down
                this->profile.start_L_down();
                Eigen::SparseMatrix<storage> down(L.rows(),L.rows());
                get_down(dim, a, down);
                // Down produces a sparse matrix so we convert to dense
                L = DenseMatrix_PL(down.template cast<coefficient_type>());

                // Time monitoring
                this->profile.stop_L_down();
                profile.durations_sum_up_down.push_back(0);
                profile.durations_L_up.push_back(0);
            
                return;
            } else if (dim > top_dim){
                // L = 0
                profile.durations_sum_up_down.push_back(0);
                profile.durations_L_up.push_back(0);
                profile.durations_L_down.push_back(0);
                L.setZero(0,0);
                return;
            }
            // Else L = L_up + L_down
            
            // Set L = L_up then add L_down later
            this->profile.start_L_up();
            get_up(dim, a, b, L);
            this->profile.stop_L_up();

            this->profile.start_L_down();
            
            // Get L_down
            Eigen::SparseMatrix<storage> down(L.rows(), L.rows());
            get_down(dim,a, down);
            this->profile.stop_L_down();
            
            if (L.size() == 0){// L_up is empty, use L = L_down
                profile.durations_sum_up_down.push_back(0);
                L = DenseMatrix_PL(down.template cast<coefficient_type>());
                return;
            } else if (down.size() == 0){// L_down is empty, use L = L_up (already done)
                profile.durations_sum_up_down.push_back(0);
                return;
            }
            assert (L.cols() == down.cols() && L.rows() == down.rows() && "up and down Laplacians must have same dimensions");

            // L = L_up + L_down,
            // but L = L_up already, so just add L_down.
            // https://libeigen.gitlab.io/docs/group__TutorialSparse.html
            // this is faster than just returning up+down
            this->profile.start_sum_up_down();
            L+= down.template cast<coefficient_type>();
            this->profile.stop_sum_up_down();
            
            return;
            // L is passed by reference to avoid a copy
        }

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
        DenseMatrix_PL Complex::get_L(int dim, filtration_type a, filtration_type b){
            DenseMatrix_PL L;
            get_L(dim, a, b, L);
            return L;
        }

        /**
         * Get a matrix with the same nonzero eigenvalues as the top-dimensional Persistent Laplacian Matrix (by reference). 
         * @param a start filtration value
         * @param[out] L the matrix (by reference) that has the same eigenvalues as the top-dimensional Persistent Laplacian
         */
        void Complex::get_L_top_dim_flipped(filtration_type a, SparseMatrix_PL &L){
            Eigen::SparseMatrix<storage> B;
            filtered_boundaries[top_dim].submatrix_at_filtration(a, B);
            // The nonzero eigenvalues of (B * B^T) are the same as (B^T * B), and its possible BB^T is faster to compute the eigenvalues of.
            L = (B*B.transpose()).template cast<coefficient_type>();
        }

        /**
         * Get a matrix with the same nonzero eigenvalues as the top-dimensional Persistent Laplacian Matrix (by value). 
         * 
         * Warning: this does a potentially expensive copy. The primary reason this function exists is to provide reasonable python
         * binding access to the matrix itself. If you do not need this matrix directly in python (e.g. you want its eigenvalues)
         * it will be more efficient to compute by reference via get_L_top_dim_flipped(filtration_type a, SparseMatrix_PL &L).
         * @param a start filtration value
         * \return A matrix with the same nonzero eigenvalues as the top-dimensional persistent Laplacian (but not the persistent Laplacian itself).
         */
        SparseMatrix_PL Complex::get_L_top_dim_flipped(filtration_type a){
            SparseMatrix_PL L;
            get_L_top_dim_flipped(a, L);
            return L;
        }
        // void get_up_standard(int up_dim, filtration_type a, filtration_type b, DenseMatrix_PL &L_up);
        

        /**
         * Get the up persistent Laplacian (by reference). The algorithm used is determined by the template parameter. See up_algorithms.hpp.
         * 
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value
         * @param[out] L_up up persistent Laplacian (by reference)
         */
        void Complex::get_up(int dim, filtration_type a, filtration_type b, DenseMatrix_PL &L_up){
            up_algorithm_func(&filtered_boundaries[dim+1],a,b,L_up);
            // up_Algorithm up_alg; // see up_algorithms.hpp
            // up_alg(&filtered_boundaries[dim+1],a,b, L_up);
        }


        /**
         * Get the up persistent Laplacian (by value). The algorithm used is determined by the template parameter. See up_algorithms.hpp.
         * 
         * @param dim dimension
         * @param a start filtration value
         * @param b end filtration value
         * \return up persistent Laplacian
         */    
        DenseMatrix_PL Complex::get_up(int dim, filtration_type a, filtration_type b){
            DenseMatrix_PL L_up;
            get_up(dim, a, b, L_up);
            return L_up;
        }

        /**
         * Get the down persistent Laplacian (by reference).
         * @param dim dimension
         * @param a start filtration value
         * @param[out] L_down the matrix (by reference) of the down presistent Laplacian  
         */
        void Complex::get_down(int dim, filtration_type a, Eigen::SparseMatrix<storage> &L_down){ 
            Eigen::SparseMatrix<storage> B;
            filtered_boundaries[dim].submatrix_at_filtration(a, B);
            L_down = B.transpose()*B;

            // Every test of the self adjoint view rank update has been slower for down Laplacian than just B.transpose()*B.
            // L_down.selfadjointView<Eigen::Lower>().rankUpdate(B.transpose());
            return;
        }

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
        DenseMatrix_PL Complex::get_down(int dim, filtration_type a){
            Eigen::SparseMatrix<storage> L_down;
            get_down(dim, a, L_down);
            return DenseMatrix_PL(L_down.cast<coefficient_type>());
        }

        /**
         * Compute the nonzero eigenvalues of a persistent Laplacian using Schur restriction with the null space, either given by Persistent Homology representatives or by computing the null space.
         * @param dim dimension
         * @param a start filtration level
         * @param b end filtration level
         * @param (optional) PH_basis Basis for the null space of the Laplacian, possibly obtained through persistent homology
         * @param use_dummy_harmonic_basis Compute the null space of the Laplacian here, then do the same projected as if we had been given the null space
         * \return sorted vector of real, positive eigenvalues
         */
        std::vector<spectra_type> Complex::nonzero_spectra(int dim, filtration_type a, filtration_type b, SparseMatrix_PL PH_basis, bool use_dummy_harmonic_basis){
            spectra_vec eigenvalues;
            this->profile.start_L();
            // no 1-simplices, L_0 has no nonzero spectra
            if (dim == 0 && filtered_boundaries.size() == 1){
                int betti0 = filtered_boundaries[0].index_of_filtration(true, a)+1;
                profile.L_rows.push_back(betti0);
                return std::vector<spectra_type>(); // all spectra zero -> empty vector
            }

            // Get number of rows in Laplacian
            int L_rows;
            if (dim == 0){
                L_rows = filtered_boundaries[1].index_of_filtration(false,a)+1;
            } else {
                L_rows = filtered_boundaries[dim].index_of_filtration(true,a) + 1;
            }
            profile.L_rows.push_back(L_rows);
            
            // Note: PH Reduction solves a linear system, producing a dense matrix. We no longer have a sparse matrix in top dimension, so no need
            //       to have a separate case for the top dimension.

            // Get Laplacian matrix
            DenseMatrix_PL L(L_rows, L_rows);
            get_L(dim,a,b, L);

            this->profile.stop_L();
            // L = 0 return trivial
            if (L.size()==0){
                eigenvalues.setZero(0);
                profile.durations_eigs.push_back(0);
                return std::vector<spectra_type>();
            }

            int m = L.rows();
            int k;
            DenseMatrix_PL change_of_basis(m,m);
            // Get size of projection matrix
            if (use_dummy_harmonic_basis){
                // compute null space inefficiently for testing
                // use essentially this answer from StackOverflow: https://stackoverflow.com/a/53598471/3727807
                std::cout << "using inneficient harmonic basis" << std::endl; 
                Eigen::CompleteOrthogonalDecomposition<DenseMatrix_PL> cod;
                cod.compute(L);
                DenseMatrix_PL V = cod.matrixZ().transpose();
                DenseMatrix_PL Null_space = V.block(0, cod.rank(),V.rows(), V.cols() - cod.rank());
                DenseMatrix_PL P = cod.colsPermutation();
                DenseMatrix_PL PH_basis_dense = P * Null_space; // Unpermute the columns
                std::cout << "PH_basis_dense computed";
                // print_full_matrix_precise(PH_basis_dense);
                // here PH_basis is dense
                int n = PH_basis_dense.cols();
                k = m - n;
                change_of_basis.rightCols(PH_basis_dense.cols()) = PH_basis_dense;
            } else{
                // Use actual PH_basis
                // here PH_basis is sparse
                int n = PH_basis.cols();
                k = m - n;
                change_of_basis.rightCols(PH_basis.cols()) = PH_basis;
            }
                
            // A random matrix will be linearly independent with probability 1
            // A more intelligent basis could be used.
            // DenseMatrix_PL::Random
            DenseMatrix_PL nonharmonic_basis = DenseMatrix_PL::Random(m,k);
            change_of_basis.leftCols(nonharmonic_basis.cols()) = nonharmonic_basis;

            // Do the change of basis
            DenseMatrix_PL temp = change_of_basis.inverse()*L*change_of_basis; // TODO: optimize this
            DenseMatrix_PL Schur = temp.topLeftCorner(k,k);

            std::cout << "reduction by PH Schur matrix = " << Schur << std::endl;

            // Compute Eigenvalues of the smaller matrix
            // TODO: convert to parameterized eigenvalue algorithm
            this->profile.start_eigs();

            // petls::eigensolver es;
            eigenvalues = EigensolverEigen(Schur);
            // Eigen::EigenSolver<DenseMatrix_PL> es(Schur,Eigen::EigenvaluesOnly);
            // eigenvalues = es.eigenvalues(Schur);
            // Eigen::BDCSVD<DenseMatrix_PL> bdcsvd(Schur);
            // eigenvalues = bdcsvd.singularValues();
            this->profile.stop_eigs();

            // Copy to std::vector and sort
            std::vector<spectra_type> std_eigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
            std::sort(std_eigenvalues.begin(),std_eigenvalues.end()); // Standard EigenSolver can return in any order; we want sorted.
            
            return std_eigenvalues;
        }


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
        std::vector<spectra_type> Complex::spectra(int dim, filtration_type a, filtration_type b){
            this->profile.start_all();
            spectra_vec eigenvalues;
            profile.dims.push_back(dim);
            profile.filtration_a.push_back(a);
            profile.filtration_b.push_back(b);

            this->profile.start_L();
            // no 1-simplices, return vector of zeros
            if (dim == 0 && filtered_boundaries.size() == 1){
                int betti0 = filtered_boundaries[0].index_of_filtration(true, a)+1;
                // for (int i = 0; i < betti0; i++){
                //     eigenvalues.push_back(0.0);
                // }
                this->profile.stop_L();
                profile.L_rows.push_back(betti0);
                this->profile.stop_all();
                return std::vector<spectra_type>(betti0,0.0);
            }

            // Get number of rows in Laplacian
            int L_rows;
            if (dim == 0){
                L_rows = filtered_boundaries[1].index_of_filtration(false,a)+1;
            } else {
                L_rows = filtered_boundaries[dim].index_of_filtration(true,a) + 1;
            }
            profile.L_rows.push_back(L_rows);
            // Top dimension might use the "flipped" technique
            // This must be called from the spectra function so 
            // that the eigenvalues vector can be 0-padded correctly 
            if (dim == top_dim && use_flipped){
                int L_rows_flipped = filtered_boundaries[dim-1].index_of_filtration(true,a) + 1;
                int zero_pad_length = L_rows - L_rows_flipped;
               
                // Only use flipped version if it will be smaller
                if (L_rows_flipped < L_rows){
                    this->profile.start_L_down();
                    SparseMatrix_PL L(L_rows_flipped, L_rows_flipped);
                    get_L_top_dim_flipped(a, L);
                    this->profile.stop_L_down();
                    this->profile.stop_L();

                    // L=0 record trivial stats
                    if (L.size()==0){
                        eigenvalues.setZero(0);
                        profile.durations_eigs.push_back(0);
                        profile.durations_sum_up_down.push_back(0);
                        profile.durations_L_up.push_back(0);
                        this->profile.stop_all();
                        return std::vector<spectra_type>();
                    }

                    // Get eigenvalues of smaller L
                    // TODO: use parameterized eigenvalue solver
                    this->profile.start_eigs();
                    spectra_vec eigs = petls::SelfAdjointEigenSparse(L);
                    // Eigen::SelfAdjointEigenSolver<SparseMatrix_spectra_PL> es(L, Eigen::EigenvaluesOnly);
                    // spectra_vec eigs = es.eigenvalues();
                    std::vector<spectra_type> std_eigs_flipped(eigs.data(), eigs.data() + eigs.size());
                    if (zero_pad_length <= 0){
                        std::cout << "zero_pag_length <= 0 (should not happen)"<< std::endl;
                        this->profile.stop_eigs();
                    
                        profile.durations_sum_up_down.push_back(0);
                        profile.durations_L_up.push_back(0);
                        this->profile.stop_all();
                        return std_eigs_flipped;
                    }

                    // Pad the eigenvalues with 0s
                    std::vector<spectra_type> zero_pad(zero_pad_length, 0.0); 
                    std::move(std_eigs_flipped.begin(), std_eigs_flipped.end(), std::back_inserter(zero_pad));
                    
                    // Record timing
                    this->profile.stop_eigs();
                    profile.durations_sum_up_down.push_back(0);
                    profile.durations_L_up.push_back(0);
                    this->profile.stop_all();
                    return zero_pad;                    
                } 
            }
            
            // Get Laplacian matrix
            DenseMatrix_PL L(L_rows, L_rows);
            get_L(dim,a,b, L);
            this->profile.stop_L();
            
            // L = 0, recover trivial stats
            if (L.size()==0){
                eigenvalues.setZero(0);
                profile.durations_eigs.push_back(0);
                this->profile.stop_all();
                return std::vector<spectra_type>();
            }

            // Compute Eigenvalues
            this->profile.start_eigs();
            eigenvalues = eigs_algorithm_func(L);
            this->profile.stop_eigs();

            std::vector<spectra_type> std_eigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
            this->profile.stop_all();
            return std_eigenvalues;
        }

        /**
         * Get all eigenvalues for all combinations of dimension and successive filtration values: a=filtrations[i] and b=filtrations[i+1]. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues.  
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> Complex::spectra(){
            // Output:
            //      A vector of tuples (dim, a, b, eigenvalues), where the eigenvalues is a sorted vector of real numbers     
            

            // Create a vector listing all dimensions 0, 1, ..., top_dim
            std::vector<int> dims(top_dim+1);
            std::iota (std::begin(dims), std::end(dims), 0);
            
            // Get all filtration values that occur in the complex
            std::vector<filtration_type> all_filtrations = get_all_filtrations();

            // Convert the vector of dimensions and filtration values into triples
            // (dim, a, b) to get L_{dim}^{a,b}
            std::vector<std::tuple<int, filtration_type, filtration_type>> requests = filtration_list_to_spectra_request(all_filtrations,dims);
            
            // Compute the PL matrices and eigenvalues
            std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> requested_spectra = spectra(requests);

            // Record results for the profiler, not necessary for the user outside of benchmarking
            for (int i = 0; i < (int) requested_spectra.size(); i++){
                std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>> current = requested_spectra[i];
                // profile.dims.push_back(std::get<0>(current));
                // profile.filtration_a.push_back(std::get<1>(current));
                // profile.filtration_b.push_back(std::get<2>(current));
                std::pair<int,spectra_type> summary = eigenvalues_summarize(std::get<3>(current));
                profile.bettis.push_back(summary.first);
                profile.lambdas.push_back(summary.second);
            }
            // Output profiler to csv
            if (filtered_boundaries.size() > 1)// dont report profile when no 1-simplices
                profile.to_csv("./profile.csv");
            return requested_spectra;
        }
        
        /**
         * This function essentially just calls spectra(dim, a, b) in a loop.
         * @param spectra_quest_list vector of tuples (dim, a, b) to compute the eigenvalues of L_{dim}^{a,b}.
         * \return vector of tuples (dim, a, b, eigenvalues), where eigenvalues it istelf a vector of real, nonnegative eigenvalues.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> Complex::spectra(std::vector<std::tuple<int,filtration_type,filtration_type>> spectra_request_list){
            // Declare variables and vectors
            int dim;
            filtration_type a;
            filtration_type b;
            std::tuple<int,filtration_type,filtration_type>  spectra_request;
            std::vector<spectra_vec> spectra_list;
            std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> requested_spectra;
            requested_spectra.reserve(spectra_request_list.size());
            
            // Loop through list of tuples (dim, a, b)
            for (unsigned long int i = 0; i < spectra_request_list.size(); i++){
                // Unpack the tuples
                spectra_request = spectra_request_list[i];
                dim = std::get<0>(spectra_request);
                a = std::get<1>(spectra_request);
                b = std::get<2>(spectra_request);    

                // Get eigenvalues and re-pack as (dim, a, b, eigenvalues)
                // this->profile.start_all();
                requested_spectra.push_back(std::make_tuple(dim, a, b, spectra(dim, a, b)));
                // this->profile.stop_all();
                if (verbose)
                    std::cout << "duration spectra for dim=" << dim << ", a= " << a << ", b=" << b << ": " << this->profile.all.duration << std::endl;
            }
            return requested_spectra;
        }

        /**
         * Get all eigenvalues for all combinations of dimension and filtration values. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues.  
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> Complex::spectra_allpairs(){
            // Output:
            //      A vector of tuples (dim, a, b, eigenvalues), where the eigenvalues is a sorted vector of real numbers     
            

            // Create a vector listing all dimensions 0, 1, ..., top_dim
            std::vector<int> dims(top_dim+1);
            std::iota (std::begin(dims), std::end(dims), 0);
            
            // Get all filtration values that occur in the complex
            std::vector<filtration_type> all_filtrations = get_all_filtrations();

            // Convert the vector of dimensions and filtration values into triples
            // (dim, a, b) to get L_{dim}^{a,b}
            std::vector<std::tuple<int, filtration_type, filtration_type>> requests = filtration_list_to_spectra_request_allpairs(all_filtrations,dims);
            
            // Compute the PL matrices and eigenvalues

            std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> requested_spectra = spectra(requests);

            // Record results for the profiler, not necessary for the user outside of benchmarking
            for (int i = 0; i < (int) requested_spectra.size(); i++){
                std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>> current = requested_spectra[i];
                // profile.dims.push_back(std::get<0>(current));
                // profile.filtration_a.push_back(std::get<1>(current));
                // profile.filtration_b.push_back(std::get<2>(current));
                std::pair<int,spectra_type> summary = eigenvalues_summarize(std::get<3>(current));
                profile.bettis.push_back(summary.first);
                profile.lambdas.push_back(summary.second);
            }
            // Output profiler to csv
            if (filtered_boundaries.size() > 1)// dont report profile when no 1-simplices
                profile.to_csv("./profile.csv");
            return requested_spectra;
        }
        


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
        std::pair<std::vector<spectra_type>,DenseMatrix_PL> Complex::eigenpairs(int dim, filtration_type a, filtration_type b){
            spectra_vec eigenvalues;
            DenseMatrix_PL eigenvectors;
            this->profile.start_L();
            // no 1-simplices, return all eigenvalues zero
            // each eigenvector is a unit vector in only one component, so the matrix of all eigenvectors is just the identity matrux
            if (dim == 0 && filtered_boundaries.size() == 1){
                int betti0 = filtered_boundaries[0].index_of_filtration(true, a)+1;
                this->profile.stop_L();
                profile.L_rows.push_back(betti0);
                std::vector<spectra_type> std_eigs = std::vector<spectra_type>(betti0,0.0);
                eigenvectors.setIdentity(betti0,betti0);
                return std::pair<std::vector<spectra_type>,DenseMatrix_PL>(std_eigs, eigenvectors);
            }

            // Get number of rows in Laplacian
            int L_rows;
            if (dim == 0){
                L_rows = filtered_boundaries[1].index_of_filtration(false,a)+1;
            } else {
                L_rows = filtered_boundaries[dim].index_of_filtration(true,a) + 1;
            }
            profile.L_rows.push_back(L_rows);
            
            // Get L
            DenseMatrix_PL L(L_rows, L_rows);
            get_L(dim,a,b, L);
            this->profile.stop_L();
            
            // L = 0, recover trivial stats
            if (L.size()==0){
                eigenvalues.setZero(0);
                profile.durations_eigs.push_back(0);
                return std::pair<std::vector<spectra_type>,DenseMatrix_PL>(std::vector<spectra_type>(),DenseMatrix_PL());
            }

            // Compute Eigenvalues
            this->profile.start_eigs();
            std::pair<spectra_vec,DenseMatrix_spectra_PL> eigenpairs = eigenpairs_algorithm_func(L);
            // eigs_Algorithm es;
            // std::pair<spectra_vec,DenseMatrix_spectra_PL> eigenpairs = es.eigenpairs(L);
            eigenvalues = eigenpairs.first;
            this->profile.stop_eigs();
            std::vector<spectra_type> std_eigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
            return std::pair<std::vector<spectra_type>,DenseMatrix_PL>(std_eigenvalues,eigenpairs.second);
        }

        /**
         * Get all eigenvalues for all combinations of dimension and successive filtration values: a=filtrations[i] and b=filtrations[i+1]. Note: the caller does not know what spectra to expect from this.       
         * \return vector of tuples (dim, a, b, eigenvalues, eigenvectors) where "eigenvalues" is a sorted vector of real, nonnegative eigenvalues and "eigenvectors" is an Eigen::MatrixXf where column i is the eigenvector for eigenvalue i
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL>> Complex::eigenpairs(){

            // Create a vector listing all dimensions 0, 1, ..., top_dim
            std::vector<int> dims(top_dim+1);
            std::iota (std::begin(dims), std::end(dims), 0);
            
            // Get all filtration values that occur in the complex
            std::vector<filtration_type> all_filtrations = get_all_filtrations();

            // Convert the vector of dimensions and filtration values into triples
            // (dim, a, b) to get L_{dim}^{a,b}
            std::vector<std::tuple<int, filtration_type, filtration_type>> requests = filtration_list_to_spectra_request(all_filtrations,dims);
            
            // Compute the PL matrices and eigenvalues
            std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL>> requested_spectra = eigenpairs(requests);
            
            // Record results for the profiler, not necessary for the user outside of benchmarking
            for (int i = 0; i < (int) requested_spectra.size(); i++){
                std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>,DenseMatrix_PL> current = requested_spectra[i];
                // profile.dims.push_back(std::get<0>(current));
                // profile.filtration_a.push_back(std::get<1>(current));
                // profile.filtration_b.push_back(std::get<2>(current));
                std::pair<int,spectra_type> summary = eigenvalues_summarize(std::get<3>(current));
                profile.bettis.push_back(summary.first);
                profile.lambdas.push_back(summary.second);
            }
            // Output profiler to csv
            if (filtered_boundaries.size() > 1)// dont report profile when no 1-simplices
                profile.to_csv("./profile.csv");
            return requested_spectra;
        }
        
        /**
         * This function essentially just calls eigenpairs(dim, a, b) in a loop.
         * @param spectra_quest_list vector of tuples (dim, a, b) to compute the eigenvalues of L_{dim}^{a,b}.
         * \return vector of tuples (dim, a, b, eigenvalues, eigenvalues), where eigenvalues it istelf a vector of real, nonnegative eigenvalues and "eigenvectors" is an Eigen::MatrixXf where column i is the eigenvector for eigenvalue i.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>, DenseMatrix_PL>> Complex::eigenpairs(std::vector<std::tuple<int,filtration_type,filtration_type>> spectra_request_list){
            // Declare variables and vectors
            int dim;
            filtration_type a;
            filtration_type b;
            std::tuple<int,filtration_type,filtration_type>  spectra_request;
            std::vector<spectra_vec> spectra_list;
            std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>, DenseMatrix_PL>> requested_spectra;
            requested_spectra.reserve(spectra_request_list.size());
            
            // Loop through list of tuples (dim, a, b)
            for (unsigned long int i = 0; i < spectra_request_list.size(); i++){
                // Unpack the tuples
                spectra_request = spectra_request_list[i];
                dim = std::get<0>(spectra_request);
                a = std::get<1>(spectra_request);
                b = std::get<2>(spectra_request);    

                // Get eigenvalues and eigenvectors and re-pack as (dim, a, b, eigenvalues, eigenvectors)
                this->profile.start_all();
                std::pair<std::vector<spectra_type>,DenseMatrix_PL> eigenpairs = this->eigenpairs(dim, a, b);
                std::vector<spectra_type> eigenvalues = eigenpairs.first;
                DenseMatrix_PL eigenvectors = eigenpairs.second;
                requested_spectra.push_back(std::make_tuple(dim, a, b, eigenvalues, eigenvectors));
                this->profile.stop_all();
                if (verbose)
                    std::cout << "duration spectra for dim=" << dim << ", a= " << a << ", b=" << b << ": " << this->profile.all.duration << std::endl;
            }
            return requested_spectra;
        }


        /**
         * Utility function to get the Betti number and least nonzero eigenvalue form a vector of eigenvalues.
         * @param eigenvalues eigenvalues
         * \return pair of (Betti number, least nonzero eigenvalue)
         */
        std::pair<int, spectra_type> Complex::eigenvalues_summarize(std::vector<spectra_type> eigenvalues){
            // Input: vector of eigenvalues
            // Output: betti number and least nonzero eigenvalue (tolerance 1e-4)
            int current_betti = 0;
            spectra_type tol = 1e-4;
            for (int k = 0; k < (int) eigenvalues.size(); k++){
                if (eigenvalues[k] > tol){ // reached a nonzero eigenvalue
                    return std::make_pair(current_betti, eigenvalues[k]);
                }
                current_betti++;
            }
            // if we reach the end of the loop, then there were no nonzero eigenvalues
            if (current_betti == 0){
                if ((int) eigenvalues.size() > 0){
                    return std::make_pair(0,eigenvalues[0]);
                }
                return std::make_pair(0,0.0);
            }
            // else never encountered a nonzero eigenvalue
            return std::make_pair(current_betti,0.0);
        }

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
        void Complex::store_L(int dim, filtration_type a, filtration_type b, std::string filename){
            int L_rows = filtered_boundaries[dim].index_of_filtration(true,a) + 1;
            DenseMatrix_PL L(L_rows, L_rows);
            
            get_L(dim,a,b, L);
            SparseMatrix_PL L_sparse = L.sparseView();
            bool success = Eigen::saveMarket(L_sparse, filename); // conversion to sparse is probably extremely expensive, but not all versions of eigen have saveMarketDense
            if (success)
                std::cout << "saved matrix to file " << filename << std::endl;
            else
                std::cout << "failed saving matrix to file " << filename << std::endl;
            return;
        }
    
        /**
         * Print all boundaries and corresponding filtrations.
         */
        void Complex::print_boundaries(){

            for (int i = 0; i <= top_dim; i++){
    
                filtered_boundaries[i].print();
    
                filtered_boundaries[i].print_range_filtration();
    
                filtered_boundaries[i].print_domain_filtration();
    
            }
        }

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
        void Complex::store_spectra(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra, std::string out_prefix){
            // Open the files
            std::vector<std::ofstream> out_streams(top_dim+1);
            for (int i = 0; i <= top_dim; i++){
                out_streams[i] = std::ofstream("./" + out_prefix + "_spectra_" + std::to_string(i) + ".txt");
            }

            // Write the eigenvalues (not assuming any ordering on the dimensions)
            for (int i = 0; i < (int) spectra.size(); i++){
                std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>> spectrum_record = spectra[i];
                int dim = std::get<0>(spectrum_record);
                std::vector<spectra_type> eigs = std::get<3>(spectrum_record);
                for (int j = 0; j < (int) eigs.size(); j++){
                    out_streams[dim] << eigs[j] << " "; 
                }
                out_streams[dim] << std::endl;
            }
            // Close the files   
            for (int i = 0; i <= top_dim; i++){
                out_streams[i].close();
            }
        }


        /**
         * Write spectra summary to files "{out_prefix}_spectra_summary.txt"
         * 
         * Each line is a space-separated list of filtrations, bettti numbers, and least nonzero eigenvalues: (filtration a) (filtration b) (betti 0) ... (betti top_dim) (lambda 0) ... (lambda top_dim)
         * @param spectra tuples (dim, a, b, eigenvalues)
         * @param out_prefix Eigenvalues will be written to "{out_prefix}_spectra_summary.txt" 
         */
        void Complex::store_spectra_summary(std::vector<std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>>> spectra, std::string out_prefix){
            std::set<std::pair<filtration_type, filtration_type>> filtration_values_set;
            
            // get all unique (a,b) filtration value pairs
            for (int i = 0 ; i < (int) spectra.size(); i++){
                filtration_type a = std::get<1>(spectra[i]);
                filtration_type b = std::get<2>(spectra[i]);
                filtration_values_set.insert(
                    std::make_pair(a,b));
            }
            int num_filtrations = filtration_values_set.size();
            std::vector<std::pair<filtration_type,filtration_type>> filtration_values_vec(filtration_values_set.begin(), filtration_values_set.end());

            // create a mapping for their index
            // unordered_map may be much faster, but cannot hash a std::pair
            std::map<std::pair<filtration_type,filtration_type>,int> filtration_index_map;
            std::vector<std::vector<spectra_type>> output_lines(num_filtrations);
            int items_per_line = 2 + 2*(top_dim+1);
            for (int i = 0; i < num_filtrations; i++){
                filtration_index_map[filtration_values_vec[i]] = i;
                
                std::vector<spectra_type> temp_line(items_per_line);
                temp_line[0] = (spectra_type) filtration_values_vec[i].first;
                temp_line[1] = (spectra_type) filtration_values_vec[i].second;
                output_lines[i] = temp_line;
            }

            for (int i = 0; i < (int) spectra.size(); i++){
                // Unpack each tuple
                std::tuple<int, filtration_type, filtration_type, std::vector<spectra_type>> spectrum_record = spectra[i];
                int dim = std::get<0>(spectrum_record);
                filtration_type a = std::get<1>(spectrum_record);
                filtration_type b = std::get<2>(spectrum_record);
                
                // compute eigenvalue summary
                std::pair<filtration_type, filtration_type> filtration_pair = std::make_pair(a,b);
                std::vector<spectra_type> eigs = std::get<3>(spectrum_record);
                std::pair<int, spectra_type> temp_pair = eigenvalues_summarize(eigs);
                spectra_type betti = (spectra_type) temp_pair.first;
                spectra_type lambda = (spectra_type) temp_pair.second;
                
                int line_index = filtration_index_map[filtration_pair];
                // a, b, betti_0, betti_1, betti_2, lambda_0, lambda_1, lambda_2
                // 0, 1,       2,       3,       4,        5,        6,        7
                output_lines[line_index][2+dim] = betti;
                output_lines[line_index][3+top_dim+dim] = lambda;
            }
            
            // open file and write headers
            std::ofstream outstream("./" + out_prefix + "_spectra_summary.txt");
            outstream << "a\tb";
            for (int i = 0; i <= top_dim; i++){
                outstream <<"\tbetti_" << i;
            }
            for (int i = 0; i <= top_dim; i++){
                outstream <<"\tlambda_" << i;
            }
            outstream << std::endl;

            // write the data to the file
            for (int i = 0; i < (int) num_filtrations; i++){
                outstream << output_lines[i][0];
                for (int j = 1; j < items_per_line; j++){
                    outstream << "\t" << output_lines[i][j];
                } 
                outstream << std::endl;
            }
            outstream.close();
        }       
        

        /**
         * Get tuples (dim, a, b) for all combinations of dimension and successive filtration values: a=filtrations[i], b=filtrations[i+1].
         * @param filtrations vector of filtration values
         * @param dims vector of dimensions
         * \return tuples (dim, a, b) for all combinations of dimension and successive filtration values: a=filtrations[i], b=filtrations[i+1].
         */
        std::vector<std::tuple<int, filtration_type, filtration_type>> Complex::filtration_list_to_spectra_request(std::vector<filtration_type> filtrations, std::vector<int> dims){            
            // Declare variables
            filtration_type a;
            filtration_type b;
            int dim;
            std::vector<std::tuple<int, filtration_type, filtration_type>> requests;
            std::tuple<int, filtration_type, filtration_type> spectra_request;
            
            // create the pairs in order (dim=0, a=0, b=1), (dim=1, a=0, b=1), ..., (dim=top_dim, a=0, b=1), (dim=0, a=1, b=2), ...
            for(unsigned long int filtration_index = 0; filtration_index < filtrations.size()-1; filtration_index++){
                a = filtrations[filtration_index];
                b = filtrations[filtration_index+1];
                for (unsigned long int dim_index = 0; dim_index < dims.size(); dim_index++){
                    dim = dims[dim_index];
                    spectra_request = std::make_tuple(dim, a, b);
                    requests.push_back(spectra_request);
                }
            }
            // add the [a,infinity) case as (dim=0, a, a), (dim=1, a, a), ... 
            a = filtrations[filtrations.size()-1];
            b = a;
            for (unsigned long int dim_index = 0; dim_index < dims.size(); dim_index++){
                dim = dims[dim_index];
                spectra_request = std::make_tuple(dim, a, b);
                requests.push_back(spectra_request);
            }
            return requests;
        }

        /**
         * Get tuples (dim, a, b) for all combinations of dimension and filtrations.
         * @param filtrations vector of filtration values
         * @param dims vector of dimensions
         * \return tuples (dim, a, b) for all combinations of dimension and filtration values.
         */
        std::vector<std::tuple<int, filtration_type, filtration_type>> Complex::filtration_list_to_spectra_request_allpairs(std::vector<filtration_type> filtrations, std::vector<int> dims){            
            // Declare variables
            filtration_type a;
            filtration_type b;
            int dim;
            std::vector<std::tuple<int, filtration_type, filtration_type>> requests;
            std::tuple<int, filtration_type, filtration_type> spectra_request;
            
            // create the pairs in order (dim=0, a=0, b=1), (dim=1, a=0, b=1), ..., (dim=top_dim, a=0, b=1), (dim=0, a=0, b=2), ...
            for(unsigned long int filtration_index_a = 0; filtration_index_a < filtrations.size()-1; filtration_index_a++){
                a = filtrations[filtration_index_a];
                for(unsigned long int filtration_index_b = filtration_index_a; filtration_index_b < filtrations.size(); filtration_index_b++){
                    b = filtrations[filtration_index_b];
                    for (unsigned long int dim_index = 0; dim_index < dims.size(); dim_index++){
                        dim = dims[dim_index];
                        spectra_request = std::make_tuple(dim, a, b);
                        requests.push_back(spectra_request);
                    }
                }
            }
            // add the [a,infinity) case as (dim=0, a, a), (dim=1, a, a), ... 
            a = filtrations[filtrations.size()-1];
            b = a;
            for (unsigned long int dim_index = 0; dim_index < dims.size(); dim_index++){
                dim = dims[dim_index];
                spectra_request = std::make_tuple(dim, a, b);
                requests.push_back(spectra_request);
            }
            return requests;
        }

        /**
         * Get all unique filtration values in the complex.
         */
        std::vector<filtration_type> Complex::get_all_filtrations(){
            std::set<filtration_type> all_filtrations;

            //add zero-th dimensional
            std::vector<filtration_type> c0 = filtered_boundaries[0].get_domain_filtrations();
            for (unsigned long int i = 0; i < c0.size(); i++){
                all_filtrations.insert(c0[i]);
            }
            //add other dimensions, which correspond to domain filtrations
            std::vector<filtration_type> temp;
            for (unsigned long int dim = 1; dim <= (unsigned long int) top_dim; dim++){
                temp = filtered_boundaries[dim].get_domain_filtrations();
                for (unsigned long int i = 0; i < temp.size(); i++){
                    all_filtrations.insert(temp[i]);
                } 
            };
            std::vector<filtration_type> all_filtrations_vector(all_filtrations.begin(), all_filtrations.end());
            return all_filtrations_vector;
        }


        petls::FilteredBoundaryMatrix<storage> Complex::dummy_d0(){
            // Construct a placeholder matrix for d0 so that filtered_boundaries[dim] actually gives d_{dim}
            Eigen::SparseMatrix<storage> dummy_d0_matrix;
            std::vector<filtration_type> domain_filtrations = {0.0};
            std::vector<filtration_type> range_filtrations = {0.0};
            petls::FilteredBoundaryMatrix<storage> dummy_d0_fbm(dummy_d0_matrix,domain_filtrations,range_filtrations);
            return dummy_d0_fbm;
        }    

}