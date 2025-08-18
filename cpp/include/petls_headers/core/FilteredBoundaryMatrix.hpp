#ifndef FBM_H
#define FBM_H

#include "../typedefs.hpp"

#include <vector>
#include <cassert>
#include <iostream>

namespace petls{
/**
 * Class for storing, manipulating, and retrieving the boundary matrix at different filtration levels 
 */
template<typename FBMcoeff = int> 
class FilteredBoundaryMatrix {
    /*********************************************************/
    /* Template: FBMcoeff: the type of entries in the matrix */
    /*********************************************************/
    private:
        /********************/
        /* Member variables */
        /********************/
        std::vector<filtration_type> domain_filtrations; ///< filtrations of n-simplices
        std::vector<filtration_type> range_filtrations; ///< filtrations of (n-1)-simplices
        Eigen::SparseMatrix<FBMcoeff, Eigen::ColMajor> matrix; ///< boundary matrix d_n
        int num_rows; ///< number of (n-1)-simplices
        int num_cols; ///<number of n-simplices

    public:
        /**
         * Constructor. Note: dimensions of filtrations and matrix must match.
         * @param _matrix boundary matrix d_n
         * @param _domain_filtrations filtrations of n-simplices
         * @param _range_filtrations filtrations of (n-1)-simplices
         */
        FilteredBoundaryMatrix(Eigen::SparseMatrix<FBMcoeff, Eigen::ColMajor> _matrix,
            std::vector<filtration_type> _domain_filtrations,
            std::vector<filtration_type> _range_filtrations){
                matrix = _matrix;
                domain_filtrations = _domain_filtrations;
                range_filtrations = _range_filtrations;                                          
                num_rows = _range_filtrations.size();
                num_cols = _domain_filtrations.size();
        }

        /***********************/
        /* Getters and setters */
        /***********************/

        /**
         * Get the filtration values of n-simplices
         */
        std::vector<filtration_type> get_domain_filtrations(){return domain_filtrations;}
        
        /**
         * Get the filtration values of (n-1)-simplices
         */
        std::vector<filtration_type> get_range_filtrations(){return range_filtrations;}

        /*********************/
        /* useful operations */
        /*********************/
        
        /**
         * Get the transpose of the boundary matrix, which swaps the domain and range filtrations 
         */
        FilteredBoundaryMatrix transpose(){
             return FilteredBoundaryMatrix(matrix.transpose(),range_filtrations,domain_filtrations);
        }

        /**
         * Get the largest index where the corresponding filtration value is not more than a.
         * If there are multiple simplices with the same filtration highest value, it will return the largest index.
         * @param use_domain_filtrations if true, iterate over the domain filtrations (n-simplices), if false over the range filtrations (n-1)-simplices
         * @param a filtration value   
         * \return The largest index where the corresponding filtration value is not more than a.
         */
        int index_of_filtration(bool use_domain_filtrations, filtration_type a){
            int index = 0;
            int filtration_size;
            if (use_domain_filtrations){
                filtration_size = domain_filtrations.size();
                while (index < filtration_size && domain_filtrations[index] <= a){ 
                    index++; 
                }
            } else{ //range_filtrations
                filtration_size = range_filtrations.size();
                while (index < filtration_size && range_filtrations[index] <= a){ 
                    index++;
                }
            }
            return index-1;
        }

        
        /**
         * Get the largest submatrix (by reference) that does not have corresponding filtration values greater than a.
         * @param a filtration value
         * @param[out] M the submatrix.
         */
        void submatrix_at_filtration(filtration_type a, Eigen::SparseMatrix<FBMcoeff,Eigen::ColMajor> &M){
            int col_index = index_of_filtration(true,a);
            int row_index = index_of_filtration(false,a);
            M = matrix.block(0,0,row_index+1,col_index+1);
        }

    
        

        /**********************************/
        /* Helpful input/output functions */
        /**********************************/

        /**
         * Print the filtration values of (n-1)-simplices.
         */
        void print_range_filtration(){
            std::cout << "\n[";
            for (unsigned long int i = 0; i < range_filtrations.size(); i++){
                std::cout << range_filtrations[i] << ",";
            }
            std::cout << "]\n";
        }

        /**
         * Print the filtration values of (n-1)-simplices.
         */
        void print_domain_filtration(){
            std::cout << "\n[";
            for (unsigned long int i = 0; i < domain_filtrations.size(); i++){
                std::cout << domain_filtrations[i] << ",";
            }
            std::cout << "]\n";
        }
        /**
         * Print the underlying matrix for d_n.
         */
        void print(){
            Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
            std::cout << Eigen::Matrix<FBMcoeff,Eigen::Dynamic, Eigen::Dynamic>(matrix).format(HeavyFmt) <<std::endl;
        }
};
}
#endif