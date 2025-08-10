#ifndef PSL_H
#define PSL_H

#include "../core/Complex.hpp"
#include "sheaf_simplex_tree.hpp"

namespace petls{
    class PersistentSheafLaplacian : public Complex {
        public:
            // PersistentSheafLaplacian(sheaf_simplex_tree& sst) : petls::Complex() { // call default constructor of parent class
            //     std::vector<FilteredBoundaryMatrix<coefficient_type>> filtered_boundaries_temp = sst.apply_restriction_function();
            //     for (auto fbm : filtered_boundaries_temp){
            //         this->filtered_boundaries.push_back(fbm);
            //     }
            //     this->top_dim = this->filtered_boundaries.size();    
            //     for (int i = 0; i < (int) this->filtered_boundaries.size(); i++){
            //         this->filtered_boundaries[i] = this->filtered_boundaries[i].transpose();
            //     }
            // }
        private:


    };
}

#endif