#ifndef PL_DF_H
#define PL_DF_H

#include "../core/Complex.hpp"
#include "../../include/Flagser_modified/src/flagser-laplacian.h"

namespace petls{
    class dFlag : public Complex{
        public:
            // dFlag(const char* filename, int max_dim); // .flag
            dFlag(const char* filename, int max_dim) : petls::Complex() { // call default constructor of parent class
                std::vector<SparseMatrixInt> boundaries;
                std::vector<std::vector<filtration_type>> filtrations;
                int min_dim = 0;
                basic_flagser(filename, min_dim, max_dim, boundaries, filtrations);
                // verbose = false;
                this->set_boundaries_filtrations(boundaries, filtrations);
            }
        private:


    };
}

#endif