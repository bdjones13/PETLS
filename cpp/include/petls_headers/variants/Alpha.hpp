#ifndef PL_A_H
#define PL_A_H
#include "../core/Complex.hpp"


namespace petls{
    class Alpha : public Complex {
        public:
            Alpha(const char* filename, int max_dim); // .OFF file 
            Alpha(std::vector<std::tuple<double,double,double>> points, int max_dim); // coordinantes
                                                                                                        // note that public API not using Gudhi type Point_3
                                                                                                        // so caller does not need knowledge of Gudhi
                                                                                                        // TODO: should pass points by reference, but pybind11 has problems

    };
}

#endif 