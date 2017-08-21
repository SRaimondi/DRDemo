//
// Created by Simon on 10.08.2017.
//

#ifndef DRDEMO_GRADIENT_DESCENT_HPP
#define DRDEMO_GRADIENT_DESCENT_HPP

#include "scalar_function.hpp"

namespace drdemo {

    /**
     * Fixed step gradient descent routine to minimize a ScalarFunction
     */
    class GradientDescent {
    public:
        // Minimize given ScalarFunction given a step size, a maximum number of iterations,
        // a tolerance for the gradient norm.
        // The last argument tells if the routine should print information about the minimization process at every step
        static void Minimize(ScalarFunctionInterface &f,
                             float step,
                             size_t max_iters,
                             float tol,
                             bool verbose = false);
    };

} // drdemo namespace

#endif //DRDEMO_GRADIENT_DESCENT_HPP
