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
        // a tolerance for the gradient norm and one for the step size, if it becomes to small, we stop
        // The last argument tells if the routine should print information about the minimization process at every step
        static void Minimize(ScalarFunctionInterface &f,
                             float step,
                             size_t max_iters,
                             float grad_tol,
                             bool verbose = false
        );
    };

    /**
     * Gradient descent using backtracking, from Numerical Optimization - Nocedal
     */
    class GradientDescentBT {
    private:
        /**
         * Compute step legth procedure
         */
        static float
        ComputeStepSize(ScalarFunctionInterface &f, float f_x, const std::vector<float> &gradient, float c,
                        float rho);

    public:
        // Minimize given scalar function using backtracking to compute step size
        // Other parameters are the maximum number of steps, the gradient tollerance to be accpeted as minuimum and the
        // c and rho value as state in Procedure 3.1 in the book, page 41-42
        // Next are a target function value and a tolerance on that, this can be used if we know the target value and we are happy with
        // a value near it in a certain margin.
        static void Minimize(ScalarFunctionInterface &f,
                             size_t max_iters,
                             float grad_tol,
                             float c,
                             float rho,
                             float step_tol,
                             bool verbose = false
        );
    };

} // drdemo namespace

#endif //DRDEMO_GRADIENT_DESCENT_HPP
