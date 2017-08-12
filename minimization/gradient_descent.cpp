//
// Created by Simon on 10.08.2017.
//

#include "gradient_descent.hpp"

namespace drdemo {

    void GradientDescent::MinimizeFunction(ScalarFunctionInterface &f,
                                           float step,
                                           size_t max_iters,
                                           float tol,
                                           bool verbose) {
        // Allocate vector passed to the function as updates and gradient
        std::vector<float> gradient;
        std::vector<float> deltas(f.InputDim(), 0.f);
        // Current value of the function
        float f_val;
        // Norm of the gradient
        float grad_norm;
        // Number of iterations
        size_t iters = 0;

        // Minimize function
        do {
            // Push current status of the automatic differentiation tape
            default_tape.Push();

            // Evaluate function
            Float result = f.Evaluate();
            f_val = result.GetValue();

            // Compute gradient of the function
            gradient = f.ComputeGradient(result);
            grad_norm = GradientNorm(gradient);

            // Compute updates for the function
            for (size_t i = 0; i < deltas.size(); ++i) {
                deltas[i] = -step * gradient[i];
            }

            // Update internal state of the function
            f.UpdateStatus(deltas);

            // Pop status from default tape
            default_tape.Pop();

            // Print information if verbose is true
            if (verbose) {
                std::cout << "Iteration: " << iters << std::endl;
                std::cout << "Function value: " << f_val << std::endl;
                std::cout << "Gradient norm: " << grad_norm << std::endl;
                std::cout << std::endl;
            }

            // Increase number of iterations
            iters++;
        } while (grad_norm > tol && iters < max_iters);
    }

} // drdemo namespace
