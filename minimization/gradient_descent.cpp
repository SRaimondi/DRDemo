//
// Created by Simon on 10.08.2017.
//

#include "gradient_descent.hpp"

namespace drdemo {

    void GradientDescent::Minimize(ScalarFunctionInterface &f,
                                   float step,
                                   size_t max_iters,
                                   float grad_tol,
                                   bool verbose) {
        // Allocate vector passed to the function as updates and gradient
        std::vector<float> gradient(f.InputDim(), 0.f);
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
            Float result = f.Evaluate(true);
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
                // Print current iteration, function value and gradient norm
                std::cout << "Iteration: " << iters << std::endl;
                std::cout << "Function value: " << f_val << std::endl;
                std::cout << "Gradient norm: " << grad_norm << std::endl;
                // Print gradient
                PrintGradient(gradient);
                std::cout << std::endl;
            }

            // Increase number of iterations
            iters++;
        } while (grad_norm > grad_tol && iters < max_iters);
    }

    static std::vector<float>
    ComputeNewStatus(const std::vector<float> &c_s, const std::vector<float> &p_k, float alpha) {
        std::vector<float> new_status(c_s.size(), 0.f);

        for (size_t i = 0; i < new_status.size(); i++) {
            new_status[i] = c_s[i] + alpha * p_k[i];
        }

        return new_status;
    }

    float
    GradientDescentBT::ComputeStepSize(ScalarFunctionInterface &f,
                                       float f_x,
                                       const std::vector<float> &gradient,
                                       float c,
                                       float rho) {
        // Name as in the book, final step size is alpha
        float alpha = 1.f; // Positive, in the book p_k is the search direction and in our case is -grad(f)

        // Store current function status
        const std::vector<float> original_status = f.GetStatus();

        // Compute gradient norm, in the book is the dot product between the gradient and the search direction,
        // which in our case reduces to the negative squared norm of the gradient since we look in that direction
        const float grad_sqrd_norm = -GradientNorm2(gradient);

        // Negate gradient (p_k)
        std::vector<float> p_k(gradient.size(), 0.f);
        for (size_t i = 0; i < p_k.size(); i++) {
            p_k[i] = -gradient[i];
        }

        // Function status to test
        std::vector<float> test_status = ComputeNewStatus(original_status, p_k, alpha);
        // Set function to test status
        f.SetStatus(test_status);

        // Value of the function at current test status
        default_tape.Push();
        float f_eval = f.Evaluate(false).GetValue();
        default_tape.Pop();

        // Iterate to find good step size
        while (f_eval > f_x + c * alpha * grad_sqrd_norm) {
            // Push tape before evaluation
            default_tape.Push();

            // Compute status to test
            test_status = ComputeNewStatus(original_status, p_k, alpha);
            // Set new status to function
            f.SetStatus(test_status);
            // Evaluate function at new status
            f_eval = f.Evaluate(false).GetValue();

            // Pop tape after
            default_tape.Pop();

            // Update alpha
            alpha *= rho;
        }

        // Reset function to original status
        f.SetStatus(original_status);

        // Check that step is actually positive
        assert(alpha > 0.f);

        std::cout << "Step size: " << alpha << std::endl;

        return alpha;
    }

    void GradientDescentBT::Minimize(ScalarFunctionInterface &f,
                                     size_t max_iters,
                                     float grad_tol,
                                     float c,
                                     float rho,
                                     float step_tol,
                                     bool verbose) {
        // Allocate vector passed to the function as updates and gradient
        std::vector<float> gradient(f.InputDim(), 0.f);
        std::vector<float> deltas(f.InputDim(), 0.f);
        // Current value of the function
        float f_val;
        // Norm of the gradient
        float grad_norm;
        // Number of iterations
        size_t iters = 0;
        // Step size
        float alpha;

        // Minimize function
        do {
            // Push current status of the automatic differentiation tape
            default_tape.Push();

            // Evaluate function
            Float result = f.Evaluate(true);
            f_val = result.GetValue();

            // Compute gradient of the function
            gradient = f.ComputeGradient(result);
            grad_norm = GradientNorm(gradient);

            // Pop status from default tape
            default_tape.Pop();

            // Print information if verbose is true
            if (verbose) {
                // Print current iteration, function value and gradient norm
                std::cout << "Iteration: " << iters << std::endl;
                std::cout << "Function value: " << f_val << std::endl;
                // Print function information
                std::cout << f.ToString() << std::endl;
                std::cout << "Gradient norm: " << grad_norm << std::endl;
                // Print gradient
                // PrintGradient(gradient);         // FIXME Removed for debug
            }

            // Compute step size
            std::cout << "Computing step size..." << std::endl;
            alpha = ComputeStepSize(f, result.GetValue(), gradient, c, rho);
            std::cout << std::endl;

            // Compute updates for the function
            for (size_t i = 0; i < deltas.size(); ++i) {
                deltas[i] = -alpha * gradient[i];
            }

            // Update internal state of the function
            f.UpdateStatus(deltas);

            // Increase number of iterations
            iters++;
        } while (grad_norm > grad_tol && alpha > step_tol && iters < max_iters);
    }

} // drdemo namespace
