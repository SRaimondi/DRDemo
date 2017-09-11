//
// Created by Simone Raimondi on 02.08.17.
//

#ifndef DRDEMO_SCALAR_FUNCTION_HPP
#define DRDEMO_SCALAR_FUNCTION_HPP

#include "rad.hpp"

namespace drdemo {

    // Define interface for Scalar functions to be used in the minimization algorithms
    class ScalarFunctionInterface {
    public:
        // Get input dimension
        virtual size_t InputDim() const = 0;

        // Evaluate the function at his current status
        virtual Float Evaluate(bool output) const = 0;

        // Compute gradient of the function, with respect to the given output variable
        virtual std::vector<float> ComputeGradient(const Float &out) const = 0;

        // Get value off the internal function status
        virtual std::vector<float> GetStatus() const = 0;

        // Update the internal state of the function given deltas (x_i = x_i + delta_i)
        virtual void UpdateStatus(const std::vector<float> &deltas) = 0;

        // Set internal status
        virtual void SetStatus(const std::vector<float> &new_status) = 0;

        // Request string description to function
        virtual std::string ToString() const = 0;
    };

    // General utility function

    // Gradient norm squared
    float GradientNorm2(const std::vector<float> &gradient);

    // Gradient norm
    float GradientNorm(const std::vector<float> &gradient);

    // Print gradient
    void PrintGradient(const std::vector<float> &gradient);

} // drdemo namespace

#endif //DRDEMO_SCALAR_FUNCTION_HPP
