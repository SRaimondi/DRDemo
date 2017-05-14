//
// Created by simon on 08.05.17.
//

#ifndef DRDEMO_DERIVATIVE_HPP
#define DRDEMO_DERIVATIVE_HPP

#include <map>
#include "rad.hpp"

namespace drdemo {

    /**
     * The Derivatives class is used in conjunction with a variable to compute the derivative of a variable with respect
     * to another one
     */
    class Derivatives {
    private:
        // Associate variables with derivatives
        std::map<Float, std::vector<float> > var_derivatives_map;

    public:
        Derivatives() = default;

        // Clear derivatives
        void Clear();

        // Compute derivatives for a given variable
        void ComputeDerivatives(Float const &var);

        // Request derivative for a given outgoing variable with respect to a given input variable
        // Basically for df/dx, var_out is f and var_in x
        float Dwrt(Float const &f, Float const &x) const;
    };

} // drdemo namespace

#endif //DRDEMO_DERIVATIVE_HPP
