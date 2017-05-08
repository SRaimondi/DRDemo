//
// Created by simon on 08.05.17.
//

#ifndef DRDEMO_DERIVATIVE_HPP
#define DRDEMO_DERIVATIVE_HPP

#include <map>
#include "rad.hpp"

namespace drdemo {

    /**
     * The Derivatives class is used in conjugation with a variable to compute the derivatives
     */
    class Derivatives {
    private:
        // Associate variables with derivatives
        std::map<Float, std::vector<float> > var_derivatives_map;

    public:
        Derivatives() = default;

        // Compute derivatives for a given variable
        void ComputeDerivatives(Float const &var);

        // Request derivative for a given outgoing variable with respect to a given input variable
        // Basically for df/dx, var_out is f and var_in x
        float DfwrtDx(Float const &f, Float const &x) const;
    };

} // drdemo namespace

#endif //DRDEMO_DERIVATIVE_HPP
