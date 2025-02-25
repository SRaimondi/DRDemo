//
// Created by simon on 13.05.17.
//

#ifndef DRDEMO_DIFF_OBJECT_HPP
#define DRDEMO_DIFF_OBJECT_HPP

#include "rad.hpp"

namespace drdemo {

    /**
     * Define differentiable object interface
     */
    class DiffObjectInterface {
    public:
        // Return a list of references to the differentiable variables in the object
        virtual void GetDiffVariables(std::vector<Float const *> &vars) const = 0;

        // Get number of differentiable variables of the object
        virtual size_t GetNumVars() const noexcept = 0;

        // Update values of differentiable variables in the object, starting index tells the object
        // from where to start to get the deltas
        virtual void UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) = 0;

        // Set value of the differentiable variables of the object
        virtual void SetDiffVariables(const std::vector<float> &vals, size_t starting_index) = 0;
    };

} // drdemo namespace

#endif //DRDEMO_DIFF_OBJECT_HPP
