//
// Created by simon on 08.05.17.
//

#include <iostream>
#include "derivative.hpp"

namespace drdemo {

    void Derivatives::Clear() {
        var_derivatives_map.clear();
    }

    void Derivatives::ComputeDerivatives(Float const &var) {
        size_t const size = default_tape.Size();
        // Check if element is already in the map
        if (var_derivatives_map.find(var) != var_derivatives_map.end()) {
            std::cout << "Variable derivatives already computed, exiting..." << std::endl;
            return;
        }
        // Create new entry
        var_derivatives_map[var] = std::vector<float>(size, 0.f);
        // Get reference to make code more clear
        std::vector<float> &derivs = var_derivatives_map[var];

        // Seed at index of the input variable
        derivs[var.NodeIndex()] = 1.f;

        // Traverse the tape in reverse
        for (size_t i = 0; i < size; i++) {
            size_t const rev_index = size - 1 - i;
            TapeNode const &node = default_tape.At(rev_index);

            // Add children contribution
            if (node.parent_i[0] == NOT_REGISTERED || node.parent_i[1] == NOT_REGISTERED) {
                std::cerr << "Warning! Trying to compute derivative that involves unregistered nodes!" << std::endl;
                exit(EXIT_FAILURE);
            }
            derivs[node.parent_i[0]] += node.weights[0] * derivs[rev_index];
            derivs[node.parent_i[1]] += node.weights[1] * derivs[rev_index];
        }
    }

    float Derivatives::Dwrt(Float const &f, Float const &x) const {
        // Check if we computed the derivatives for the given out variable
        auto it = var_derivatives_map.find(f);
        if (it == var_derivatives_map.end()) {
            std::cout << "Could not find given output variable!" << std::endl;
            return 0.f;
        }
        return it->second[x.NodeIndex()];
    }

} // drdemo namespace