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

//            if (std::isnan(derivs[rev_index])) {
//                int a = 1;
//            }
//
//            // DEBUG ASSERTIONS
////            assert(!std::isinf(node.weights[0]));
////            assert(!std::isinf(node.weights[1]));
////            assert(!std::isnan(derivs[node.parent_i[0]]));
////            assert(!std::isnan(derivs[node.parent_i[1]]));
//
//            // Need to check for the case of doing Inf * 0.f -> Nan
//            float new_value;
//            if ((std::isinf(node.weights[0]) && derivs[rev_index] == 0.f) ||
//                (node.weights[0] == 0.f && std::isinf(derivs[rev_index]))) {
//                // Add derivative to 0.f
//                derivs[node.parent_i[0]] += 0.f; // TODO Fix this for clarity
//            } else {
//                float temp = node.weights[0] * derivs[rev_index];
//                if (std::isnan(temp)) { // FIXME
//                    int bla = 1;
//                }
//                new_value = derivs[node.parent_i[0]] + temp;
//                // derivs[node.parent_i[0]] += temp;
//            }
//            if(std::isnan(new_value)) {
//                int bla = 1;
//            }
//            derivs[node.parent_i[0]] = new_value;
//
//            if ((std::isinf(node.weights[1]) && derivs[rev_index] == 0.f) ||
//                (node.weights[1] == 0.f && std::isinf(derivs[rev_index]))) {
//                // Add derivative to 0.f
//                derivs[node.parent_i[1]] += 0.f; // TODO Fix this for clarity
//            } else {
//                float temp = node.weights[1] * derivs[rev_index];
//                if (std::isnan(temp)) { // FIXME
//                    int bla = 1;
//                }
//                derivs[node.parent_i[1]] += temp;
//            }
//            if(std::isnan(derivs[node.parent_i[1]])) {
//                int bla = 1;
//            }

            // Add children contribution
            derivs[node.parent_i[0]] += node.weights[0] * derivs[rev_index];
            derivs[node.parent_i[1]] += node.weights[1] * derivs[rev_index];
        }

        // DEBUG ASSERTIONS
//        for (size_t i = 0; i < derivs.size(); ++i) {
//            assert(!std::isnan(derivs[i]));
//        }
    }

    float Derivatives::Dwrt(Float const &f, Float const &x) const {
        // Check if we computed the derivatives for the given out variable
        auto it = var_derivatives_map.find(f);
        if (it == var_derivatives_map.end()) {
            std::cout << "Could not find given output variable!" << std::endl;
            return 0.f;
        } else {
            return it->second[x.NodeIndex()];
        }
    }

} // drdemo namespace