//
// Created by Simone Raimondi on 07.08.17.
//

#include "scalar_function.hpp"

namespace drdemo {

    float GradientNorm2(const std::vector<float> &gradient) {
        float norm2 = 0.f;

        for (auto const v : gradient) {
            norm2 += v * v;
        }

        return norm2;
    }

    float GradientNorm(const std::vector<float> &gradient) {
        return std::sqrt(GradientNorm2(gradient));
    }

} // drdemo namespace