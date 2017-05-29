//
// Created by Simon on 22.05.2017.
//

#include "directional_light.hpp"

namespace drdemo {

    DirectionalLight::DirectionalLight(Vector3F const &d, Spectrum const &i)
            : direction(Normalize(d)), intensity(i) {}

    Spectrum
    DirectionalLight::SampleLi(Interaction const &, float, float,
                               Vector3F *const wi, Float *pdf) const {
        *wi = direction;
        *pdf = 1.f;

        return intensity;
    }

} // drdemo namespace
