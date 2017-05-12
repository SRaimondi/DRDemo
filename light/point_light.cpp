//
// Created by simon on 12.05.17.
//

#include "point_light.hpp"

namespace drdemo {

    PointLight::PointLight(Vector3F const &p, Spectrum const &i)
            : position(p), intensity(i) {}

    Spectrum
    PointLight::SampleLi(Interaction const &interaction, float, float, Vector3F *const wi,
                         Float *pdf) const {
        *wi = Normalize(position - interaction.p);
        *pdf = 1.f;

        return (intensity / LengthSquared(position - interaction.p));
    }

} // drdemo namespace
