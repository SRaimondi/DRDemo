//
// Created by Simon on 28.08.2017.
//

#include "ambient_light.hpp"

namespace drdemo {

    AmbientLight::AmbientLight(const Spectrum &i) : intensity(i) {}

    Spectrum
    AmbientLight::SampleLi(Interaction const &interaction, float u0, float u1, Vector3F *wi, Float *pdf) const {

        // Set as interaction normal
        wi->x = interaction.n.x.GetValue();
        wi->y = interaction.n.y.GetValue();
        wi->z = interaction.n.z.GetValue();

        *pdf = 1.f;
        return intensity;
    }

} // drdemo namespace