//
// Created by Simon on 22.05.2017.
//

#ifndef DRDEMO_DIRECTIONAL_LIGHT_HPP
#define DRDEMO_DIRECTIONAL_LIGHT_HPP

#include "light.hpp"

namespace drdemo {

    /**
     * Define DirectionalLight class
     */
    class DirectionalLight : public LightInterface {
    private:
        // Point light direction
        Vector3f direction;
        // Intensity
        Spectrum intensity;

    public:
        DirectionalLight(Vector3F const &d, Spectrum const &i);

        Spectrum SampleLi(Interaction const &interaction, float u0, float u1,
                          Vector3F *const wi, Float *pdf) const override;
    };

} // drdemo namespace

#endif //DRDEMO_DIRECTIONAL_LIGHT_HPP
