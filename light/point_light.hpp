//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_POINT_LIGHT_HPP
#define DRDEMO_POINT_LIGHT_HPP

#include "light.hpp"

namespace drdemo {

    /**
     * Define PointLight class
     */
    class PointLight : public LightInterface {
    private:
        // Point light position
        Vector3F position;
        // Intensity
        Spectrum intensity;

    public:
        PointLight(Vector3F const &p, Spectrum const &i);

        Spectrum
        SampleLi(Interaction const &interaction, float u0, float u1, Vector3F *const wi, Float *pdf) const override;

    };

} // drdemo namespace

#endif //DRDEMO_POINT_LIGHT_HPP
