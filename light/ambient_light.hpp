//
// Created by Simon on 28.08.2017.
//

#ifndef DRDEMO_AMBIENT_LIGHT_HPP
#define DRDEMO_AMBIENT_LIGHT_HPP

#include "light.hpp"

namespace drdemo {

    /**
     * Define ambient light class
     */
    class AmbientLight : public LightInterface {
    private:
        // Ambient intensity
        Spectrum intensity;

    public:
        explicit AmbientLight(const Spectrum &i);

        Spectrum SampleLi(Interaction const &interaction, float u0, float u1,
                          Vector3F *wi, Float *pdf) const override;
    };

} // drdemo namespace

#endif //DRDEMO_AMBIENT_LIGHT_HPP
