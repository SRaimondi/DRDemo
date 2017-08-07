//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_LIGHT_HPP
#define DRDEMO_LIGHT_HPP

#include "spectrum.hpp"
#include "interaction.hpp"

namespace drdemo {

    /**
     * Define LightInterface, all lights implementations must inherit from this one
     */
    class LightInterface {
    public:
        // Sample incoming light at a given Interaction, returns incoming radiance and fills sampling parameters
        virtual Spectrum
        SampleLi(Interaction const &interaction, float u0, float u1, Vector3F *wi, Float *pdf) const = 0;
    };

} // drdemo namespace


#endif //DRDEMO_LIGHT_HPP
