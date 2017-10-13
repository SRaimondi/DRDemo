//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_INTERACTION_HPP
#define DRDEMO_INTERACTION_HPP

#include "geometry.hpp"
#include "spectrum.hpp"

namespace drdemo {

    /**
     * Define Interaction class, used to move around the information of the interaction between a Ray and a Shape
     */
    class Interaction {
    public:
        Interaction() = default;

        // Hit point
        Vector3F p;
        // Normal
        Vector3F n;
        // Intersection parameter
        Float t;
        // Outgoing direction, in world space
        Vector3F wo;
        // Albedo value
        Spectrum albedo;
    };

} // drdemo namespace

#endif //DRDEMO_INTERACTION_HPP
