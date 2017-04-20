//
// Created by simon on 18.04.17.
//

#ifndef DRDEMO_RAY_HPP
#define DRDEMO_RAY_HPP

/**
 * This file contains the definition of the Ray class used in the rest of the project
 */

#include "diffobject.hpp"

namespace rt {

    class Ray {
    public:
        // Ray origin and direction
        ad::Vec3F origin;
        ad::Vec3F direction;

        Ray(ad::Vec3F const &o, ad::Vec3F const &d);

        // Compute ray point at given parameter
        inline ad::Vec3F operator()(float t) const {
            return origin + t * direction;
        }
    };

    Ray::Ray(ad::Vec3F const &o, ad::Vec3F const &d)
            : origin(o), direction(d) {}

} // rt namespace

#endif //DRDEMO_RAY_HPP
