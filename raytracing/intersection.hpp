//
// Created by simon on 20.04.17.
//

#ifndef DRDEMO_INTERSECTION_HPP
#define DRDEMO_INTERSECTION_HPP

/**
 * This file contains the definition of the Intersection struct, which is used as a pass through of the information
 * of the intersection of a Shape with a Ray
 */

#include "diffobject.hpp"

namespace rt {

    struct Intersection {
        // Intersection point
        ad::Vec3F p;
        // Intersection normal
        ad::Vec3F n;
        // Intersection paramter
        ad::Float t_hit;
    };

}

#endif //DRDEMO_INTERSECTION_HPP
