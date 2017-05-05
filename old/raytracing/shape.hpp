//
// Created by simon on 20.04.17.
//

#ifndef DRDEMO_SHAPE_HPP
#define DRDEMO_SHAPE_HPP

/**
 * This file defines the base interface for a Shape that can be rendered
 */

#include "ray.hpp"
#include "intersection.hpp"

namespace rt {

    class Shape : public ad::RADDiffObjectInterface<float> {
    public:
        virtual ~Shape() {}

        // Compute intersection of a Ray with the Shape
        virtual bool Intersect(Ray const &ray, Intersection *const intersection) const = 0;
    };

} // rt namespace

#endif //DRDEMO_SHAPE_HPP
