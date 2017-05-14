//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_SHAPE_HPP
#define DRDEMO_SHAPE_HPP

#include "interaction.hpp"
#include "diff_object.hpp"

namespace drdemo {

    /**
     * Define base Shape class, all other shapes must inherit from it
     */
    class Shape : public DiffObjectInterface {
    public:
        virtual ~Shape() {}

        // Compute intersection between a Ray and the shape
        virtual bool Intersect(Ray const &ray, Interaction *const interaction) const = 0;

        // Check for intersection between Ray and the shape
        virtual bool IntersectP(Ray const &ray) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_SHAPE_HPP
