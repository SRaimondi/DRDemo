//
// Created by simon on 18.04.17.
//

#ifndef DRDEMO_DIFFOBJECT_HPP
#define DRDEMO_DIFFOBJECT_HPP

/**
 * This file defines the interface for a differentiable object.
 * The idea behind the interface is that this object is involved in some computations and it holds some variable
 * for which we may want to compute a derivative with respect to it at the end, so we must be able to access them later
 * and they must be initialized thorough a tape in case of reverse automatic differentiation
 */

#include "vector.hpp"

namespace ad {

    /**
     * Define the interface of an object that can be differentiated, using reverse automatic differentiation
     */
    template<typename T>
    class RADDiffObjectInterface {
    public:
        // Return a std::vector with all the variables that can be used to compute a derivative
        virtual std::vector<rad::Variable<T> > GetDifferentiableVariables() const = 0;
    };

    /**
     * Define custom Float and Double
     */
    using Float = rad::Variable<float>;
    using Double = rad::Variable<double>;

    /**
     * Define custom Vector3 types
     */
    using Vec3F = utils::Vector3<Float>;
    using Vec3f = utils::Vector3<float>;
    using Vec3D = utils::Vector3<Double>;
    using Vec3d = utils::Vector3<double>;

} // ad namespace

#endif //DRDEMO_DIFFOBJECT_HPP
