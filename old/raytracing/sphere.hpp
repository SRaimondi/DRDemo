//
// Created by simon on 20.04.17.
//

#ifndef DRDEMO_SPHERE_HPP
#define DRDEMO_SPHERE_HPP

/**
 * Define Sphere that we can use to compute the derivative with respect to his center and radius
 */

#include "shape.hpp"

namespace rt {

    class Sphere : public Shape {
    private:
        // Sphere center
        ad::Vec3F center;
        // Radius
        ad::Float radius;

    public:
        // Constructor, needs a tape that is used to initialize the variable we can compute the derivative with respect to
        Sphere(float x_c, float y_c, float z_c, float r, rad::Tape<float> &tape);

        std::vector <rad::Variable<float>> GetDifferentiableVariables() const override;

        bool Intersect(Ray const &ray, Intersection *const intersection) const override;
    };

}

#endif //DRDEMO_SPHERE_HPP
