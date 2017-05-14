//
// Created by simon on 11.05.17.
//

#ifndef DRDEMO_SPHERE_HPP
#define DRDEMO_SPHERE_HPP

#include "shape.hpp"

namespace drdemo {

    /**
     * Define Sphere shape class
     */
    class Sphere : public Shape {
    private:
        // Sphere position
        Vector3F center;
        // Sphere radius
        Float radius;

    public:
        Sphere(Vector3F const &c, Float const &r);

        // Shape methods
        bool Intersect(Ray const &ray, Interaction *const interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        // Differentiable object methods
        std::vector<Float const *> GetDiffVariables() const override;

        void UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index = 0);
    };

} // drdemo namespace

#endif //DRDEMO_SPHERE_HPP
