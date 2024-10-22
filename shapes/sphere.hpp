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
        bool Intersect(Ray const &ray, Interaction *interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        BBOX BBox() const override;

        Vector3f Centroid() const override;

        std::string ToString() const override;

        // Differentiable object methods
        void GetDiffVariables(std::vector<Float const *> &vars) const override;

        size_t GetNumVars() const noexcept override;

        void UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) override;

        void SetDiffVariables(const std::vector<float> &vals, size_t starting_index) override;
    };

} // drdemo namespace

#endif //DRDEMO_SPHERE_HPP
