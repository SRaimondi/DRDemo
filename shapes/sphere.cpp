//
// Created by simon on 11.05.17.
//

#include "sphere.hpp"

namespace drdemo {

    Sphere::Sphere(Vector3F const &c, Float const &r)
            : center(c), radius(r) {}

    bool Sphere::Intersect(Ray const &ray, Interaction *const interaction) const {
        Vector3F oc = ray.o - center;
        // Compute terms for quadratic form
        Float a = Dot(ray.d, ray.d);
        Float b = 2.f * Dot(oc, ray.d);
        Float c = Dot(oc, oc) - radius * radius;
        // Compute discriminant
        Float discr = b * b - 4.f * a * c;
        if (discr < 0.f) { return false; }
        // Compute root
        Float root = Sqrt(discr);
        Float q;
        if (b < 0.f) {
            q = -0.5f * (b - root);
        } else {
            q = -0.5f * (b + root);
        }
        // Find the two roots
        Float t0 = q / a;
        Float t1 = c / q;
        if (t0 > t1) {
            Float temp = t0;
            t0 = t1;
            t1 = temp;
        }
        // Check boundaries
        if (t0 > ray.t_max || t1 < ray.t_min) { return false; }
        Float t_hit = t0;
        if (t_hit < ray.t_min) {
            t_hit = t1;
            if (t_hit < ray.t_max) { return false; }
        }

        // Update new ray maximum value
        ray.t_max = t_hit;

        // Fill interaction
        interaction->p = ray(t_hit);
        interaction->n = Normalize(interaction->p - center);
        interaction->t = t_hit;
        interaction->wo = Normalize(-ray.d);

        return true;
    }

    bool Sphere::IntersectP(Ray const &ray) const {
        Vector3F oc = ray.o - center;
        // Compute terms for quadratic form
        Float a = Dot(ray.d, ray.d);
        Float b = 2.f * Dot(oc, ray.d);
        Float c = Dot(oc, oc) - radius * radius;
        // Compute discriminant
        Float discr = b * b - 4.f * a * c;
        if (discr < 0.f) { return false; }
        // Compute root
        Float root = Sqrt(discr);
        Float q;
        if (b < 0.f) {
            q = -0.5f * (b - root);
        } else {
            q = -0.5f * (b + root);
        }
        // Find the two roots
        Float t0 = q / a;
        Float t1 = c / q;
        if (t0 > t1) {
            Float temp = t0;
            t0 = t1;
            t1 = temp;
        }
        // Check boundaries
        if (t0 > ray.t_max || t1 < ray.t_min) { return false; }
        Float t_hit = t0;
        if (t_hit < ray.t_min) {
            t_hit = t1;
            if (t_hit < ray.t_max) { return false; }
        }

        return true;
    }

    void Sphere::GetDiffVariables(std::vector<Float const *> & vars) const {
        // Push sphere center variables
        vars.push_back(&(center.x));
        vars.push_back(&(center.y));
        vars.push_back(&(center.z));
        // Push radius variable
        vars.push_back(&(radius));
    }

    size_t Sphere::GetNumVars() const noexcept {
        return 4;
    }

    void Sphere::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {
        // Order is supposed to be sphere's center (x,y,z) and radius
        center.x.SetValue(center.x.GetValue() + delta[starting_index]);
        center.y.SetValue(center.y.GetValue() + delta[starting_index + 1]);
        center.z.SetValue(center.z.GetValue() + delta[starting_index + 2]);
        radius.SetValue(radius.GetValue() + delta[starting_index + 3]);
    }

} // drdemo namespace