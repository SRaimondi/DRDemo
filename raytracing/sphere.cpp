//
// Created by simon on 20.04.17.
//

#include "sphere.hpp"

rt::Sphere::Sphere(float x_c, float y_c, float z_c, float r, rad::Tape<float> &tape) {
    // Set center
    center[0] = tape.NewVariable(x_c);
    center[1] = tape.NewVariable(y_c);
    center[2] = tape.NewVariable(z_c);
    // Set radius
    radius = tape.NewVariable(r);
}

std::vector<rad::Variable<float>> rt::Sphere::GetDifferentiableVariables() const {
    std::vector<rad::Variable<float>> vars;
    // Push into the vector the center and the radius
    vars.push_back(center[0]);
    vars.push_back(center[1]);
    vars.push_back(center[2]);
    vars.push_back(radius);

    return vars;
}

bool rt::Sphere::Intersect(const Ray &ray, Intersection *const intersection) const {
    using Float = ad::Float;
    // Compute terms for quadratic form
    Float a = utils::Dot(ray.direction, ray.direction);
    Float b = 2.f * utils::Dot(ray.origin, ray.direction);
    Float c = utils::Dot(ray.origin, ray.origin) - radius * radius;
    Float discr = b * b - 4.f * a * c;
    if (discr < 0.f) { return false; }
    Float root = rad::Sqrt(discr);
    Float q;
    if (b < 0.f) {
        q = -0.5f * (b - root);
    } else {
        q = -0.5f * (b + root);
    }
    // Find the two roots
    Float t0 = q / a;
    Float t1 = c / q;
    // Find nearest
    Float t_hit = t0;
    if (t0 > t1) { t_hit = t1; }
    // Fill intersection
    intersection->p = ray(t_hit.Value());
    intersection->n = intersection->p - center;
    intersection->n /= rad::Sqrt(utils::Length2(intersection->n));
    intersection->t_hit = t_hit;

    return true;
}
