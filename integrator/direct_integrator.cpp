//
// Created by simon on 13.05.17.
//

#include "direct_integrator.hpp"

namespace drdemo {

    Spectrum DirectIntegrator::IncomingRadiance(Ray const &ray, Scene const &scene, size_t depth) const {
        // Final ray incoming radiance
        Spectrum L;

        // Find closes interaction
        Interaction interaction;
        if (!scene.Intersect(ray, &interaction)) { return L; }

        // Compute shading at hit point
        for (auto const &light : scene.GetLights()) {
            Vector3F wi;
            Float pdf;
            const Spectrum Li = light->SampleLi(interaction, 0.f, 0.f, &wi, &pdf);
            L += Li * Abs(Dot(interaction.n, wi)) / pdf;
        }

        return L;
    }

} // drdemo namespace