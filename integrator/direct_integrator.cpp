//
// Created by simon on 13.05.17.
//

#include "direct_integrator.hpp"

namespace drdemo {

    Spectrum DirectIntegrator::IncomingRadiance(Ray const &ray, Scene const &scene, size_t) const {
        // Final ray incoming radiance
        Spectrum L;

        // Find closes interaction
        Interaction interaction;
        if (!scene.Intersect(ray, &interaction)) { return L; }

        // Compute shading at hit point
        // TODO: Add different material support
        for (auto const &light : scene.GetLights()) {
            Vector3F wi;
            Float pdf;
            // FIXME: Add random number generation
            Spectrum Li = light->SampleLi(interaction, 0.f, 0.f, &wi, &pdf);
            // L += Li * Clamp(Dot(interaction.n, wi), Float(0.f), Float(1.f)) / pdf;
//            Float dot = Dot(interaction.n, wi);
//            assert(!std::isnan(dot.GetValue()));
//            assert(!std::isinf(dot.GetValue()));
//            assert(dot != 0.f);
//            assert(pdf != 0.f);
            L += Li * Abs(Dot(interaction.n, wi)) / pdf;
        }

        return L;
    }

} // drdemo namespace