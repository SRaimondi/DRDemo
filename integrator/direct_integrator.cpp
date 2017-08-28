//
// Created by simon on 13.05.17.
//

#include "direct_integrator.hpp"

namespace drdemo {

    Spectrum DirectIntegrator::IncomingRadiance(Ray const &ray, Scene const &scene, const CameraInterface &camera,
                                                size_t depth) const {
        // Final ray incoming radiance
        Spectrum L;

        // Find closes interaction
        Interaction interaction;
        if (!scene.Intersect(ray, &interaction)) { return L; }

            // TODO Testing if we get better result when lights come from the camera
            const Float n_dot_l = Abs(Dot(interaction.n, -camera.LookDir()));
            L = Spectrum(n_dot_l, n_dot_l, n_dot_l);

        // Compute shading at hit point
            // for (auto const &light : scene.GetLights()) {
            //     if (light->IsEnabled()) {
            // Vector3F wi;
            // Float pdf;
            // const Spectrum Li = light->SampleLi(interaction, 0.f, 0.f, &wi, &pdf);
            // L += Li * Abs(Dot(interaction.n, wi)) / pdf;
            //     }
            // }

        return L;
    }

} // drdemo namespace