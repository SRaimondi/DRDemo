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

        // TODO Testing if we get better result when lights comes from the camera
        const Float n_dot_l = Abs(Dot(interaction.n, -camera.LookDir()));
        if (scene.GetLights().empty()) {
            // No lights, only take albedo into account
            L = n_dot_l * interaction.albedo;
        } else {
            // Take albedo and light color into account, we suppose here that we are only dealing with the ambient color light
            // TODO This should be general for all lights!!!!
            Vector3F wi;
            Float pdf;
            const Spectrum Li = scene.GetLights()[0]->SampleLi(interaction, 0.f, 0.f, &wi, &pdf);
            L = n_dot_l * interaction.albedo * Li / pdf;
        }

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