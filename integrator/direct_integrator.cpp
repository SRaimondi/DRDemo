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
        if (scene.GetLights().empty()) {
            // No lights, only take albedo into account
            const Float n_dot_l = Clamp(Dot(interaction.n, -camera.LookDir()), Float(0.f), Float(1.f));
            L = n_dot_l * interaction.albedo;
        } else {
            for (const auto &light : scene.GetLights()) {
                for (int s = 0; s < light->NumSamples(); ++s) {
                    Vector3F wi;
                    Float pdf;
                    const Spectrum Li = scene.GetLights()[0]->SampleLi(interaction, 0.f, 0.f, &wi, &pdf);
                    const Float n_dot_l = Clamp(Dot(interaction.n, wi), Float(0.f), Float(1.f));
                    if (!Li.IsBlack() && pdf != 0.f && n_dot_l > 0.f) {
                        L += n_dot_l * interaction.albedo * Li / pdf;
                    }
                }
                // Scale given the number of samples
                L = L / (float) light->NumSamples();
            }
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