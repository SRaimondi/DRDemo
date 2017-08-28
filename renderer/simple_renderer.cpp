//
// Created by simon on 12.05.17.
//

#include <iostream>
#include "simple_renderer.hpp"

// Ray passes thorough the center of the pixel
const float s_x = 0.5f;
const float s_y = 0.5f;

namespace drdemo {

    SimpleRenderer::SimpleRenderer(std::shared_ptr<const SurfaceIntegratorInterace> const &s_i)
            : surface_integrator(s_i) {}

    void SimpleRenderer::RenderImage(Film *const film, Scene const &scene,
                                     CameraInterface const &camera) const {
        // Current Ray
        Ray ray;
        // Incoming radiance
        Spectrum Li;

        for (size_t i = 0; i < film->Width(); i++) {
            for (size_t j = 0; j < film->Height(); j++) {
                // Generate ray
                ray = camera.GenerateRay(i, j, s_x, s_y);
                // Compute incoming radiance
                Li = surface_integrator->IncomingRadiance(ray, scene, camera, 0);
                // Add sample
                if (!film->AddSample(Li, i, j, s_x, s_y)) {
                    std::cerr << "Error adding sample to film!" << std::endl;
                }
            }
        }
    }

} // drdemo namespace
