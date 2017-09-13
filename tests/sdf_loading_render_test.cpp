//
// Created by simon on 14.09.17.
//

#include <grid.hpp>
#include <scene.hpp>
#include <clamp_tonemapper.hpp>
#include <box_film.hpp>
#include <direct_integrator.hpp>
#include <simple_renderer.hpp>
#include <pinhole_camera.hpp>
#include "sdf_loading_render_test.hpp"

namespace drdemo {

    void LoadAndTestSDF(const std::string &sdf_file_name, size_t w, size_t h) {
        // Create SDF from file
        auto sdf_grid = std::make_shared<SignedDistanceGrid>(sdf_file_name);

        // Create scene
        Scene scene;
        scene.AddShape(sdf_grid);

        // Create camera
        auto camera = PinholeCamera(Vector3F(1.f, 2.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60, w, h);

        // Disable tape
        default_tape.Disable();

        // Create renderer class with direct illumination integrator
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Create film and tonemapper
        BoxFilterFilm target(w, h);
        ClampTonemapper tonemapper;

        // Render camera
        render->RenderImage(&target, scene, camera);

        // Output image
        tonemapper.Process("sdf_render_test.png", target);

        // Re-enable tape
        default_tape.Enable();
    }

} // drdemo namespace