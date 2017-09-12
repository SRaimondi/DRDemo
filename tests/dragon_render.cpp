//
// Created by Simon on 12.09.2017.
//

#include <triangle_mesh.hpp>
#include <scene.hpp>
#include <camera.hpp>
#include <box_film.hpp>
#include <clamp_tonemapper.hpp>
#include <direct_integrator.hpp>
#include <simple_renderer.hpp>
#include <pinhole_camera.hpp>
#include "dragon_render.hpp"
#include "test_common.hpp"

namespace drdemo {

    void RenderDragonImages(const std::string &viewpoints_file, size_t w, size_t h) {
        // Disable tape
        default_tape.Disable();

        // Load dragon mesh
        auto mesh = std::make_shared<TriangleMesh>("../objs/dragon.obj");
        // Create scene and add mesh
        Scene scene;
        scene.AddShape(mesh);

        // Create target cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;
        // Load viewpoints file
        std::vector<Vector3f> cameras_viewpoints = LoadViewPoints(viewpoints_file, "%f %f %f", 430, 479);
        for (auto const &p : cameras_viewpoints) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(p.x, p.y, p.z), Vector3F(0.f, 0.f, 0.f),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, w, h));
        }

        // Create renderer class with direct illumination integrator
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Create film and tonemapper
        BoxFilterFilm target(w, h);
        ClampTonemapper tonemapper;

        // Render target image
        for (int i = 0; i < cameras.size(); ++i) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("view_" + std::to_string(i) + ".png", target);
        }

        // Re-enable tape
        default_tape.Enable();
    }

} // drdemo namespace