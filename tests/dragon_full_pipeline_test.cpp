//
// Created by simon on 14.09.17.
//

#include <triangle_mesh.hpp>
#include <grid.hpp>
#include <scene.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <clamp_tonemapper.hpp>
#include <box_film.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include <reconstruction_energy_opt.hpp>
#include <gradient_descent.hpp>
#include "dragon_full_pipeline_test.hpp"
#include "test_common.hpp"

namespace drdemo {

    void FullPipelineTestDragon(size_t w, size_t h) {
        // Maximum gradient descent iterations
        const size_t MAX_ITERS = 250;

        // Load target mesh (the dragon file is NOT in the obj folder in the git repo)
        auto dragon_mesh = std::make_shared<TriangleMesh>("../objs/dragon.obj");

        // Load sdf
        auto sdf_grid = std::make_shared<SignedDistanceGrid>("../sdfs/dragon_mvs_output.sdf");

        // Create scene and add mesh
        Scene scene;
        scene.AddShape(dragon_mesh);

        // Create vector of target cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;
        // Load viewpoints from file
        std::vector<Vector3f> cameras_viewpoints = LoadViewPoints("../camera_points/up_hemi_dragon.txt", "%f,%f,%f");
        for (auto const &p : cameras_viewpoints) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(p.x, p.y, p.z), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, w, h));
        }

        // Create renderer class with direct illumination integrator
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());


        // Disable tape
        default_tape.Disable();

        BoxFilterFilm target(w, h);
        ClampTonemapper tonemapper;

        // Create target image
        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("target_" + std::to_string(i) + ".png", target);
        }

        // Raw target views
        std::vector<std::vector<float> > raw_views;
        // default_tape.Push();        // 2
        for (auto const &camera : cameras) {
            render->RenderImage(&target, scene, *camera);
            raw_views.push_back(target.Raw());
        }

        // Remove mesh and add SDF
        scene.ClearShapes();
        scene.AddShape(sdf_grid);

        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("start_" + std::to_string(i) + ".png", target);
        }

        default_tape.Enable();

        // Create energy
        auto energy = ReconstructionEnergyOpt(scene, sdf_grid, raw_views, cameras, render, 1.f, w, h);

        // Do first minimisation
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-10f, true);

        default_tape.Disable();
        // Render final SDF status
        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("final_ " + std::to_string(i) + ".png", target);
        }
        default_tape.Enable();
    }

} // drdemo namespace