//
// Created by Simon on 24.08.2017.
//

#include <grid.hpp>
#include <scene.hpp>
#include <sphere.hpp>
#include <directional_light.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include <box_film.hpp>
#include <clamp_tonemapper.hpp>
#include <reconstruction_energy.hpp>
#include <gradient_descent.hpp>
#include "sdf_sphere.hpp"

namespace drdemo {

    void RadiusSphereTest(int grid_res, float radius) {
        // Target render size
        const size_t WIDTH = 256;
        const size_t HEIGHT = 256;

        // Maximum iterations
        const size_t MAX_ITERS = 100;

        /**
         * Create new SDF grid and initialise it as a sphere of radius 1.f
         */
        int grid_dims[3] = {grid_res, grid_res, grid_res};
        auto grid = std::make_shared<SignedDistanceGrid>(grid_dims[0], grid_dims[1], grid_dims[2],
                                                         BBOX(Vector3f(-2.f, -2.f, -2.f), Vector3f(2.f, 2.f, 2.f)));

        // Initialize grid using sphere of radius 1 as SDF
        for (int z = 0; z < grid_dims[2]; z++) {
            for (int y = 0; y < grid_dims[1]; y++) {
                for (int x = 0; x < grid_dims[0]; x++) {
                    // Compute point coordinates
                    const Vector3f p = grid->CoordsAt(x, y, z);
                    // Set SDF value
                    grid->operator()(x, y, z) = Length(p) - 1.f;
                }
            }
        }

        /**
         * Create scene and add target sphere of different radius
         */
        Scene scene;

        // Add sphere
        scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(radius)));
        // Add lights
        scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));
        //scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(1.f, 0.f, 0.f), Spectrum(0.9f)));

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        // Add first camera
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));

        // Add second camera
//        cameras.push_back(
//                std::make_shared<const PinholeCamera>(Vector3F(-5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
//                                                      60.f,
//                                                      WIDTH, HEIGHT));

        /**
         * Create renderer class with a simple direct illumination integrator
         */
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Push status of tape before rendering target image
        default_tape.Push();

        // Render target images
        BoxFilterFilm target(WIDTH, HEIGHT);

        // Create target image
        render->RenderImage(&target, scene, *cameras[0]);
        ClampTonemapper tonemapper;
        tonemapper.Process("target.png", target);

        std::vector<std::vector<float> > raw_views;
        // Render views
        for (auto const &camera : cameras) {
            render->RenderImage(&target, scene, *camera);
            raw_views.push_back(target.Raw());
        }

        // Use SDF
        scene.ClearShapes();
        scene.AddShape(grid);

        render->RenderImage(&target, scene, *cameras[0]);
        tonemapper.Process("start.png", target);

        // Remove from tape rendering variables
        default_tape.Pop();

        // Test with new energy
        auto energy = ReconstructionEnergy(scene, grid, raw_views, cameras, render, 1.f, WIDTH, HEIGHT);

        // Minimise energy
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-10f, true);

        render->RenderImage(&target, scene, *cameras[0]);
        tonemapper.Process("final.png", target);
    }

    void MoveSphereTest(int grid_res) {

    }

}
