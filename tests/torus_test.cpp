//
// Created by Simon on 07.09.2017.
//

#include <cstdlib>
#include <triangle_mesh.hpp>
#include <grid.hpp>
#include <scene.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <direct_integrator.hpp>
#include <simple_renderer.hpp>
#include <box_film.hpp>
#include <clamp_tonemapper.hpp>
#include <reconstruction_energy.hpp>
#include <gradient_descent.hpp>
#include <reconstruction_energy_opt.hpp>
#include "torus_test.hpp"
#include "test_common.hpp"

namespace drdemo {

    void TorusTest(int start_resolution, float res_multiplier, int ref_steps) {
        // Target render size
        const size_t WIDTH = 512;
        const size_t HEIGHT = 512;

        // Maximum gradient descent iterations
        const size_t MAX_ITERS = 200;

        // Load torus mesh
        auto mesh = std::make_shared<TriangleMesh>("../objs/torus.obj");

        // Create SDF grid
        int grid_dims[3] = {start_resolution, start_resolution, start_resolution};
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

        // Create scene and add mesh
        Scene scene;
        scene.AddShape(mesh);

        // Create vector of target cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;
        // Load viewpoints from file
        std::vector<Vector3f> cameras_viewpoints = LoadViewPoints("../camera_points/normal_camera.txt", "%f,%f,%f");
        for (auto const &p : cameras_viewpoints) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(p.x, p.y, p.z), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
        }

        // Create renderer class with direct illumination integrator
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Disable tape
        default_tape.Disable();

        BoxFilterFilm target(WIDTH, HEIGHT);
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
        scene.AddShape(grid);

        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("start_" + std::to_string(i) + ".png", target);
        }

        default_tape.Enable();

        // Create energy
        // auto energy = ReconstructionEnergy(scene, grid, raw_views, cameras, render, 1.f, WIDTH, HEIGHT);
        auto energy = ReconstructionEnergyOpt(scene, grid, raw_views, cameras, render, 1.f, WIDTH, HEIGHT);

        // Do first minimisation
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-10f, true);

        // Start refinement
        int new_dims[3];
        for (int step = 0; step < ref_steps; step++) {
            std::cout << "Starting refinement step " << std::to_string(step + 1) << " of " << std::to_string(ref_steps)
                      << std::endl;
            // Compute new grid resolution
            for (int i = 0; i < 3; i++) { new_dims[i] = (int) (grid->Size(i) * res_multiplier); }
            std::cout << "Grid resolution: " << new_dims[0] << "x" << new_dims[1] << "x" << new_dims[2] << std::endl;
            // Refine grid
            grid->Refine(new_dims);
            // Rebind variables
            energy.RebindVars();
            // Minimise energy again
            GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-10f, true);
        }

        // Render final SDF status
        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("final_ " + std::to_string(i) + ".png", target);
        }
    }

} // drdemo namespace
