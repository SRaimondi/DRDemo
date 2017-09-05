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
#include <triangle_mesh.hpp>
#include "sdf_sphere.hpp"

namespace drdemo {

    void RadiusSphereTest(int grid_res, float start_radius, float target_radius) {
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
                    grid->operator()(x, y, z) = Length(p) - start_radius;
                }
            }
        }

        /**
         * Create scene and add target sphere of different radius
         */
        Scene scene;

        // Add sphere
        scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(target_radius)));
        // Add lights
        scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        // Add first camera
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, -5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));


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
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.7f, 10e-7f, true);

        std::cout << "Done with first minimization, refining grid and restarting minimization..." << std::endl;
        int new_dims[3] = {grid_dims[0] * 2, grid_dims[1] * 2, grid_dims[2] * 2};
        grid->Refine(new_dims);

        // Rebind variables
        energy.RebindVars();

        // Minimise energy again
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-7f, true);

        render->RenderImage(&target, scene, *cameras[0]);
        tonemapper.Process("final.png", target);

        // At this point, since we know what should be the exact final values in the SDF, we can compute the error
        float final_error = 0.f;
        for (int z = 0; z < grid_dims[2]; z++) {
            for (int y = 0; y < grid_dims[1]; y++) {
                for (int x = 0; x < grid_dims[0]; x++) {
                    // Compute point coordinates
                    const Vector3f p = grid->CoordsAt(x, y, z);
                    // Get correct value
                    const float target_value = Length(p) - target_radius;
                    // Get initial value
                    const float initial_value = Length(p) - start_radius;

                    // Accumulate error
                    final_error += std::pow(target_value - grid->operator()(x, y, z).GetValue(), 2.f);
                }
            }
        }

        std::cout << "Final error on reconstruction: " << final_error << std::endl;
    }

    void RadiusSphereTestMR(int start_grid_res, float start_radius, float target_radius, int ref_steps) {
        // Target render size
        const size_t WIDTH = 256;
        const size_t HEIGHT = 256;

        // Maximum iterations
        const size_t MAX_ITERS = 100;

        /**
         * Create new SDF grid and initialise it as a sphere of radius 1.f
         */
        int grid_dims[3] = {start_grid_res, start_grid_res, start_grid_res};
        auto grid = std::make_shared<SignedDistanceGrid>(grid_dims[0], grid_dims[1], grid_dims[2],
                                                         BBOX(Vector3f(-2.f, -2.f, -2.f), Vector3f(2.f, 2.f, 2.f)));

        // Initialize grid using sphere of radius 1 as SDF
        for (int z = 0; z < grid_dims[2]; z++) {
            for (int y = 0; y < grid_dims[1]; y++) {
                for (int x = 0; x < grid_dims[0]; x++) {
                    // Compute point coordinates
                    const Vector3f p = grid->CoordsAt(x, y, z);
                    // Set SDF value
                    grid->operator()(x, y, z) = Length(p) - start_radius;
                }
            }
        }

        /**
         * Create scene and add target sphere of different radius
         */
        Scene scene;

        // Add sphere
        scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(target_radius)));

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        // Add cameras
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, -5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(-5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));


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
        int new_dims[3];

        // Minimise energy
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.7f, 10e-7f, true);

        std::cout << "Starting with additional refinement steps..." << std::endl;
        for (int step = 0; step < ref_steps; step++) {
            // Compute new grid resolution, 2 * current
            for (int i = 0; i < 3; i++) { new_dims[i] = grid->Size(i) * 2; }
            // Refine grid
            grid->Refine(new_dims);
            // Rebind variables
            energy.RebindVars();
            // Minimise energy again
            GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-7f, true);
        }

        render->RenderImage(&target, scene, *cameras[0]);
        tonemapper.Process("final.png", target);

        // At this point, since we know what should be the exact final values in the SDF, we can compute the error
        float final_error = 0.f;
        for (int z = 0; z < grid_dims[2]; z++) {
            for (int y = 0; y < grid_dims[1]; y++) {
                for (int x = 0; x < grid_dims[0]; x++) {
                    // Compute point coordinates
                    const Vector3f p = grid->CoordsAt(x, y, z);
                    // Get correct value
                    const float target_value = Length(p) - target_radius;
                    // Get initial value
                    const float initial_value = Length(p) - start_radius;

                    // Accumulate error
                    final_error += std::pow(target_value - grid->operator()(x, y, z).GetValue(), 2.f);
                }
            }
        }

        std::cout << "Final error on reconstruction: " << final_error << std::endl;
    }

    void MoveSphereTest(int grid_res) {
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

        // Add sphere, but move it slightly
        scene.AddShape(std::make_shared<Sphere>(Vector3F(0.5f, 0.f, 0.f), Float(1.f)));
        // Add lights
        scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        // Add first camera
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, -5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));


        // Add second camera
//        cameras.push_back(
//                std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
//                                                      60.f,
//                                                      WIDTH, HEIGHT));
//
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

    void OBJRenderTest(int grid_res, const std::string &file_name) {
        // Target render size
        const size_t WIDTH = 256;
        const size_t HEIGHT = 256;

        // Maximum iterations
        const size_t MAX_ITERS = 100;

        // Load mesh
        auto mesh = std::make_shared<TriangleMesh>(file_name);

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
         * Create scene and add mesh to render target views
         */
        Scene scene;
        // Add sphere
        scene.AddShape(mesh);
        // Add lights
        scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.6f, 0.8f, 1.f), Spectrum(0.9f)));

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(-3.54f, 0.f, -3.54f), Vector3F(),
                                                      Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));
        cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(3.06f, 2.5f, 3.06f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f, WIDTH, HEIGHT));

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

    void OBJRenderTestMR(int start_grid_res, const std::string &file_name, int ref_steps, const CAM_CONFIG &config) {
        // Target render size
        const size_t WIDTH = 256;
        const size_t HEIGHT = 256;

        // Maximum iterations
        const size_t MAX_ITERS = 200;

        // Load mesh
        auto mesh = std::make_shared<TriangleMesh>(file_name);

        /**
         * Create new SDF grid and initialise it as a sphere of radius 1.f
         */
        int grid_dims[3] = {start_grid_res, start_grid_res, start_grid_res};
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
         * Create scene and add mesh to render target views
         */
        Scene scene;
        // Add sphere
        scene.AddShape(mesh);

        // Create target_cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;

        if (config == CAM_CONFIG::NORMAL) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, -5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
        } else if (config == CAM_CONFIG::UP_DOWN) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(1.7678f, 3.5355f, 3.0619f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(3.0619f, -3.5355f, -1.7678f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-1.7678f, 3.5355f, -3.0619f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-3.0619f, -3.5355f, 1.7678f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
        } else if (config == CAM_CONFIG::ALL) {
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(1.7678f, 3.5355f, 3.0619f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(3.0619f, -3.5355f, -1.7678f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, -5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-1.7678f, 3.5355f, -3.0619f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
            cameras.push_back(
                    std::make_shared<const PinholeCamera>(Vector3F(-3.0619f, -3.5355f, 1.7678f), Vector3F(),
                                                          Vector3F(0.f, 1.f, 0.f),
                                                          60.f, WIDTH, HEIGHT));
        } else {
            std::cerr << "No valid camera configuration!" << std::endl;
        }

        /**
         * Create renderer class with a simple direct illumination integrator
         */
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Push status of tape before rendering target image
        default_tape.Push();

        // Render target images
        BoxFilterFilm target(WIDTH, HEIGHT);

        ClampTonemapper tonemapper;

        // Create target image
        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("target_" + std::to_string(i) + ".png", target);
        }

        std::vector<std::vector<float> > raw_views;
        default_tape.Push();
        // Render views
        for (auto const &camera : cameras) {
            render->RenderImage(&target, scene, *camera);
            raw_views.push_back(target.Raw());
        }
        default_tape.Pop();

        // Use SDF
        scene.ClearShapes();
        scene.AddShape(grid);

        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("start_" + std::to_string(i) + ".png", target);
        }

        // Remove from tape rendering variables
        default_tape.Pop();

        // Test with new energy
        auto energy = ReconstructionEnergy(scene, grid, raw_views, cameras, render, 1.f, WIDTH, HEIGHT);

        // Minimise energy
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.7f, 10e-8f, true);
        int new_dims[3];

        std::cout << "Starting with additional refinement steps..." << std::endl;
        for (int step = 0; step < ref_steps; step++) {
            // Compute new grid resolution, 2 * current
            for (int i = 0; i < 3; i++) { new_dims[i] = grid->Size(i) * 2; }
            // Refine grid
            grid->Refine(new_dims);
            // Rebind variables
            energy.RebindVars();
            // Minimise energy again
            GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-8f, true);
        }

        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("final_ " + std::to_string(i) + ".png", target);
        }
    }

} // drdemo namespace
