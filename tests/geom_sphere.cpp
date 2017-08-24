//
// Created by Simon on 24.08.2017.
//

#include <scene.hpp>
#include <sphere.hpp>
#include <directional_light.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include <box_film.hpp>
#include <clamp_tonemapper.hpp>
#include <multi_view_energy.hpp>
#include <gradient_descent.hpp>
#include "geom_sphere.hpp"

namespace drdemo {

    void GeometricSphereTest() {

        // Create scene
        Scene scene;

        // Target render size
        const size_t WIDTH = 256;
        const size_t HEIGHT = 256;

        // Maximum iterations
        const size_t MAX_ITERS = 20;

        /**
         * Create the starting scene that we render to create the target images for our problem
         */
        // Add sphere
        scene.AddShape(std::make_shared<Sphere>(Vector3F(1.f, 0.f, 0.f), Float(2.f)));
        // Add lights
        scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));

        /**
         * Create the cameras we use in our rendering process
         */
        std::vector<std::shared_ptr<const CameraInterface>> target_cameras;
        // Add first camera
        target_cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(-10.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));
        // Add second camera
        target_cameras.push_back(
                std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                                      60.f,
                                                      WIDTH, HEIGHT));

        /**
         * Create renderer class with a simple direct illumination integrator
         */
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Push status of tape before rendering target image
        default_tape.Push();

        // Create film
        BoxFilterFilm target(WIDTH, HEIGHT);

        // Tonemapper used to process the film and create images
        ClampTonemapper tonemapper;

        // Render one target image
        render->RenderImage(&target, scene, *target_cameras[1]);
        tonemapper.Process("target.png", target);

        // Render target images
        std::vector<std::vector<float>> raw_views;
        // Render views
        for (auto const &camera : target_cameras) {
            render->RenderImage(&target, scene, *camera);
            raw_views.push_back(target.Raw());
        }

        // Remove from tape rendering variables
        default_tape.Pop();

        // Change sphere position to center, reduce radius and try to match the images
        scene.ClearShapes();
        scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(1.f)));

        // Render start image
        render->RenderImage(&target, scene, *target_cameras[1]);
        tonemapper.Process("start.png", target);

        // Create multi-view energy
        auto energy = MultiViewEnergy(scene, raw_views, target_cameras, render, WIDTH, HEIGHT);

        // Minimise energy
        GradientDescentBT::Minimize(energy, MAX_ITERS, 0.5f, 0.5f, 0.8f, 10e-10f, true);

        // Render final view
        render->RenderImage(&target, scene, *target_cameras[1]);
        tonemapper.Process("final.png", target);

        // Print data of the sphere at the end, shoudl match the one at line 36
        std::cout << "Final sphere data" << std::endl;
        std::cout << scene.GetShapes()[0]->ToString() << std::endl;
    }

} // drdemo namespace
