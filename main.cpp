//
// Created by simon on 08.05.17.
//

#include <scene.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <sphere.hpp>
#include <point_light.hpp>
#include <renderer.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include <derivative.hpp>
#include "clamp_tonemapper.hpp"
#include "box_film.hpp"

#define WIDTH 512
#define HEIGHT 512

int main(void) {

    // Set namespace used
    using namespace drdemo;

    // Disable tape
    default_tape.Disable();

    // Create target film and starting
    BoxFilterFilm target(WIDTH, HEIGHT);
    BoxFilterFilm start(WIDTH, HEIGHT);

    // Create tonemapper
    ClampTonemapper tonemapper;

    // Create scene
    Scene scene;

    // Create camera
    std::shared_ptr<CameraInterface> camera =
            std::make_shared<PinholeCamera>(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                            60.f, WIDTH, HEIGHT);

    // Add spheres
    scene.AddShape(std::make_shared<Sphere>(Vector3F(2.f, 0.f, 0.f), Float(2.f)));
    // Add lights
    scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, 3.f), Spectrum(10.f)));
    scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, 3.f), Spectrum(10.f)));

    // Create renderer
    std::shared_ptr<RendererInterface> renderer = std::make_shared<SimpleRenderer>(
            std::make_shared<DirectIntegrator>());

    // Render target image
    renderer->RenderImage(&target, scene, *camera.get());

    // Clear scene's spheres
    scene.ClearShapes();
    // Add new sphere
    scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f), Float(2.f)));
    // Render starting image
    renderer->RenderImage(&start, scene, *camera.get());

    // Compute image difference
    BoxFilterFilm difference = target - start;
    // print Squared Norm of difference
    std::cout << difference.SquaredNorm() << std::endl;
    // Compute abs
    difference.Abs();

    // Create images
    tonemapper.Process("target.ppm", target);
    tonemapper.Process("start.ppm", start);
    tonemapper.Process("difference.ppm", difference);


    // Create derivatives
    // Derivatives derivatives;
    // derivatives.ComputeDerivatives(squared_norm);

    // std::cout << default_tape.Size() << std::endl;

    return EXIT_SUCCESS;
}
