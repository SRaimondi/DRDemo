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

    default_tape.Enable();

    // Create film
    BoxFilterFilm film(WIDTH, HEIGHT);
    // Create camera
    std::shared_ptr<CameraInterface> camera =
            std::make_shared<PinholeCamera>(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                            60.f, WIDTH, HEIGHT);

    // Create scene
    Scene scene;
    // Add spheres
    scene.AddShape(std::make_shared<Sphere>(Vector3F(2.f, 0.f, 0.f), Float(2.f)));
    // Add light
    // scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, -5.f), Spectrum(10.f)));
    scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, 3.f), Spectrum(10.f)));
    scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, 3.f), Spectrum(10.f)));
    // scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, -5.f), Spectrum(10.f)));

    // Create renderer
    std::shared_ptr<RendererInterface> renderer = std::make_shared<SimpleRenderer>(
            std::make_shared<DirectIntegrator>());

    // Render image
    renderer->RenderImage(&film, scene, *camera.get());

    Float squared_norm = film.SquaredNorm();
    std::cout << squared_norm << std::endl;

    // Create derivatives
    Derivatives derivatives;
    derivatives.ComputeDerivatives(squared_norm);

    std::cout << default_tape.Size() << std::endl;

    // Create image
    ClampTonemapper tonemapper;
    tonemapper.Process("test.ppm", film);

    return EXIT_SUCCESS;
}
