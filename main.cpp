//
// Created by simon on 08.05.17.
//

//#include <iostream>
//#include "derivative.hpp"
//#include "geometry.hpp"

#include <scene.hpp>
#include <camera.hpp>
#include <pinhole_camera.hpp>
#include <sphere.hpp>
#include <point_light.hpp>
#include <renderer.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include "clamp_tonemapper.hpp"
#include "box_film.hpp"

#define WIDTH 512
#define HEIGHT 512

int main(void) {

    // Set namespace used
    using namespace drdemo;

//    Float x(1.f);
//    Float y(2.f);
//
//    // Compute result
//    Float z = Exp(x) * Sin(y);
//    Float t(5.f);
//    Float f = t * t * z;
//
//    // Compute derivatives
//    Derivatives derivatives;
//    derivatives.ComputeDerivatives(z);
//    std::cout << "dz/dx: " << derivatives.Dwrt(z, x) << std::endl;
//    std::cout << "dz/dy: " << derivatives.Dwrt(z, y) << std::endl;
//    std::cout << std::endl;
//
//    derivatives.ComputeDerivatives(f);
//    std::cout << "df/dx: " << derivatives.Dwrt(f, x) << std::endl;
//    std::cout << "df/dy: " << derivatives.Dwrt(f, y) << std::endl;
//    std::cout << "df/dt: " << derivatives.Dwrt(f, t) << std::endl;
//    std::cout << std::endl;

    // Create film
    BoxFilterFilm film(WIDTH, HEIGHT);
    // Create camera
    std::shared_ptr<CameraInterface> camera =
            std::make_shared<PinholeCamera>(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
                                            60.f, WIDTH, HEIGHT);

    // Create scene
    Scene scene;
    // Add sphere
    scene.AddShape(std::make_shared<Sphere>(Vector3F(), Float(2.f)));
    // Add light
    // scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, -5.f), Spectrum(10.f)));
    scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, 5.f), Spectrum(10.f)));
    scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, 5.f), Spectrum(10.f)));
    // scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, -5.f), Spectrum(10.f)));

    // Create renderer
    std::shared_ptr<RendererInterface> renderer = std::make_shared<SimpleRenderer>(
            std::make_shared<DirectIntegrator>());

    // Render image
    renderer->RenderImage(&film, scene, *camera.get());

    // Create image
    ClampTonemapper tonemapper;
    tonemapper.Process(std::string("test.ppm"), film);

    return EXIT_SUCCESS;
}
