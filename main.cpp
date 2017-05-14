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

// Compute gradient norm, considering the gradient as a float vector
float GradNorm(std::vector<float> const &grad) {
    float norm = 0.f;
    for (auto const &e : grad) {
        norm += e * e;
    }

    return std::sqrt(norm);
}

int main(void) {

    // Set namespace used
    using namespace drdemo;

    // Disable tape
    // default_tape.Disable();

    // Derivatives computation class
    Derivatives derivatives;

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

    // Create target image
    tonemapper.Process("target.ppm", target);

    // Clear scene's spheres
    scene.ClearShapes();
    // Add new sphere, in a different position with respect to the one used in the target image generation
    std::shared_ptr<Shape> sphere = std::make_shared<Sphere>(Vector3F(1.9f, 0.f, 0.f), Float(2.f));
    // Get variables we can use to compute the gradient
    std::vector<Float const *> vars = sphere->GetDiffVariables();
    scene.AddShape(sphere);

    // Create gradient storage vector
    std::vector<float> gradient(4, 0.f);
    std::vector<float> deltas(4, 0.f);
    size_t iters = 0;
    float delta = 0.0001f;

    // Save index of tape to clear after
    size_t clear_index = default_tape.Size();

    do {
        derivatives.Clear();
        // Render new image
        renderer->RenderImage(&start, scene, *camera.get());
        // Compute image difference
        BoxFilterFilm difference = target - start;
        // Compute squared norm of images differences
        Float sq_norm = difference.SquaredNorm();
        std::cout << "Difference squared norm: " << sq_norm << std::endl;
        // Compute derivatives
        derivatives.ComputeDerivatives(sq_norm);
        // Compute gradient
        for (int i = 0; i < 4; i++) {
            gradient[i] = derivatives.Dwrt(sq_norm, *vars[i]);
            // Compute deltas to update variables
            deltas[i] = -delta * gradient[i];
        }
        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
        // Update sphere
        sphere->UpdateDiffVariables(deltas);
        // Clear tape
        default_tape.Clear(clear_index);
        // Increase number of iterations
        iters++;
        // Create difference image
        tonemapper.Process(std::string("difference_") + std::to_string(iters) + std::string(".ppm"), difference);
    } while (GradNorm(gradient) > 0.01f && iters < 100);

    std::cout << "Iterations: " << iters << std::endl;

    // Create images
    // tonemapper.Process("start.ppm", start);
    // tonemapper.Process("difference.ppm", difference);


    // Create derivatives
    // Derivatives derivatives;
    // derivatives.ComputeDerivatives(squared_norm);

    // std::cout << default_tape.Size() << std::endl;

    return EXIT_SUCCESS;
}
