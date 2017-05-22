//
// Created by simon on 08.05.17.
//

#include <scene.hpp>
#include <derivative.hpp>
#include <sphere.hpp>
#include <pinhole_camera.hpp>
#include <simple_renderer.hpp>
#include <direct_integrator.hpp>
#include <box_film.hpp>
#include <clamp_tonemapper.hpp>
#include <directional_light.hpp>

#define WIDTH   512
#define HEIGHT  512

// Compute gradient norm, considering the gradient as a float vector
float GradNorm(std::vector<float> const &grad) {
    float norm = 0.f;
    for (auto const &e : grad) {
        norm += e * e;
    }

    return std::sqrt(norm);
}


//// Spherical function, minimization test
//drdemo::Float Spherical(std::vector<drdemo::Float> const &x) {
//    drdemo::Float res = x[0] * x[0];
//    for (size_t i = 1; i < x.size(); i++) {
//        res += x[i] * x[i];
//    }
//
//    return res;
//}

//// Matyas function
//Float Matyas(std::vector<Float> const &x) {
//    return 0.26f * (x[0] * x[0] + x[1] * x[1]) - 0.48f * x[0] * x[1];
//}

int main(void) {

    // Set namespace used
    using namespace drdemo;

//    // Minimization test
//    std::vector<Float> x(10);
//
//    // Initialize with some random data
//    for (size_t i = 0; i < x.size(); i++) {
//        x[i] = 6.f;
//    }
//
//    Derivatives derivatives;
//
//    std::vector<float> gradient(x.size(), 0.f);
//    size_t iters = 0;
//    float delta = 0.1f;
//
//    // Get index of the variable we need to keep
//    // size_t clear_index = default_tape.Size();

//    do {
//        // Push current index of tape to remove nodes after
//        default_tape.Push();
//        // Clear derivatives
//        derivatives.Clear();
//        // Compute spherical function value
//        Float y = Spherical(x);
//        // Float y = Matyas(x);
//        // Compute derivatives
//        derivatives.ComputeDerivatives(y);
//        // Compute gradient and update
//        for (int i = 0; i < x.size(); i++) {
//            gradient[i] = derivatives.Dwrt(y, x[i]);
//            x[i].SetValue(x[i].GetValue() - delta * gradient[i]);
//        }
//        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//        // Increase number of iterations
//        iters++;
//        // Clear tape
//        // default_tape.Clear(clear_index);
//        default_tape.Pop();
//        std::cout << "Tape size: " << default_tape.Size() << std::endl;
//    } while (GradNorm(gradient) > 0.0001f && iters < 10000);
//
//    std::cout << std::endl << "Final gradient norm: " << GradNorm(gradient) << std::endl;
//    std::cout << "Iterations: " << iters << std::endl;
//
//    std::cout << "#### Final x ####" << std::endl;
//    for (auto const & x_i : x) {
//        std::cout << x_i << " ";
//    }
//    std::cout << std::endl;

    // Derivatives computation class
    Derivatives derivatives;

    // Disable tape to render target image
    default_tape.Disable();

    // Create scene
    Scene scene;

    // Add sphere
    scene.AddShape(std::make_shared<Sphere>(Vector3F(2.f, 0.f, 0.f), Float(2.f)));
    // Add lights
    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));

    // Create camera
    auto camera = PinholeCamera(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f, WIDTH, HEIGHT);

    // Create renderer
    auto render = SimpleRenderer(std::make_shared<DirectIntegrator>());

    // Render target image
    BoxFilterFilm target(WIDTH, HEIGHT);
    render.RenderImage(&target, scene, camera);

    // Create tonemapper and process target image
    ClampTonemapper tonemapper;
    tonemapper.Process("target.ppm", target);

    // Re-enable tape
    default_tape.Enable();

    // Change sphere position
    scene.ClearShapes();
    scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(2.f)));
    //scene.AddShape(std::make_shared<Sphere>(Vector3F(-2.f, 0.f, 0.f), Float(2.f)));

    // Render start image
    BoxFilterFilm x(WIDTH, HEIGHT);
    render.RenderImage(&x, scene, camera);

    // Compute squared norm of image
    Float x_2_norm = x.SquaredNorm();

    derivatives.ComputeDerivatives(x_2_norm);

    std::vector<Float const *> vars;

    for (auto const & s : scene.GetShapes()) {
        s->GetDiffVariables(vars);
    }

    for (auto const & var : vars) {
        std::cout << derivatives.Dwrt(x_2_norm, *var) << std::endl;
    }

    // Process target image
    tonemapper.Process("x.ppm", x);



//    // Disable tape
//    // default_tape.Disable();
//
//    // Derivatives computation class
//    Derivatives derivatives;
//
//    // Create target film and starting
//    BoxFilterFilm target(WIDTH, HEIGHT);
//    BoxFilterFilm start(WIDTH, HEIGHT);
//
//    // Create tonemapper
//    ClampTonemapper tonemapper;
//
//    // Create scene
//    Scene scene;
//
//    // Create camera
//    std::shared_ptr<CameraInterface> camera =
//            std::make_shared<PinholeCamera>(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f),
//                                            60.f, WIDTH, HEIGHT);
//
//    // Add spheres
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(2.f, 0.f, 0.f), Float(2.f)));
//    // Add lights
//    scene.AddLight(std::make_shared<PointLight>(Vector3F(-5.f, 5.f, 3.f), Spectrum(10.f)));
//    scene.AddLight(std::make_shared<PointLight>(Vector3F(5.f, 5.f, 3.f), Spectrum(10.f)));
//
//    // Create renderer
//    std::shared_ptr<RendererInterface> renderer = std::make_shared<SimpleRenderer>(
//            std::make_shared<DirectIntegrator>());
//
//    default_tape.Push();
//    // Render target image
//    renderer->RenderImage(&target, scene, *camera.get());
//
//    // Create target image
//    tonemapper.Process("target.ppm", target);
//
//    default_tape.Pop();
//
//    // Clear scene's spheres
//    scene.ClearShapes();
//    // Add new sphere, in a different position with respect to the one used in the target image generation
//    std::shared_ptr<Shape> sphere = std::make_shared<Sphere>(Vector3F(1.f, 0.f, 0.f), Float(2.f));
//    // Get variables we can use to compute the gradient
//    std::vector<Float const *> vars = sphere->GetDiffVariables();
//    scene.AddShape(sphere);
//
//    // Render starting image
//    // renderer->RenderImage(&start, scene, *camera.get());
//    // tonemapper.Process("start.ppm", start);
//
//    // Create gradient storage vector
//    std::vector<float> gradient(4, 0.f);
//    std::vector<float> deltas(4, 0.f);
//    size_t iters = 0;
//    float delta = 0.001f;
//
//    // Save index of tape to clear after
//    // size_t clear_index = default_tape.Size();
//
//    do {
//        // default_tape.Push();
//        derivatives.Clear();
//        // Render new image
//        renderer->RenderImage(&start, scene, *camera.get());
//        // Compute image difference
//        BoxFilterFilm difference = target - start;
//        // Compute squared norm of images differences
//        Float sq_norm = difference.SquaredNorm();
//        std::cout << "Difference squared norm: " << sq_norm << std::endl;
//        // Compute derivatives
//        derivatives.ComputeDerivatives(sq_norm);
//        // Compute gradient
//        for (int i = 0; i < 4; i++) {
//            gradient[i] = derivatives.Dwrt(sq_norm, *vars[i]);
//            // Compute deltas to update variables
//            deltas[i] = -delta * gradient[i];
//        }
//        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//        // Update sphere
//        sphere->UpdateDiffVariables(deltas);
//        // Increase number of iterations
//        iters++;
//        // Create difference image
//        difference.Abs();
//        tonemapper.Process(std::string("difference_") + std::to_string(iters) + std::string(".ppm"), difference);
//        // Clear tape
//        // default_tape.Clear(clear_index);
//        // default_tape.Pop();
//        std::cout << "Tape size: " << default_tape.Size() << std::endl;
//    } while (GradNorm(gradient) > 0.01f);
//
//    std::cout << "Iterations: " << iters << std::endl;
//    std::cout << "Final norm: " << GradNorm(gradient) << std::endl;
//
//    // Create final image
//    BoxFilterFilm final(WIDTH, HEIGHT);
//    // Render new image
//    renderer->RenderImage(&final, scene, *camera.get());
//    tonemapper.Process("final.ppm", final);
    // tonemapper.Process("difference.ppm", difference);


    return EXIT_SUCCESS;
}
