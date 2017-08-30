//
// Created by simon on 08.05.17.
//

#include <scene.hpp>
#include <sdf_sphere.hpp>

#define WIDTH   512
#define HEIGHT  512

#define MAX_ITERS 10

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

int main() {

    // Set namespace used
    using namespace drdemo;

    // Geometric sphere test
    // GeometricSphereTest();

    // Radius sphere SDF test
    // RadiusSphereTest(10, 1.f, 1.2f);

    // RadiusSphereTestMR(5, 1.f, 1.5f, 2);

    // Move sphere test
    // MoveSphereTest(11);

    // Ellipse render test
    // OBJRenderTestMR(5, "../objs/ellipse.obj", 2);
    // OBJRenderTestMR(5, "../objs/sphere1_2.obj", 2);

    // Two sphere test
    // OBJRenderTestMR(7, "../objs/spheres2.obj", 3, false);   // See two sphere folder

    // Torus test
    OBJRenderTestMR(7, "../objs/torus.obj", 3);

    /**
     * Triangle mesh loading + simple minimization against black image
     */

//    // Derivatives computation class
//    Derivatives derivatives;
//
    // Try to load sphere mesh
//    auto mesh = std::make_shared<TriangleMesh>("../objs/cube.obj");
//
//    // Get list of triangles
//    // std::vector<std::shared_ptr<Shape> > triangles;
//    // mesh.CreateTriangles(triangles);
//
//    // Build BVH
//    // auto bvh = std::make_shared<BVH>(triangles);
//
//    // Create scene
//    Scene scene;
//
//    // Add sphere
//    scene.AddShape(mesh);
//    // Add lights
//    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.6f, 0.8f, 1.f), Spectrum(0.9f)));
//
//    // Create camera
//    auto camera = PinholeCamera(Vector3F(3.06f, 2.5f, 3.06f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f, WIDTH, HEIGHT);
//
//    // Create renderer
//    auto render = SimpleRenderer(std::make_shared<DirectIntegrator>());
//
//    // Render target image
//    BoxFilterFilm test(WIDTH, HEIGHT);
//    render.RenderImage(&test, scene, camera);
//
//    // Create tone-mapper and process target image
//    ClampTonemapper tonemapper;
//    tonemapper.Process("obj_test.png", test);

    // return EXIT_SUCCESS;
//
//    // Gradient
//    std::vector<float> gradient(mesh->GetNumVars(), 0.f);
//    std::vector<float> delta(mesh->GetNumVars(), 0.f);
//
//    // Number of iterations
//    size_t iters = 0;
//    // Energy value
//    float energy;
//
//    // Try to minimize squared norm of image
//    do {
//        // Clear derivatives
//        derivatives.Clear();
//        // Store current variables of the tape
//        default_tape.Push();
//
//        // Render current image
//        BoxFilterFilm x(WIDTH, HEIGHT);
//        render.RenderImage(&x, scene, camera);
//
//        // Output image of current rendering
//        tonemapper.Process("iters_" + std::to_string(iters) + ".ppm", x);
//
//        // Compute difference
//        // BoxFilterFilm difference = x - raw_target;
//
//        // Compute squared norm of difference
//        Float x_2_norm = x.Norm();
//        energy = x_2_norm.GetValue();
//        std::cout << "Energy: " << energy << std::endl;
//
//        // Create difference image
//        // difference.Abs();
//        // tonemapper.Process("iters_" + std::to_string(iters) + "_difference.ppm", difference);
//
//        // Compute derivatives
//        derivatives.ComputeDerivatives(x_2_norm);
//
//        // Get differentiable variables from scene's shapes, hardcoded TODO Fix this
//        std::vector<Float const *> vars;
//        scene.GetShapes()[0]->GetDiffVariables(vars);
//
//        // Compute gradient and deltas
//        for (size_t i = 0; i < vars.size(); i++) {
//            gradient[i] = derivatives.Dwrt(x_2_norm, *vars[i]);
//            delta[i] = -0.000001f * gradient[i];    // Learning rate
//        }
//
//        std::cout << "Iteration: " << iters << std::endl;
//        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//        // std::cout << "Gradient values: ";
//        // PrintGradient(gradient);
//
//        // Update scene vars, hardcoded for the moment
//        scene.GetShapes()[0]->UpdateDiffVariables(delta);
//
//        std::cout << "Tape size before pop: " << default_tape.Size() << std::endl;
//        // Pop variables
//        default_tape.Pop();
//
//        // Print tape size
//        std::cout << "Tape size after pop: " << default_tape.Size() << std::endl;
//        std::cout << std::endl;
//
//        // std::cout << "Sphere data" << std::endl;
//        // std::cout << scene.GetShapes()[0]->ToString() << std::endl << std::endl;
//        iters++;
//        // Stop when gradient is almost zero, "energy" is almost zero or maximum iterations reached
//    } while (GradNorm(gradient) > 0.001f && energy > 0.01f && iters < MAX_ITERS);
//
//    std::cout << "Final energy: " << energy << std::endl;
//    std::cout << "Total iterations: " << iters << std::endl;
//    std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//
//    // default_tape.Disable();
//    BoxFilterFilm final(WIDTH, HEIGHT);
//
//    render.RenderImage(&final, scene, camera);
//    tonemapper.Process("final.ppm", final);




    /**
     * Function minimisation test
     */
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
//
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


    /**
     * Geometric sphere description minimisation against target sphere image
     */

//    // Derivatives computation class
//    Derivatives derivatives;
//
//    // Disable tape to render target image
//    // default_tape.Disable();
//
//    // Create scene
//    Scene scene;
//
//    // Add sphere
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(1.f, 0.f, 0.f), Float(2.f)));
//    // Add lights
//    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));
//
//    // Create camera
//    auto camera = PinholeCamera(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f, WIDTH, HEIGHT);
//
//    // Create renderer
//    auto render = SimpleRenderer(std::make_shared<DirectIntegrator>());
//
//    // Push status of tape before rendering target image
//    default_tape.Push();
//
//    // Render target image
//    BoxFilterFilm target(WIDTH, HEIGHT);
//    render.RenderImage(&target, scene, camera);
//
//    // Convert image to raw
//    std::vector<float> raw_target = target.Raw();
//
//    // Create tone-mapper and process target image
//    ClampTonemapper tonemapper;
//    tonemapper.Process("target.png", target);
//
//    // Re-enable tape to compute derivatives
//    // default_tape.Enable();
//
//    // Remove from tape rendering variables
//    default_tape.Pop();
//
//    // Change sphere position to center and try to match the images
//    scene.ClearShapes();
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(1.f)));
//
//    // Gradient
//    std::vector<float> gradient(4, 0.f);
//    std::vector<float> delta(4, 0.f);
//
//    // Get differentiable variables from scene's shapes
//    std::vector<Float const *> vars;
//    scene.GetShapes()[0]->GetDiffVariables(vars);
//
//    // Number of iterations
//    size_t iters = 0;
//    // Energy value
//    float energy;
//
//    // Try to minimize squared norm of image
//    do {
//        // Clear derivatives
//        derivatives.Clear();
//        // Store current variables of the tape
//        default_tape.Push();
//
//        // Render current image
//        BoxFilterFilm x(WIDTH, HEIGHT);
//        render.RenderImage(&x, scene, camera);
//
//        // Output image of current rendering
//        tonemapper.Process("iters_" + std::to_string(iters) + ".png", x);
//
//        // Compute difference
//        BoxFilterFilm difference = x - raw_target;
//
//        // Compute squared norm of difference
//        Float x_2_norm = difference.Norm();
//        energy = x_2_norm.GetValue();
//        std::cout << "Energy: " << energy << std::endl;
//
//        // Create difference image
//        difference.Abs();
//        tonemapper.Process("iters_" + std::to_string(iters) + "_difference.ppm", difference);
//
//        // Compute derivatives
//        derivatives.ComputeDerivatives(x_2_norm);
//
//        // Compute gradient and deltas
//        for (size_t i = 0; i < vars.size(); i++) {
//            gradient[i] = derivatives.Dwrt(x_2_norm, *vars[i]);
//            delta[i] = -0.000001f * gradient[i];    // Learning rate
//        }
//
//        std::cout << "Iteration: " << iters << std::endl;
//        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//        std::cout << "Gradient values: ";
//        PrintGradient(gradient);
//
//        // Update scene vars, hardcoded for the moment
//        scene.GetShapes()[0]->UpdateDiffVariables(delta, 0);
//
//        std::cout << "Tape size before pop: " << default_tape.Size() << std::endl;
//        // Pop variables
//        default_tape.Pop();
//
//        // Print tape size
//        std::cout << "Tape size after pop: " << default_tape.Size() << std::endl;
//
//        std::cout << "Sphere data" << std::endl;
//        std::cout << scene.GetShapes()[0]->ToString() << std::endl << std::endl;
//        iters++;
//        // Stop when gradient is almost zero, "energy" is almost zero or maximum iterations reached
//    } while (GradNorm(gradient) > 0.001f && energy > 0.01f && iters < MAX_ITERS);
//
//    std::cout << "Final energy: " << energy << std::endl;
//    std::cout << "Total iterations: " << iters << std::endl;
//    std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//
//    // default_tape.Disable();
//    BoxFilterFilm final(WIDTH, HEIGHT);
//
//    render.RenderImage(&final, scene, camera);
//    tonemapper.Process("final.png", final);
//
//    std::cout << "Final sphere data" << std::endl;
//    std::cout << scene.GetShapes()[0]->ToString() << std::endl;

    /*
     * Signed distance grid rendering test
     */
//
//    // Create new grid
//    int grid_dims[3] = {10, 10, 10};
//    auto grid = std::make_shared<SignedDistanceGrid>(grid_dims[0], grid_dims[1], grid_dims[2],
//                                                     BBOX(Vector3f(-2.f, -2.f, -2.f), Vector3f(2.f, 2.f, 2.f)));
//
//    float delta = 4.f / static_cast<float>(grid_dims[0] - 1);
//
//    // Initialize grid using sphere of radius 1 as SDF
//    for (int z = 0; z < grid_dims[2]; z++) {
//        for (int y = 0; y < grid_dims[1]; y++) {
//            for (int x = 0; x < grid_dims[0]; x++) {
//                // Compute point coordinates
//                Vector3f p(-2.f + delta * x, -2.f + delta * y, -2.f + delta * z);
//                // Use ellipse equation
//                grid->operator()(x, y, z) = Length(p) - 1.f;
//            }
//        }
//    }
//
//    int new_dims[3] = {30, 30, 30};
//    grid->Refine(new_dims);


//    // Create scene
//    Scene scene;
////
////    // Add grid
//    scene.AddShape(grid);
////
////    // Add sphere
////    // scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(1.2f)));
////    // Add lights
//    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));
//
//    // Create target_cameras
//    std::vector<std::shared_ptr<const CameraInterface> > cameras;
//
//////    // Add first camera
//////    cameras.push_back(
//////            std::make_shared<const PinholeCamera>(Vector3F(5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f,
//////                                                  WIDTH, HEIGHT));
//////    // Add second camera
//////    cameras.push_back(
//////            std::make_shared<const PinholeCamera>(Vector3F(-5.f, 0.f, 0.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f,
//////                                                  WIDTH, HEIGHT));
//
//    cameras.push_back(
//            std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 5.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f,
//                                                  WIDTH, HEIGHT));
//
//    // Create renderer
//    auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());
//
//////    // Push status of tape before rendering target image
//////    default_tape.Push();
////
//    // Render target images
//    BoxFilterFilm target(WIDTH, HEIGHT);
//
//    // Create target image
//    render->RenderImage(&target, scene, *cameras[0]);
//    ClampTonemapper tonemapper;
//    tonemapper.Process("render_test.png", target);
//
//    // Refine and render again
//    int new_dims[3] = {30, 30, 30};
//    grid->Refine(new_dims);
//    render->RenderImage(&target, scene, *cameras[0]);
//    tonemapper.Process("render_test_refine.png", target);

//    std::vector<std::vector<float> > raw_views;
//    // Render views
//    for (auto const &camera : cameras) {
//        render->RenderImage(&target, scene, *camera);
//        raw_views.push_back(target.Raw());
//    }
//
//    // Use SDF
//    scene.ClearShapes();
//    scene.AddShape(grid);
//
//    render->RenderImage(&target, scene, *cameras[2]);
//    tonemapper.Process("start.png", target);
//
//    // Remove from tape rendering variables
//    default_tape.Pop();
//
//    // Create multi-view energy
//    // auto energy = MultiViewEnergy(scene, raw_views, cameras, render, WIDTH, HEIGHT);
//
//    // Test with new energy
//    auto energy = ReconstructionEnergy(scene, grid, raw_views, cameras, render, 1.f, WIDTH, HEIGHT);
//
//    // Minimise energy
//    // GradientDescent::Minimize(energy, 0.0002f, MAX_ITERS, 1.f, true);
//    GradientDescentBT::Minimize(energy, MAX_ITERS, 50.f, 0.5f, 0.8f, true, true, 100.f, 50.f);
//
//    render->RenderImage(&target, scene, *cameras[2]);
//    tonemapper.Process("final.png", target);


//    /**
//     * SDF sphere description minimisation against target sphere image
//     */
//
//    // Derivatives computation class
//    Derivatives derivatives;
//
//    // Create scene
//    Scene scene;
//
//    // Add sphere
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(1.1f)));
//    // Add lights
//    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));
//
//    // Create camera
//    auto camera = PinholeCamera(Vector3F(0.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f, WIDTH, HEIGHT);
//
//    // Create renderer
//    auto render = SimpleRenderer(std::make_shared<DirectIntegrator>());
//
//    // Push status of tape before rendering target image
//    default_tape.Push();
//
//    // Render target image
//    BoxFilterFilm target(WIDTH, HEIGHT);
//    render.RenderImage(&target, scene, camera);
//
//    // Convert image to raw
//    std::vector<float> raw_target = target.Raw();
//
//    // Create tone-mapper and process target image
//    ClampTonemapper tonemapper;
//    tonemapper.Process("target.png", target);
//
//    // Remove from tape rendering variables
//    default_tape.Pop();
//
//    // Change sphere position to center and try to match the images
//    scene.ClearShapes();
//
//    // Create new grid
//    int grid_dims[3] = {100, 100, 100};
//    auto grid = std::make_shared<SignedDistanceGrid>(grid_dims[0], grid_dims[1], grid_dims[2],
//                                                     BBOX(Vector3f(2.f, 2.f, 2.f), Vector3f(-2.f, -2.f, -2.f)));
//
//    // FIXME Hardcoded for testing
//    float delta_s = 4.f / 99.f;
//    // Initialize grid using sphere of radius 1 as SDF
//    for (int x = 0; x < grid_dims[0]; x++) {
//        for (int y = 0; y < grid_dims[1]; y++) {
//            for (int z = 0; z < grid_dims[2]; z++) {
//                // Compute point coordinates
//                Vector3f p(-2.f + delta_s * x, -2.f + delta_s * y, -2.f + delta_s * z);
//                grid->operator()(x, y, z) = Length(p) - 1.f;
//            }
//        }
//    }
//
//    // Print starting grid values
////    std::cout << "Starting grid data: " << std::endl;
////    std::cout << grid->ToString() << std::endl;
//
//    scene.AddShape(grid);
//
//    // Gradient
//    std::vector<float> gradient(grid->GetNumVars(), 0.f);
//    std::vector<float> delta(grid->GetNumVars(), 0.f);
//
//    // Get differentiable variables from scene's shapes
//    std::vector<Float const *> vars;
//    scene.GetShapes()[0]->GetDiffVariables(vars);
//
//    // Number of iterations
//    size_t iters = 0;
//    // Energy value
//    float energy;
//
//    // Try to minimize squared norm of image
//    do {
//        // Clear derivatives
//        derivatives.Clear();
//        // Store current variables of the tape
//        default_tape.Push();
//
//        // Render current image
//        BoxFilterFilm x(WIDTH, HEIGHT);
//        render.RenderImage(&x, scene, camera);
//
//        // Output image of current rendering
//        tonemapper.Process("iters_" + std::to_string(iters) + ".png", x);
//
//        // Compute difference
//        BoxFilterFilm difference = x - raw_target;
//
//        // Compute squared norm of difference
//        Float x_2_norm = difference.Norm();
//        energy = x_2_norm.GetValue();
//        std::cout << "Energy: " << energy << std::endl;
//
//        // Create difference image
//        difference.Abs();
//        tonemapper.Process("iters_" + std::to_string(iters) + "_difference.ppm", difference);
//
//        // Compute derivatives
//        derivatives.ComputeDerivatives(x_2_norm);
//
//        // Compute gradient and deltas
//        for (size_t i = 0; i < vars.size(); i++) {
//            gradient[i] = derivatives.Dwrt(x_2_norm, *vars[i]);
//            delta[i] = -0.0001f * gradient[i];    // Learning rate
//        }
//
//        std::cout << "Iteration: " << iters << std::endl;
//        std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//        // std::cout << "Gradient values: " << std::endl;
//        // PrintGradient(gradient);
//
//        // Update scene vars, hardcoded for the moment
//        scene.GetShapes()[0]->UpdateDiffVariables(delta, 0);
//
//        // Grid data
//        // std::cout << "Grid data: " << std::endl;
//        // std::cout << grid->ToString() << std::endl;
//
//        std::cout << "Tape size before pop: " << default_tape.Size() << std::endl;
//        // Pop variables
//        default_tape.Pop();
//
//        // Print tape size
//        std::cout << "Tape size after pop: " << default_tape.Size() << std::endl << std::endl;
//
//        iters++;
//        // Stop when gradient is almost zero, "energy" is almost zero or maximum iterations reached
//    } while (GradNorm(gradient) > 0.001f && energy > 0.01f && iters < MAX_ITERS);
//
//    std::cout << "Final energy: " << energy << std::endl;
//    std::cout << "Total iterations: " << iters << std::endl;
//    std::cout << "Gradient norm: " << GradNorm(gradient) << std::endl;
//
//    // default_tape.Disable();
//    BoxFilterFilm final(WIDTH, HEIGHT);
//
//    render.RenderImage(&final, scene, camera);
//    tonemapper.Process("final.png", final);



//    /**
//     * Geometric sphere description multi-view test
//     */
//
//    // Create scene
//    Scene scene;
//
//    // Add sphere
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(1.f, 0.f, 0.f), Float(2.f)));
//    // Add lights
//    scene.AddLight(std::make_shared<DirectionalLight>(Vector3F(0.f, 0.f, 1.f), Spectrum(0.9f)));
//
//    // Create target_cameras
//    std::vector<std::shared_ptr<const CameraInterface> > target_cameras;
//    // Add first camera
//    target_cameras.push_back(
//            std::make_shared<const PinholeCamera>(Vector3F(1.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f,
//                                                  WIDTH, HEIGHT));
//    // Add second camera
//    target_cameras.push_back(
//            std::make_shared<const PinholeCamera>(Vector3F(-1.f, 0.f, 10.f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60.f,
//                                                  WIDTH, HEIGHT));
//
//    // Create renderer
//    auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());
//
//    // Push status of tape before rendering target image
//    default_tape.Push();
//
//    // Render target images
//    BoxFilterFilm target(WIDTH, HEIGHT);
//    std::vector<std::vector<float> > raw_views;
//    // Render views
//    for (auto const &camera : target_cameras) {
//        render->RenderImage(&target, scene, *camera);
//        raw_views.push_back(target.Raw());
//    }
//
//    // Remove from tape rendering variables
//    default_tape.Pop();
//
//    // Change sphere position to center, reduce radius and try to match the images
//    scene.ClearShapes();
//    scene.AddShape(std::make_shared<Sphere>(Vector3F(0.f, 0.f, 0.f), Float(1.f)));
//
//    // Create multi-view energy
//    auto energy = MultiViewEnergy(scene, raw_views, target_cameras, render, WIDTH, HEIGHT);
//
//    // Minimise energy
//    // GradientDescent::Minimize(energy, 0.000001f, MAX_ITERS, 1.f, true);
//    GradientDescentBT::Minimize(energy, MAX_ITERS, 1.f, 0.5f, 0.8f, true, true, 0.f, 8.f);
//
//    std::cout << "Final sphere data" << std::endl;
//    std::cout << scene.GetShapes()[0]->ToString() << std::endl;

    return EXIT_SUCCESS;
}
