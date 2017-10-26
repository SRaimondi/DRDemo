//
// Created by Simon on 02.10.2017.
//

#include <cstdlib>
#include <box_film.hpp>
#include <iofile.hpp>
#include <memory>
#include <camera.hpp>
#include <perspective_camera.hpp>
#include <grid.hpp>
#include <scene.hpp>
#include <direct_integrator.hpp>
#include <simple_renderer.hpp>
#include <clamp_tonemapper.hpp>
#include <reconstruction_energy_opt.hpp>
#include <gradient_descent.hpp>
#include <triangle_mesh.hpp>
#include <mac_grid.hpp>
#include <ambient_light.hpp>
#include <reconstruction_energy_light.hpp>
#include "dino_test.hpp"

namespace drdemo {

    void DinoTest(float res_multiplier, int ref_steps) {
        // Maximum gradient descent iterations
        const size_t MAX_ITERS = 250;

        // Disable tape, load target images
        default_tape.Disable();
        std::vector<std::vector<float> > raw_views;
        // Target images width and height
        size_t width, height;
        for (int i = 1; i <= 16; i += 2) {
            // Create file names
            std::string file_name("../dinoSparseRing/dinoSR00");
            if (i < 10) {
                file_name += "0";
            }
            file_name += std::to_string(i) + ".png";
            // Load file
            BoxFilterFilm loaded_image = BoxFilterFilm::FromPNG(file_name);
            // Initialize image size once
            if (i == 1) {
                width = loaded_image.Width();
                height = loaded_image.Height();
            }
            raw_views.push_back(loaded_image.Raw());
        }

        // Enable tape
        default_tape.Enable();

        // Load cameras
        std::vector<std::shared_ptr<const CameraInterface> > cameras;
        std::vector<std::string> p_inv_file;
        // if (ReadFile("../dinoSparseRing/dinoSR_pseudo_inv.txt", p_inv_file)) {
        if (ReadFile("../dinoSparseRing/dinoSR_inv_mat.txt", p_inv_file)) {
            // float p_inv[12];
            float m_inv[9];
            float c_w[3];
            for (int i = 0; i < 8; ++i) {
//                if (sscanf(p_inv_file[2 * i].c_str(),
//                           "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
//                           &p_inv[0], &p_inv[1], &p_inv[2],
//                           &p_inv[3], &p_inv[4], &p_inv[5],
//                           &p_inv[6], &p_inv[7], &p_inv[8],
//                           &p_inv[9], &p_inv[10], &p_inv[11],
//                           &c_w[0], &c_w[1], &c_w[2]) == 15) {
//                    cameras.push_back(std::make_shared<const PerspectiveCamera>(p_inv, c_w, width, height));
//                } else {
//                    std::cerr << "Error parsing camera parameters!" << std::endl;
//                    exit(EXIT_FAILURE);
//                }
                if (sscanf(p_inv_file[2 * i].c_str(),
                           "%f %f %f %f %f %f %f %f %f %f %f %f",
                           &m_inv[0], &m_inv[1], &m_inv[2],
                           &m_inv[3], &m_inv[4], &m_inv[5],
                           &m_inv[6], &m_inv[7], &m_inv[8],
                           &c_w[0], &c_w[1], &c_w[2]) == 12) {
                    cameras.push_back(std::make_shared<const PerspectiveCamera>(m_inv, c_w, width, height));
                } else {
                    std::cerr << "Error parsing camera parameters!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }

        // cameras.push_back(std::make_shared<const PinholeCamera>(Vector3F(0.f, 0.f, 0.2f), Vector3F(), Vector3F(0.f, 1.f, 0.f), 60, width, height));

        /**
         * The (tight) bounding box for the dino model is
         * (-0.061897 -0.018874 -0.057845)
         * (0.010897 0.068227 0.015495)
         */
//        BBOX sdf_bbox(Vector3f(-0.1f, -0.1f, -0.1f), Vector3f(0.1f, 0.1f, 0.1f));
//        int grid_dims[3] = {start_resolution, start_resolution, start_resolution};
//        auto grid = std::make_shared<SignedDistanceGrid>(grid_dims[0], grid_dims[1], grid_dims[2], sdf_bbox);
//
//        // Initialize grid using sphere of radius 0.05f
//        for (int z = 0; z < grid_dims[2]; z++) {
//            for (int y = 0; y < grid_dims[1]; y++) {
//                for (int x = 0; x < grid_dims[0]; x++) {
//                    // Compute point coordinates
//                    const Vector3f p = grid->CoordsAt(x, y, z);
//                    // Set SDF value
//                    grid->operator()(x, y, z) = Length(p) - 0.05f;
//                }
//            }
//        }

        // Load starting SDF from file
        auto grid = std::make_shared<SignedDistanceGrid>("../sdfs/dino_watertight.sdf");
        // auto grid = std::make_shared<SignedDistanceGrid>("../sdfs/dino_watertight_low.sdf");
        // auto new_grid = std::make_shared<MACGrid>("../sdfs/dino_watertight.sdf");

        // Load .obj of dino
        // auto mesh = std::make_shared<TriangleMesh>("../objs/dino_watertight_n.obj");

        // Create scene and add grid
        Scene scene;
        // scene.AddShape(new_grid);
        scene.AddShape(grid);
        // scene.AddShape(mesh);

        // Create and add ambient light
        auto light = std::make_shared<AmbientLight>(Spectrum(1.f, 1.f, 1.f));
        scene.AddLight(light);

        // Create renderer class with direct illumination integrator
        auto render = std::make_shared<SimpleRenderer>(std::make_shared<DirectIntegrator>());

        // Render start images
        default_tape.Disable();

        BoxFilterFilm target(width, height);
        ClampTonemapper tonemapper;

        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("start_" + std::to_string(i) + ".png", target);
        }

        // Enable tape again
        default_tape.Enable();

        // Create energy
        // auto energy = ReconstructionEnergyOpt(scene, grid, raw_views, cameras, render, 1.f, width, height);
        auto energy = ReconstructionEnergyLight(scene, grid, light, raw_views, cameras, render, 1.f, width, height);

        // Do first minimisation
        GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-15f, true);

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
            GradientDescentBT::Minimize(energy, MAX_ITERS, 10.f, 0.5f, 0.8f, 10e-15f, true);
        }

        default_tape.Disable();
        // Render final SDF status
        for (int i = 0; i < cameras.size(); i++) {
            render->RenderImage(&target, scene, *cameras[i]);
            tonemapper.Process("final_ " + std::to_string(i) + ".png", target);
        }
        default_tape.Enable();
    }

} // drdemo namespace