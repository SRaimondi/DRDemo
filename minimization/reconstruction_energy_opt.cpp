//
// Created by Simon on 08.09.2017.
//

#include "reconstruction_energy_opt.hpp"
#include "box_film.hpp"

namespace drdemo {

    ReconstructionEnergyOpt::ReconstructionEnergyOpt(Scene &scene,
                                                     const std::shared_ptr<SignedDistanceGrid> &grid,
                                                     const std::vector<std::vector<float> > &views,
                                                     const std::vector<std::shared_ptr<const CameraInterface> > &c,
                                                     const std::shared_ptr<RendererInterface> &r,
                                                     float lambda,
                                                     size_t w, size_t h)
            : target_scene(scene), grid(grid), target_views(views), target_cameras(c),
              renderer(r), gradient(grid->GetNumVars(), 0.f), lambda(lambda), width(w), height(h),
              evaluations(0) {
        // Check that the size of the target render and the target_cameras is the same
        assert(target_views.size() == target_cameras.size());

        // Get pointer to all the differentiable variables of the grid
        this->grid->GetDiffVariables(diff_variables);
    }

    void ReconstructionEnergyOpt::RebindVars() {
        // Clear variables we need to compute the derivative with respect to
        diff_variables.clear();
        // Rebind
        grid->GetDiffVariables(diff_variables);
        // Change gradient size according to new number of variables
        gradient.resize(diff_variables.size());
    }

    size_t ReconstructionEnergyOpt::InputDim() const {
        return grid->GetNumVars();
    }

    /**
     * The energy here is computed as the average of the difference between the rendered images and the target,
     * plus the sum of the squared norm of each normal minus one, since we want them to stay as close as one as possible
     *
     * @return Final energy
     */
    Float ReconstructionEnergyOpt::Evaluate(bool output) const {
        // First energy term that contains the sum of the difference between the rendered images and the targets
        Float E_images;
        // Image energy term computed currently
        Float E_image_t;

        // bool gradient_to_zero = false;
        // Clear derivatives
        derivatives.Clear();

        // Film to render the image on
        BoxFilterFilm render(width, height);
        BoxFilterFilm difference(width, height);

        // Current evaluated energy term gradient
        // std::vector<float> image_term_grad(gradient.size(), 0.f);

        // Loop over all target target_cameras
        for (size_t target_index = 0; target_index < target_cameras.size(); ++target_index) {
            // Push where we are before rendering current image
            default_tape.Push();

            // Render scene for current camera
            renderer->RenderImage(&render, target_scene, *target_cameras[target_index]);

            // If we are at view zero and output is true, output image
            if (output && target_index == 0) {
                tonemapper.Process("iterations_" + std::to_string(evaluations) + ".png", render);
            }

            // Compute difference between rendering and target
            difference = render - target_views[target_index];
            // Compute single image energy
            E_image_t = difference.SquaredNorm();

            // Check if we need to compute the gradient
            if (default_tape.IsEnabled()) {
                // Compute gradient
                if (target_index == 0) {
                    // Reset gradient
                    for (auto &v : gradient) { v = 0.f; }
                }
                // Compute derivatives for current image term
                derivatives.Clear();
                derivatives.ComputeDerivatives(E_image_t);
                // Compute gradient for current term
                for (size_t i = 0; i < gradient.size(); ++i) {
                    gradient[i] += derivatives.Dwrt(E_image_t, *diff_variables[i]);
                }
                // Add contribution to final gradient
//                for (size_t i = 0; i < gradient.size(); ++i) {
//                    gradient[i] += image_term_grad[i];
//                }
            }

            // Sum current rendering difference to total energy
            E_images += E_image_t;      // FIXME We can not use this value anymore to compute derivateives!!!
            default_tape.Pop();
        }

        // Second energy term that contains the sum of the squared norms of the normals minus 1 (each one)
        Float E_normals;

        default_tape.Push();
        // Loop over all internal points of the grid, compute the normal and sum up
        for (int z = 0; z < grid->Size(2); z++) {
            for (int y = 0; y < grid->Size(1); y++) {
                for (int x = 0; x < grid->Size(0); x++) {
                    // Get normal
                    const Vector3F n = grid->NormalAtPoint(x, y, z);
                    // Add to energy
                    E_normals += Pow(LengthSquared(n) - 1.f, 2.f);
                }
            }
        }
        if (default_tape.IsEnabled()) {
            // Add final contribution of normal term to gradient
            derivatives.Clear();
            derivatives.ComputeDerivatives(E_normals);
            // Compute gradient for current term
//            for (size_t i = 0; i < gradient.size(); ++i) {
//                image_term_grad[i] = derivatives.Dwrt(E_normals, *diff_variables[i]);
//            }
            // Add contribution to final gradient
            for (size_t i = 0; i < gradient.size(); ++i) {
                // Lambda * gradient !!!
                gradient[i] += lambda * derivatives.Dwrt(E_normals, *diff_variables[i]);
            }
        }
        // We have our gradient, clear derivatives
        // derivatives.Clear();
        // Pop tape
        default_tape.Pop();

        // Increase number of evaluations if we used the ouput
        if (output) { evaluations++; }
        // Store last computed values for the terms
        image_term = E_images.GetValue();
        normal_term = E_normals.GetValue();

        // Return final energy = E_images + lambda * E_normals
        return E_images + lambda * E_normals;
    }

    std::vector<float> ReconstructionEnergyOpt::ComputeGradient(const Float &) const {
        // Return copy of current gradient
        return gradient;
    }

    void ReconstructionEnergyOpt::UpdateStatus(const std::vector<float> &deltas) {
        // Here we assume the only thing to be updates is the grid
        grid->UpdateDiffVariables(deltas, 0);
    }

    void ReconstructionEnergyOpt::SetStatus(const std::vector<float> &new_status) {
        // Only set the status of the grid
        grid->SetDiffVariables(new_status, 0);
    }

    std::vector<float> ReconstructionEnergyOpt::GetStatus() const {
        // Status is just the current value of all the grid values we can differentiate with respect to
        std::vector<float> status(diff_variables.size(), 0.f);
        for (size_t i = 0; i < status.size(); i++) {
            status[i] = diff_variables[i]->GetValue();
        }

        return status;
    }

    std::string ReconstructionEnergyOpt::ToString() const {
        return "Image term: " + std::to_string(image_term) + "\n" + "Normal term: " + std::to_string(normal_term);
    }

} // drdemo namespace
