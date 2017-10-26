//
// Created by simon on 25.10.17.
//

#include "reconstruction_energy_light.hpp"
#include "box_film.hpp"

namespace drdemo {

    ReconstructionEnergyLight::ReconstructionEnergyLight(Scene &scene,
                                                         const std::shared_ptr<SignedDistanceGrid> &grid,
                                                         const std::shared_ptr<AmbientLight> &light,
                                                         const std::vector<std::vector<float> > &views,
                                                         const std::vector<std::shared_ptr<const CameraInterface> > &c,
                                                         const std::shared_ptr<RendererInterface> &r,
                                                         float lambda,
                                                         size_t w, size_t h)
            : target_scene(scene), grid(grid), light(light), target_views(views), target_cameras(c),
              renderer(r), gradient(grid->GetNumVars() + light->GetNumVars(), 0.f), lambda(lambda), width(w), height(h),
              evaluations(0) {
        // Check that the size of the target render and the target_cameras is the same
        assert(target_views.size() == target_cameras.size());

        // Get pointer to all the differentiable variables of the grid
        // The gradient stores first the SDF values and after the light parameters
        this->grid->GetDiffVariables(diff_variables);
        this->light->GetDiffVariables(diff_variables);
    }

    void ReconstructionEnergyLight::RebindVars() {
        // Clear variables we need to compute the derivative with respect to
        diff_variables.clear();
        // Rebind
        grid->GetDiffVariables(diff_variables);
        light->GetDiffVariables(diff_variables);
        // Change gradient size according to new number of variables
        gradient.resize(diff_variables.size());
    }

    size_t ReconstructionEnergyLight::InputDim() const {
        return grid->GetNumVars() + light->GetNumVars();
    }

    /**
     * The energy here is computed as the sum of the difference between the rendered images and the target,
     * plus the sum of the squared norm of each normal minus one, since we want them to stay as close as one as possible
     *
     * @return Final energy
     */
    Float ReconstructionEnergyLight::Evaluate(bool output) const {
        // First energy term that contains the sum of the difference between the rendered images and the targets
        Float E_images;
        // Image energy term computed currently
        Float E_image_t;

        // Clear derivatives
        derivatives.Clear();

        // Film to render the image on
        BoxFilterFilm render(width, height);
        BoxFilterFilm difference(width, height);

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

            // Add contribution to final gradient (-3 because the contribution of the light intensity is null here)
            for (size_t i = 0; i < gradient.size() - 3; ++i) {
                // Lambda * gradient !!!
                gradient[i] += lambda * derivatives.Dwrt(E_normals, *diff_variables[i]);
            }
        }
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

    std::vector<float> ReconstructionEnergyLight::ComputeGradient(const Float &) const {
        // Return copy of current gradient
        return gradient;
    }

    void ReconstructionEnergyLight::UpdateStatus(const std::vector<float> &deltas) {
        // Here we assume the only thing to be updates is the grid
        grid->UpdateDiffVariables(deltas, 0);
        light->UpdateDiffVariables(deltas, grid->GetNumVars());
    }

    void ReconstructionEnergyLight::SetStatus(const std::vector<float> &new_status) {
        // Only set the status of the grid
        grid->SetDiffVariables(new_status, 0);
        light->SetDiffVariables(new_status, grid->GetNumVars());
    }

    std::vector<float> ReconstructionEnergyLight::GetStatus() const {
        // Status is just the current value of all the grid values we can differentiate with respect to
        std::vector<float> status(diff_variables.size(), 0.f);
        for (size_t i = 0; i < status.size(); i++) {
            status[i] = diff_variables[i]->GetValue();
        }

        return status;
    }

    std::string ReconstructionEnergyLight::ToString() const {
        return "Image term: " + std::to_string(image_term) + "\n" +
               "Normal term: " + std::to_string(normal_term) + "\n" +
               "Light color: [" + std::to_string(light->intensity.r.GetValue()) + ", " +
               std::to_string(light->intensity.g.GetValue()) + ", " +
               std::to_string(light->intensity.b.GetValue()) + "]";
    }

} // drdemo namespace
