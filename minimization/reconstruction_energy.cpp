//
// Created by Simon on 21.08.2017.
//

#include "reconstruction_energy.hpp"
#include "box_film.hpp"

namespace drdemo {

    ReconstructionEnergy::ReconstructionEnergy(Scene &scene,
                                               const std::shared_ptr<SignedDistanceGrid> &grid,
                                               const std::vector<std::vector<float> > &views,
                                               const std::vector<std::shared_ptr<const CameraInterface> > &c,
                                               const std::shared_ptr<RendererInterface> &r,
                                               float lambda,
                                               size_t w, size_t h)
            : target_scene(scene), grid(grid), target_views(views), target_cameras(c),
              renderer(r), lambda(lambda), width(w), height(h), evalutations(0) {
        // Check that the size of the target render and the target_cameras is the same
        assert(target_views.size() == target_cameras.size());

        // Get pointer to all the differentiable variables of the grid
        this->grid->GetDiffVariables(diff_variables);
    }

    void ReconstructionEnergy::RebindVars() {
        // Clear variables we need to compute the derivative with respect to
        diff_variables.clear();
        // Rebind
        grid->GetDiffVariables(diff_variables);
    }

    size_t ReconstructionEnergy::InputDim() const {
        return grid->GetNumVars();
    }

    /**
     * The energy here is computed as the average of the differnce between the rendered images and the target,
     * plus the sum of the squared norm of each normal minus one, since we want them to stay as close as one as possible
     *
     * @return Final energy
     */
    Float ReconstructionEnergy::Evaluate(bool output) const {
        // First energy term that contains the sum of the difference between the rendered images and the targets
        Float E_images;

        // Loop over all target target_cameras
        for (size_t target_index = 0; target_index < target_cameras.size(); ++target_index) {
            // Film to render the image on
            BoxFilterFilm render(width, height);
            // Render scene for current camera
            renderer->RenderImage(&render, target_scene, *target_cameras[target_index]);

            // If we are at view zero and output is true, output image
            if (output && target_index == 0) {
                tonemapper.Process("iterations_" + std::to_string(evalutations) + ".png", render);
            }

            // Compute difference between rendering and target
            BoxFilterFilm difference = render - target_views[target_index];
            // Compute energy of difference
            Float difference_norm = difference.SquaredNorm();

            // Sum current rendering difference to total energy
            E_images += difference_norm;
        }

        // Second energy term that contains the sum of the squared norms of the normals minus 1 (each one)
        Float E_normals;

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

        // Increase number of evaluations if we used the ouput
        if (output) { evalutations++; }
        // Store last computed values for the terms
        image_term = E_images.GetValue();
        normal_term = E_normals.GetValue();


//        std::cout << "Images energy: " << E_images << std::endl;
//        std::cout << "Normal energy: " << lambda * E_normals << std::endl;

        // Return final energy = E_images + lambda * E_normals
        return E_images + lambda * E_normals;
    }

    std::vector<float> ReconstructionEnergy::ComputeGradient(const Float &out) const {
        // Clear current derivatives
        derivatives.Clear();
        // Compute derivatives for out variable
        derivatives.ComputeDerivatives(out);

//        // Check if we need to rebind the differentiable variables
//        if (grid->GetNumVars() != input_dim) {
//            diff_variables.clear();
//            grid->GetDiffVariables(diff_variables);
//            input_dim = grid->GetNumVars();
//        }

        // Create and compute gradient
        std::vector<float> gradient(diff_variables.size(), 0.f);
        for (size_t v = 0; v < diff_variables.size(); ++v) {
            gradient[v] = derivatives.Dwrt(out, *diff_variables[v]);
        }

        // Clear derivatives
        derivatives.Clear();

        return gradient;
    }

    void ReconstructionEnergy::UpdateStatus(const std::vector<float> &deltas) {
        // Here we assume the only thing to be updates is the grid
        grid->UpdateDiffVariables(deltas, 0);
    }

    void ReconstructionEnergy::SetStatus(const std::vector<float> &new_status) {
        // Only set the status of the grid
        grid->SetDiffVariables(new_status, 0);
    }

    std::vector<float> ReconstructionEnergy::GetStatus() const {
        // Status is just the current value of all the gird values we can differentiate with respect to
        std::vector<float> status(diff_variables.size(), 0.f);
        for (size_t i = 0; i < status.size(); i++) {
            status[i] = diff_variables[i]->GetValue();
        }

        return status;
    }

    std::string ReconstructionEnergy::ToString() const {
        return "Image term: " + std::to_string(image_term) + "\n" + "Normal term: " + std::to_string(normal_term);
    }

} // drdemo namespace
