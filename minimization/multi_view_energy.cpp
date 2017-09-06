//
// Created by Simone Raimondi on 08.08.17.
//

#include "box_film.hpp"
#include "multi_view_energy.hpp"

namespace drdemo {

    MultiViewEnergy::MultiViewEnergy(Scene &scene, const std::vector<std::vector<float> > &views,
                                     const std::vector<std::shared_ptr<const CameraInterface> > &c,
                                     const std::shared_ptr<RendererInterface> &r,
                                     size_t w, size_t h)
            : target_scene(scene), target_views(views), cameras(c), renderer(r), width(w), height(h), evalutations(0) {
        // Check that the size of the target render and the target_cameras is the same
        assert(target_views.size() == cameras.size());

        // Get pointer to all the differentiable variables of the Shapes in the scene
        for (auto const &shape : target_scene.GetShapes()) {
            shape->GetDiffVariables(diff_variables);
        }
    }

    Float drdemo::MultiViewEnergy::Evaluate(bool output) const {
        // Final energy, sum of all the single rendering differences
        Float energy;

        // Loop over all target target_cameras
        for (size_t target_index = 0; target_index < cameras.size(); ++target_index) {
            // Film to render the image on
            BoxFilterFilm render(width, height);
            // Render scene for current camera
            renderer->RenderImage(&render, target_scene, *cameras[target_index]);
            // Compute difference between rendering and target
            BoxFilterFilm difference = render - target_views[target_index];
            // Compute energy of difference
            Float difference_norm = difference.SquaredNorm();

            // Sum current rendering difference to total energy
            energy += difference_norm;
        }

        // Increase number of evaluations
        evalutations++;

        return energy;
    }

    std::vector<float> MultiViewEnergy::ComputeGradient(const Float &out) const {
        // Clear current derivatives
        derivatives.Clear();
        // Compute derivatives for out variable
        derivatives.ComputeDerivatives(out);

        // Create and compute gradient
        std::vector<float> gradient(diff_variables.size(), 0.f);
        for (size_t v = 0; v < diff_variables.size(); ++v) {
            gradient[v] = derivatives.Dwrt(out, *diff_variables[v]);
        }

        return gradient;
    }

    void MultiViewEnergy::UpdateStatus(const std::vector<float> &deltas) {
        size_t starting_index = 0;
        for (auto &shape : target_scene.GetShapes()) {
            // Update variables
            shape->UpdateDiffVariables(deltas, starting_index);
            // Increase starting index by number of used variables
            starting_index += shape->GetNumVars();
        }
    }

    void MultiViewEnergy::SetStatus(const std::vector<float> &new_status) {
        size_t starting_index = 0;
        for (auto &shape : target_scene.GetShapes()) {
            // Update variables
            shape->SetDiffVariables(new_status, starting_index);
            // Increase starting index by number of used variables
            starting_index += shape->GetNumVars();
        }
    }

    size_t MultiViewEnergy::InputDim() const {
        return diff_variables.size();
    }

    std::vector<float> MultiViewEnergy::GetStatus() const {
        // Status is just the current value of all the value we can differentiate with respect to
        std::vector<float> status(diff_variables.size(), 0.f);
        for (size_t i = 0; i < status.size(); i++) {
            status[i] = diff_variables[i]->GetValue();
        }

        return status;
    }

    std::string MultiViewEnergy::ToString() const {
        return std::string("");
    }

} // drdemo namespace
