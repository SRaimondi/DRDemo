//
// Created by Simon on 21.08.2017.
//

#ifndef DRDEMO_RECONSTRUCTION_ENERGY_HPP
#define DRDEMO_RECONSTRUCTION_ENERGY_HPP

#include "clamp_tonemapper.hpp"
#include "camera.hpp"
#include "renderer.hpp"
#include "derivative.hpp"
#include "grid.hpp"
#include "scene.hpp"
#include "scalar_function.hpp"

namespace drdemo {

    /**
     * This file contains the implementation of a class that computes the energy term to minimize
     * to match our target images from our SDF
     */
    class ReconstructionEnergy : public ScalarFunctionInterface {
    private:
        // Reference to the Scene to use in the rendering, the scene should only consist of lights and a SDF grid
        Scene &target_scene;
        // Pointer ot the SDF grid of the scene
        std::shared_ptr<SignedDistanceGrid> grid;
        // Reference to the list of target render view
        const std::vector<std::vector<float> > &target_views;
        // Reference to the list of camera used to render the view, order MUST be the same
        const std::vector<std::shared_ptr<const CameraInterface> > &target_cameras;
        // Renderer to be used
        const std::shared_ptr<RendererInterface> renderer;

        // Vector of pointers to all the differentiable variables
        std::vector<Float const *> diff_variables;
        // Class to compute derivatives
        mutable Derivatives derivatives;

        // Regularization term
        const float lambda;

        // Target render resolution
        const size_t width, height;

        // Tonemapper to create images
        ClampTonemapper tonemapper;
        // Current number of function evaluations
        mutable size_t evalutations;

    public:
        // Constructor
        ReconstructionEnergy(Scene &scene,                                                  // Target scene
                             const std::shared_ptr<SignedDistanceGrid> &grid,               // SDF grid in the scene
                             const std::vector<std::vector<float> > &views,                 // Target views
                             const std::vector<std::shared_ptr<const CameraInterface> > &c, // Cameras to use
                             const std::shared_ptr<RendererInterface> &r,                   // Renderer
                             float lambda,                                                  // Normal regularization parameter
                             size_t w, size_t h);                                           // Render image size

        // Rebind differentibale variables
        void RebindVars();

        // Scalar function methods
        size_t InputDim() const override;

        Float Evaluate() const override;

        std::vector<float> ComputeGradient(const Float &out) const override;

        std::vector<float> GetStatus() const override;

        void UpdateStatus(const std::vector<float> &deltas) override;

        void SetStatus(const std::vector<float> &new_status) override;
    };

} // drdemo namespace

#endif //DRDEMO_RECONSTRUCTION_ENERGY_HPP
