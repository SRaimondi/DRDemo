//
// Created by Simone Raimondi on 08.08.17.
//

#ifndef DRDEMO_MULTI_VIEW_ENERGY_HPP
#define DRDEMO_MULTI_VIEW_ENERGY_HPP

#include "renderer.hpp"
#include "derivative.hpp"
#include "scene.hpp"
#include "camera.hpp"
#include "scalar_function.hpp"

namespace drdemo {

    /**
     * This file contains the main class of the project
     * It takes a Scene, a list of target images and the respective cameras
     *
     * It then computes the energy as the sum of all the differences between the current view rendering
     * and the target image
     *
     * This can then be used in a minimization procedure
     *
     * Note that this class only provides the gradient for the Shapes in the Scene. It does not provide access to the
     * Cameras and the Light differentiable variables
     */
    class MultiViewEnergy : ScalarFunctionInterface {
    private:
        // Reference to the Scene to use in the rendering
        const Scene &target_scene;
        // Reference to the list of target render view
        const std::vector<std::vector<float> > &target_views;
        // Reference to the list of camera used to render the view, order MUST be the same
        const std::vector<std::shared_ptr<const CameraInterface> > &cameras;
        // Renderer to be used
        const std::shared_ptr<RendererInterface> renderer;

        // Vector of pointers to all the differentiable variables
        std::vector<Float const *> diff_variables;
        // Class to compute derivatives
        mutable Derivatives derivatives;

        // Target render resolution
        const size_t width, height;

    public:
        // Constructor
        MultiViewEnergy(Scene &scene,
                        const std::vector<std::vector<float> > &views,
                        const std::vector<std::shared_ptr<const CameraInterface> > &c,
                        const std::shared_ptr<RendererInterface> &r,
                        size_t w, size_t h);

        // Scalar function methods
        Float Evaluate() const override;

        std::vector<float> ComputeGradient(const Float &out) const override;

        void UpdateStatus(const std::vector<float> &deltas) override;

        void SetStatus(const std::vector<float> &new_status) override;
    };

} // drdemo namespace

#endif //DRDEMO_MULTI_VIEW_ENERGY_HPP
