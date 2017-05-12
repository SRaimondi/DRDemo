//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_SIMPLE_RENDERER_HPP
#define DRDEMO_SIMPLE_RENDERER_HPP

#include "renderer.hpp"
#include "integrator.hpp"

namespace drdemo {

    /**
     * Define SimpleRenderer class, integrates one ray per pixel at the center
     */
    class SimpleRenderer : public RendererInterface {
    private:
        // Surface integrator
        const std::shared_ptr<const SurfaceIntegratorInterace> surface_integrator;

    public:
        SimpleRenderer(std::shared_ptr<const SurfaceIntegratorInterace> const &s_i);

        void RenderImage(Film *const film, Scene const &scene, CameraInterface const &camera) const override;
    };

} // dredemo namespace

#endif //DRDEMO_SIMPLE_RENDERER_HPP
