//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_RENDERER_HPP
#define DRDEMO_RENDERER_HPP

#include "film.hpp"
#include "scene.hpp"
#include "camera.hpp"

namespace drdemo {

    /**
     * Define Renderer interface class
     */
    class RendererInterface {
    public:
        virtual ~RendererInterface() {}

        // Render scene given a Film, a Scene and a Camera
        virtual void RenderImage(Film *const film, Scene const &scene, CameraInterface const &camera) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_RENDERER_HPP
