//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_INTEGRATOR_HPP
#define DRDEMO_INTEGRATOR_HPP

#include "scene.hpp"
#include "spectrum.hpp"
#include "camera.hpp"

namespace drdemo {

    /**
     * Define Integrator interface
     */
    class IntegratorInterface {

    };

    /**
     * Define SurfaceIntegrator base class
     */
    class SurfaceIntegratorInterace : public IntegratorInterface {
    public:
        // Compute incoming radiance for given ray
        virtual Spectrum
        IncomingRadiance(Ray const &ray, Scene const &scene, const CameraInterface &camera, size_t depth) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_INTEGRATOR_HPP
