//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_INTEGRATOR_HPP
#define DRDEMO_INTEGRATOR_HPP

#include "scene.hpp"
#include "spectrum.hpp"

namespace drdemo {

    /**
     * Define Integrator interface
     */
    class IntegratorInterface {
    public:
        virtual ~IntegratorInterface() {}
    };

    /**
     * Define SurfaceIntegrator base class
     */
    class SurfaceIntegratorInterace : public IntegratorInterface {
    public:
        virtual ~SurfaceIntegratorInterace() {}

        // Compute incoming radiance for given ray
        virtual Spectrum IncomingRadiance(Ray const &ray, Scene const &scene, uint32_t depth) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_INTEGRATOR_HPP
