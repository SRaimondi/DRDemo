//
// Created by simon on 13.05.17.
//

#ifndef DRDEMO_DIRECT_INTEGRATOR_HPP
#define DRDEMO_DIRECT_INTEGRATOR_HPP

#include "integrator.hpp"

namespace drdemo {

    /**
     * Define DirectIntegrator class, which is a simple integrator that computes direct illumination the scene
     */
    class DirectIntegrator : public SurfaceIntegratorInterace {
    public:
        DirectIntegrator() = default;

        Spectrum IncomingRadiance(Ray const &ray, Scene const &scene, const CameraInterface &camera,
                                  size_t depth) const override;
    };

} // drdemo namespace

#endif //DRDEMO_DIRECT_INTEGRATOR_HPP
