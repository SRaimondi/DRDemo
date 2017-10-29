//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_LIGHT_HPP
#define DRDEMO_LIGHT_HPP

#include "spectrum.hpp"
#include "interaction.hpp"

namespace drdemo {

    /**
     * Define LightInterface, all lights implementations must inherit from this one
     */
    class LightInterface {
    protected:
        // Bool flag to tell if the light is enabled or not
        // bool is_enabled{true};

        // Number of samples for the light
        int num_samples;

    public:
//        // Check if light is enabled
//        inline bool IsEnabled() const { return is_enabled; }
//
//        // Enable / disable light
//        inline void Enable() { is_enabled = true; }
//
//        inline void Disable() { is_enabled = false; }

        explicit LightInterface(int ns = 1)
                : num_samples(ns) {}

        // Get num,ber of samples
        inline int NumSamples() const { return num_samples; }

        // Sample incoming light at a given Interaction, returns incoming radiance and fills sampling parameters
        virtual Spectrum
        SampleLi(const Interaction &interaction, float u0, float u1, Vector3F *wi, Float *pdf) const = 0;
    };

} // drdemo namespace


#endif //DRDEMO_LIGHT_HPP
