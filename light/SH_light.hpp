//
// Created by Simon on 28.10.2017.
//

#ifndef DRDEMO_SH_LIGHT_HPP
#define DRDEMO_SH_LIGHT_HPP

#include <functional>
#include "light.hpp"

namespace drdemo {

    /**
     * This file defines the first Spherical Harmonics implementation attempt for the project
     * Explanations from http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
     */

    // SH sample struct
    struct SHSample {
        Vector3f sph;
        Vector3f dir;
        std::vector<float> coeff;

        SHSample() = default;

        explicit SHSample(size_t num_coeff);
    };

    // Spherical function typedef
    using SphericalFunction = std::function<float(float, float)>;

    class SHLight : public LightInterface {
    private:
        // Vector of samples
        std::vector<SHSample> samples;
        // Number of bands
        int num_bands;
        // Used number of sample
        int used_samples;

        // Evaluate P function
        float P(int l, int m, float x) const;

        // Renormalisation constant for SH function
        float K(int l, int m) const;

        // SH function
        float SH(int l, int m, float theta, float phi) const;

    public:
        SHLight(int num_bands, int sqrt_num_samples);

        Spectrum SampleLi(const Interaction &interaction, float u0, float u1,
                          Vector3F *wi, Float *pdf) const override;

        // Initialise the SH given a function
        void Initialise(const SphericalFunction &func);
    };

}

#endif //DRDEMO_SH_LIGHT_HPP
