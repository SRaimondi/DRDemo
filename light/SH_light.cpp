//
// Created by Simon on 28.10.2017.
//

#include <random>
#include "SH_light.hpp"

namespace drdemo {

    // First 33 factorial
    static const float factorial_data[] = {
            1.f, 2.f, 6.f, 24.f, 120.f, 720.f, 5040.f, 40320.f, 362880.f, 3628800.f, 39916800.f, 479001600.f,
            6227020800.f, 87178291200.f, 1307674368000.f, 20922789888000.f, 355687428096000.f, 6402373705728000.f,
            121645100408832000.f, 2432902008176640000.f, 51090942171709440000.f, 1124000727777607680000.f,
            25852016738884976640000.f, 620448401733239439360000.f, 15511210043330985984000000.f,
            403291461126605635584000000.f, 10888869450418352160768000000.f, 304888344611713860501504000000.f,
            8841761993739701954543616000000.f, 265252859812191058636308480000000.f,
            8222838654177922817725562880000000.f, 263130836933693530167218012160000000.f,
            8683317618811886495518194401280000000.f
    };

    static inline float factorial(int index) {
        assert(index >= 0 && index < 33);
        return factorial_data[index];
    }

    SHSample::SHSample(size_t num_coeff)
            : coeff(num_coeff, 0.f) {}

    float SHLight::P(int l, int m, float x) const {
        float pmm = 1.f;
        if (m > 0) {
            float somx2 = std::sqrt((1.f - x) * (1.f + x));
            float fact = 1.f;
            for (int i = 1; i <= m; i++) {
                pmm *= (-fact) * somx2;
                fact += 2.f;
            }
        }

        if (l == m) { return pmm; }

        float pmmp1 = x * (2.f * m + 1.f) * pmm;
        if (l == m + 1) { return pmmp1; }

        float pll = 0.f;
        for (int ll = m + 2; ll <= l; ++ll) {
            pll = ((2.f * ll - 1.f) * x * pmmp1 - (ll + m - 1.f) * pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pll;
    }

    float SHLight::K(int l, int m) const {
        float temp = ((2.f * l + 1.f) * factorial(l - m)) / (4.f * PI * factorial(l + m));
        return std::sqrt(temp);
    }

    float SHLight::SH(int l, int m, float theta, float phi) const {
        // Return a point sample of a Spherical Harmonic basis function
        // l is the band, range [0..N]
        // m in the range [-l..l]
        // theta in the range [0..Pi]
        // phi in the range [0..2*Pi]
        const float sqrt_2 = std::sqrt(2.f);
        const float cos_theta = std::cos(theta);

        if (m == 0) { return K(l, 0) * P(l, m, cos_theta); }
        else if (m > 0) { return sqrt_2 * K(l, m) * std::cos(m * phi) * P(l, m, cos_theta); }
        else { return sqrt_2 * K(l, -m) * std::sin(-m * phi) * P(l, -m, cos_theta); }
    }

    SHLight::SHLight(int num_bands, int sqrt_num_samples)
            : LightInterface(sqrt_num_samples * sqrt_num_samples),
              samples(num_samples),
              num_bands(num_bands), num_coeff(num_bands * num_bands), coefficients(new Float[num_coeff]),
              used_samples(0) {
        // Array index
        int i = 0;
        // Random number generator
        std::mt19937 generator;
        std::uniform_real_distribution<float> dist(0.f, 1.f);

        float one_over_n = 1.f / (float) sqrt_num_samples;
        for (int a = 0; a < sqrt_num_samples; a++) {
            for (int b = 0; b < sqrt_num_samples; b++) {
                float x = (a + dist(generator)) * one_over_n;
                float y = (b + dist(generator)) * one_over_n;
                // Compute theta and phi
                float theta = 2.f * std::acos(std::sqrt(1.f - x));
                float phi = 2.f * PI * y;
                samples[i].sph = Vector3f(theta, phi, 1.f);
                // Convert to spherical coordinates
                samples[i].dir = Vector3f(std::sin(theta) * std::cos(phi), std::cos(theta),
                                          std::sin(theta) * std::sin(phi));
                // Reserve space
                samples[i].coeff.reserve(num_coeff);
                // Pre-compute SH coefficients for this sample
                for (int l = 0; l < num_bands; ++l) {
                    for (int m = -l; m <= l; ++m) {
                        int index = l * (l + 1) + m;
                        samples[i].coeff[index] = SH(l, m, theta, phi);
                    }
                }
                // Increment linear index
                ++i;
            }
        }
    }

    Spectrum SHLight::SampleLi(const Interaction &interaction, float, float, Vector3F *wi, Float *pdf) const {
//        // Set pdf to 1
//        *pdf = 1.f;
//        // Set light direction as normal, we do the dot scaling directly here
//        *wi = interaction.n;
//        Float sh_value;
//        // Loop over all samples
//        for (const auto &sample : samples) {
//            const Float n_dot_l = interaction.n.x * sample.dir.x +
//                                  interaction.n.y * sample.dir.y +
//                                  interaction.n.z * sample.dir.z;
//            // Check if we are in the same hemisphere
//            if (n_dot_l > 0.f) {
//                for (int i = 0; i < num_coeff; ++i) {
//                    sh_value += n_dot_l * coefficients[i] * sample.coeff[i];
//                }
//            }
//        }
//
//        // Scale value by number of samples
//        sh_value = sh_value / (float) num_samples;
//
//        // Return light contribution
//        return {sh_value, sh_value, sh_value};

        // New method, only return one sample
        *pdf = 1.f; // (4.f * PI);
        // Get current sample
        const auto &sample = samples[used_samples];
        // Set light direction
        wi->x = sample.dir.x;
        wi->y = sample.dir.y;
        wi->z = sample.dir.z;
        // Compute SH value
        Float sh_value = coefficients[0] * sample.coeff[0];
        for (int i = 1; i < num_coeff; ++i) {
            sh_value += coefficients[i] * sample.coeff[i];
        }
        // Increase number of used samples
        used_samples = (used_samples == num_samples - 1) ? 0 : (used_samples + 1);

        return {sh_value, sh_value, sh_value};
    }

    void SHLight::Initialise(const SphericalFunction &func) {
        // Compute the base coefficients of our SH representation given the spherical function
        const float weight = 4.f * PI;
        // Loop over all of our samples
        for (const auto &sample : samples) {
            // Evaluate our coefficients
            for (int n = 0; n < num_coeff; ++n) {
                // We can just set the value here, no need to keep track of the derivatives
                coefficients[n].SetValue(
                        coefficients[n].GetValue() + func(sample.sph.x, sample.sph.y) * sample.coeff[n]);
            }
        }
        // Divide the result by weight and number of samples
        const float factor = weight / (float) num_samples;
        for (int i = 0; i < num_coeff; ++i) {
            coefficients[i].SetValue(coefficients[i].GetValue() * factor);
        }
    }

} // drdemo namespace