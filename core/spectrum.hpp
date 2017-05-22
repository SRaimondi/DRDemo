//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_SPECTRUM_HPP
#define DRDEMO_SPECTRUM_HPP

#include "rad.hpp"

namespace drdemo {

    /**
     * Define Spectrum class used for radiance computations
     */
    class Spectrum {
    public:
        Float r, g, b;

        // Constructors
        Spectrum() : r(), g(), b() {}

        // TODO Review this, might cause some problems
//        explicit Spectrum(Float const &v)
//                : r(v), g(v), b(v) {}


        explicit Spectrum(float v)
                : r(v), g(v), b(v) {}

        Spectrum(Float const &r, Float const &g, Float const &b)
                : r(r), g(g), b(b) {}

        Spectrum(float r, float g, float b)
                : r(r), g(g), b(b) {}

        // Math operators on itself
        inline Spectrum &operator+=(Spectrum const &s) {
            r += s.r;
            g += s.g;
            b += s.b;
            return *this;
        }

        inline Spectrum &operator-=(Spectrum const &s) {
            r -= s.r;
            g -= s.g;
            b -= s.b;
            return *this;
        }

        // Spectrum sum
        inline Spectrum operator+(Spectrum const &s) const {
            return Spectrum(r + s.r, g + s.g, b + s.b);
        }

        // Spectrum difference
        inline Spectrum operator-(Spectrum const &s) const {
            return Spectrum(r - s.r, g - s.g, b - s.b);
        }

        // Spectrum scaling
        inline Spectrum operator*(Spectrum const &s) const {
            return Spectrum(r * s.r, g * s.g, b * s.b);
        }

        template<typename U>
        inline Spectrum operator*(U const &t) const {
            return Spectrum(t * r, t * g, t * b);
        }

        template<typename U>
        inline Spectrum operator/(U const &t) const {
            return Spectrum(r / t, g / t, b / t);
        }

        // Check if spectrum is black
        inline bool IsBlack() const {
            return r == 0.f && g == 0.f && b == 0.f;
        }
    };


    // Abs of spectrum
    inline Spectrum Abs(Spectrum const &s) {
        return Spectrum(Abs(s.r), Abs(s.g), Abs(s.b));
    }

    template<typename U>
    inline Spectrum operator*(U const &t, Spectrum const &s) {
        return Spectrum(t * s.r, t * s.g, t * s.b);
    }

} // drdemo namespace

#endif //DRDEMO_SPECTRUM_HPP
