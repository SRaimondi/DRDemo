//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_FILM_HPP
#define DRDEMO_FILM_HPP

#include "spectrum.hpp"
#include "geometry.hpp"
#include "vector.hpp"

namespace drdemo {

    /**
     * Define Base Film class
     */
    class Film {
    protected:
        // Film size
        const size_t width, height;

    public:
        Film(size_t w, size_t h);

        // Access width and height
        inline size_t Width() const noexcept { return width; }

        inline size_t Height() const noexcept { return height; }

        // Add sample to the film
        virtual bool AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) = 0;

        // Get final film color at given pixel
        virtual Spectrum const &At(size_t i, size_t j) const = 0;

        // Compute the squared norm of the image
        virtual Float SquaredNorm() const = 0;

        // Convert Film to flat r, g, b data
        virtual Vector<float> Raw() const = 0;

        // Compute absolute value of the film
        virtual void Abs() = 0;
    };

} // drdemo namespace

#endif //DRDEMO_FILM_HPP
