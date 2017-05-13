//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_FILM_HPP
#define DRDEMO_FILM_HPP

#include "spectrum.hpp"
#include "geometry.hpp"

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

        virtual ~Film();

        // Access width and height
        inline size_t Width() const noexcept { return width; }

        inline size_t Height() const noexcept { return height; }

        // Add sample to the film
        virtual bool AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) = 0;

        // Get final film color at given pixel
        virtual Spectrum At(size_t i, size_t j) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_FILM_HPP
