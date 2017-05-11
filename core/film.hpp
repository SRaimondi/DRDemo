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
        const uint32_t width, height;

    public:
        Film(uint32_t w, uint32_t h);

        virtual ~Film();

        // Access width and height
        inline uint32_t Width() const noexcept { return width; }

        inline uint32_t Height() const noexcept { return height; }

        // Add sample to the film
        virtual bool AddSample(Spectrum const &s, uint32_t i, uint32_t j, float s_x, float s_y) = 0;

        // Get final film color at given pixel
        virtual Spectrum At(uint32_t i, uint32_t j) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_FILM_HPP
