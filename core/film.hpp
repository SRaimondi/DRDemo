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
        const unsigned width, height;

    public:
        Film(unsigned w, unsigned h);

        virtual ~Film();

        // Access width and height
        inline unsigned Width() const noexcept { return width; }

        inline unsigned Height() const noexcept { return height; }

        // Add sample to the film
        virtual bool AddSample(Spectrum const &s, unsigned i, unsigned j, float s_x, float s_y) = 0;

        // Get final film color at given pixel
        virtual Spectrum At(unsigned i, unsigned j) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_FILM_HPP
