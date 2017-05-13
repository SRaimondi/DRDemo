//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_BOX_FILM_HPP
#define DRDEMO_BOX_FILM_HPP

#include "film.hpp"

namespace drdemo {

    class BoxFilterFilm : public Film {
    private:
        // Raster
        std::vector<Spectrum> raster;
        // Samples added per-pixel
        std::vector<Float> num_samples;

    public:
        BoxFilterFilm(size_t w, size_t h);

        bool AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) override;

        Spectrum At(size_t i, size_t j) const override;
    };

} // drdemo namespace

#endif //DRDEMO_BOX_FILM_HPP
