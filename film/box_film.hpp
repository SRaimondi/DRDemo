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
        BoxFilterFilm(uint32_t w, uint32_t h);

        bool AddSample(Spectrum const &s, uint32_t i, uint32_t j, float s_x, float s_y) override;

        Spectrum At(uint32_t i, uint32_t j) const override;
    };

} // drdemo namespace

#endif //DRDEMO_BOX_FILM_HPP
