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
        BoxFilterFilm(u_int32_t w, u_int32_t h);

        bool AddSample(Spectrum const &s, u_int32_t i, u_int32_t j, float s_x, float s_y) override;

        Spectrum At(u_int32_t i, u_int32_t j) const override;
    };

} // drdemo namespace

#endif //DRDEMO_BOX_FILM_HPP
