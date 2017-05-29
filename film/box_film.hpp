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
        // std::vector<Float> num_samples;

    public:
        BoxFilterFilm(size_t w, size_t h);

        BoxFilterFilm &operator=(BoxFilterFilm const &other);

        bool AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) override;

        Spectrum const &At(size_t i, size_t j) const override;

        // Compute difference between this film and another one
        BoxFilterFilm operator-(BoxFilterFilm const &other) const;

        Float SquaredNorm() const override;

        void Abs() override;
    };

} // drdemo namespace

#endif //DRDEMO_BOX_FILM_HPP
