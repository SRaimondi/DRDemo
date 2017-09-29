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

        BoxFilterFilm &operator=(const BoxFilterFilm &other);

        bool AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) override;

        Spectrum const &At(size_t i, size_t j) const override;

        // Compute difference between this film and another one
        BoxFilterFilm operator-(BoxFilterFilm const &other) const;

        // Compute difference between this film and other film's raw data
        BoxFilterFilm operator-(const std::vector<float> &raw_other) const;

        Float SquaredNorm() const override;

        std::vector<float> Raw() const override;

        void Abs() override;

        // Create BoxFilterFilm from given input .png
        static BoxFilterFilm FromPNG(const std::string &file_name);
    };

} // drdemo namespace

#endif //DRDEMO_BOX_FILM_HPP
