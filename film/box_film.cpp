//
// Created by simon on 09.05.17.
//

#include <cstdint>
#include "box_film.hpp"

namespace drdemo {

    BoxFilterFilm::BoxFilterFilm(uint32_t w, uint32_t h)
            : Film(w, h), raster(width * height, Spectrum()), num_samples(width * height) {}

    bool BoxFilterFilm::AddSample(Spectrum const &s, uint32_t i, uint32_t j, float s_x, float s_y) {
        // Check we are inside pixel boundaries
        if (s_x < 0.f || s_x > 1.f || s_y < 0.f || s_y > 1.f) { return false; }
        // Add sample
        raster[j * width + i] += s;
        // Increase number of samples of the pixel
        num_samples[j * width + i]++;

        return true;
    }

    Spectrum BoxFilterFilm::At(uint32_t i, uint32_t j) const {
        return (raster[j * width + i] / num_samples[j * width + i]);
    }

} // drdemo namespace
