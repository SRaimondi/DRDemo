//
// Created by simon on 09.05.17.
//

#include <cstdint>
#include "box_film.hpp"

namespace drdemo {

    BoxFilterFilm::BoxFilterFilm(size_t w, size_t h)
            : Film(w, h), raster(width * height, Spectrum()), num_samples(width * height) {}

    BoxFilterFilm &BoxFilterFilm::operator=(BoxFilterFilm const &other) {
        if (this != &other) {
            raster = other.raster;
            num_samples = other.num_samples;
        }

        return *this;
    }

    bool BoxFilterFilm::AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) {
        // Check we are inside pixel boundaries
        if (s_x < 0.f || s_x > 1.f || s_y < 0.f || s_y > 1.f) { return false; }
        // Add sample
        raster[j * width + i] += s;
        // Increase number of samples of the pixel
        ++num_samples[j * width + i];

        return true;
    }

    Spectrum BoxFilterFilm::At(size_t i, size_t j) const {
        return (raster[j * width + i] / num_samples[j * width + i]);
    }

    BoxFilterFilm BoxFilterFilm::operator-(BoxFilterFilm const &other) const {
        BoxFilterFilm difference(width, height);

        // Compute difference between the two films
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; ++i) {
                difference.raster[j * width + i] = At(i, j) - other.At(i, j);
                difference.num_samples[j * width + i] = 1.f;
            }
        }

        return difference;
    }

    Float BoxFilterFilm::SquaredNorm() const {
        // Float squared_norm = At(0, 0).r * At(0, 0).r + At(0,0).g * At(0, 0).g + At(0, 0).b * At(0, 0).b;
        Float squared_norm;
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; i++) {
                squared_norm += At(i, j).r * At(i, j).r + At(i, j).g * At(i, j).g + At(i, j).b * At(i, j).b;
            }
        }

        return squared_norm;
    }

    void BoxFilterFilm::Abs() {
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; i++) {
                raster[j * width + i] = drdemo::Abs(raster[j * width + i]);
            }
        }
    }

} // drdemo namespace
