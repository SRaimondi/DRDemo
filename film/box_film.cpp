//
// Created by simon on 09.05.17.
//

#include <cstdint>
#include "box_film.hpp"

namespace drdemo {

    BoxFilterFilm::BoxFilterFilm(size_t w, size_t h)
    // DON'T USE THE CONSTRUCTOR OR THE SAME SPECTRUM GETS COPIED
            : Film(w, h), raster(width * height/* , Spectrum() */) /*, num_samples(width * height) */ {}

    BoxFilterFilm &BoxFilterFilm::operator=(BoxFilterFilm const &other) {
        if (this != &other) {
            raster = other.raster;
            // num_samples = other.num_samples;
        }

        return *this;
    }

    bool BoxFilterFilm::AddSample(Spectrum const &s, size_t i, size_t j, float s_x, float s_y) {
        // Check we are inside pixel boundaries
        if (s_x < 0.f || s_x > 1.f || s_y < 0.f || s_y > 1.f) { return false; }
        // Add sample
        // raster[j * width + i] += s;
        raster[j * width + i] = s;
        // Increase number of samples of the pixel by one
        // num_samples[j * width + i] += 1.f;

        return true;
    }

    Spectrum const &BoxFilterFilm::At(size_t i, size_t j) const {
        // DEBUG ASSERTION
        // assert(num_samples[j * width + i] != 0.f);

        // return (raster[j * width + i] / num_samples[j * width + i]);
        return raster[j * width + i];
    }

    BoxFilterFilm BoxFilterFilm::operator-(BoxFilterFilm const &other) const {
        BoxFilterFilm difference(width, height);

        // Compute difference between the two films
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; ++i) {
                difference.raster[j * width + i] = At(i, j) - other.At(i, j);
                // difference.num_samples[j * width + i] = 1.f;
            }
        }

        return difference;
    }

    Float BoxFilterFilm::SquaredNorm() const {
        Float squared_norm = Norm(At(0, 0));
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; i++) {
                if (j != 0 && i != 0) {
                    squared_norm += Norm(At(i, j));
                }
            }
        }

        return squared_norm;
    }

    void BoxFilterFilm::Abs() {
        for (size_t j = 0; j < height; ++j) {
            for (size_t i = 0; i < width; i++) {
                raster.at(j * width + i) = drdemo::Abs(raster.at(j * width + i));
            }
        }
    }

} // drdemo namespace
