//
// Created by simon on 11.05.17.
//

#include "clamp_tonemapper.hpp"
#include <fstream>

namespace drdemo {

    void ClampTonemapper::Process(std::string const &file_name, Film const &film) const {
        // Create output file
        std::fstream file = std::fstream(file_name, std::fstream::out);
        // Write ppm header
        file << "P6\n";
        file << film.Width() << " " << film.Height() << "\n";
        file << "255\n";

        Spectrum final_spectrum;
        // Write colors to file
        for (int j = static_cast<int>(film.Height() - 1); j >= 0; j--) {
            for (unsigned i = 0; i < film.Width(); i++) {
                final_spectrum = film.At(i, j);

                // Clamp color and output to file
                unsigned char r = static_cast<unsigned char>(Clamp(final_spectrum.r.Value() * 255.f, 0.f, 255.f));
                unsigned char g = static_cast<unsigned char>(Clamp(final_spectrum.g.Value() * 255.f, 0.f, 255.f));
                unsigned char b = static_cast<unsigned char>(Clamp(final_spectrum.b.Value() * 255.f, 0.f, 255.f));

                // Output spectrum to file
                file << r << g << b;
            }
        }
        // Close file
        file.close();
    }

} // drdemo namespace
