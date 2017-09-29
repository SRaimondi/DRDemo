//
// Created by simon on 11.05.17.
//

#include "clamp_tonemapper.hpp"
#include "lodepng.hpp"

namespace drdemo {

    void ClampTonemapper::Process(std::string const &file_name, Film const &film) const {
//        // Create output file
//        std::fstream file = std::fstream(file_name, std::fstream::out);
//        // Write ppm header
//        file << "P6\n";
//        file << film.Width() << " " << film.Height() << "\n";
//        file << "255\n";
//
//        Spectrum final_spectrum;
//        // Write colors to file
//        for (int j = static_cast<int>(film.Height() - 1); j >= 0; j--) {
//            for (unsigned i = 0; i < film.Width(); ++i) {
//                final_spectrum = film.At(i, j);
//
//                // Clamp color and output to file
//                unsigned char r = static_cast<unsigned char>(Clamp(final_spectrum.r.GetValue() * 255.f, 0.f, 255.f));
//                unsigned char g = static_cast<unsigned char>(Clamp(final_spectrum.g.GetValue() * 255.f, 0.f, 255.f));
//                unsigned char b = static_cast<unsigned char>(Clamp(final_spectrum.b.GetValue() * 255.f, 0.f, 255.f));
//
//                // Output spectrum to file
//                file << r << g << b;
//            }
//        }
//        // Close file
//        file.close();

        // Create .png image

        // Convert raster to unsigned char vector
        std::vector<unsigned char> image;
        image.reserve(film.Width() * film.Height() * 4);

//        for (int j = static_cast<int>(film.Height()) - 1; j >= 0; j--) {
//            for (unsigned int i = 0; i < film.Width(); i++) {
//                const Spectrum &c = film.At(i, j); // this->operator()(i, j).Clamp(0, 1);
//                image.push_back(static_cast<unsigned char>(Clamp(c.r.GetValue() * 255.f, 0.f, 255.f)));
//                image.push_back(static_cast<unsigned char>(Clamp(c.g.GetValue() * 255.f, 0.f, 255.f)));
//                image.push_back(static_cast<unsigned char>(Clamp(c.b.GetValue() * 255.f, 0.f, 255.f)));
//                // Alpha
//                image.push_back(static_cast<unsigned char>(255.f));
//            }
//        }

        for (unsigned int j = 0; j < film.Height(); ++j) {
            for (unsigned int i = 0; i < film.Width(); ++i) {
                const Spectrum &c = film.At(i, j); // this->operator()(i, j).Clamp(0, 1);
                image.push_back(static_cast<unsigned char>(Clamp(c.r.GetValue() * 255.f, 0.f, 255.f)));
                image.push_back(static_cast<unsigned char>(Clamp(c.g.GetValue() * 255.f, 0.f, 255.f)));
                image.push_back(static_cast<unsigned char>(Clamp(c.b.GetValue() * 255.f, 0.f, 255.f)));
                // Alpha
                image.push_back(static_cast<unsigned char>(255.f));
            }
        }

        // Encode the image
        const unsigned error = lodepng::encode(file_name, image, film.Width(), film.Height());
        // Check if there is an error
        if (error) {
            std::cerr << "Encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
        }
    }

} // drdemo namespace
