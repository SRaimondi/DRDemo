//
// Created by simon on 18.04.17.
//

#ifndef DRDEMO_IMAGE_HPP
#define DRDEMO_IMAGE_HPP

/**
 * This file defines the final image that is created from the rendering system
 */

#include <fstream>
#include "diffobject.hpp"

namespace rt {

    template<size_t WIDTH, size_t HEIGHT>
    class Image {
    private:
        // Vector containing the RGB elements of the image
        utils::Vector<ad::Vec3F, WIDTH * HEIGHT> content;

    public:
        // Create new empty image, set all content to black
        Image();

        // Access element at given indices
        inline ad::Vec3F const &operator()(size_t i, size_t j) const {
            return content[i * WIDTH + j];
        }

        inline ad::Vec3F &operator()(size_t i, size_t j) {
            return content[i * WIDTH + j];
        }

        // Access image content as vector
        inline utils::Vector<ad::Vec3F, WIDTH * HEIGHT> const &Content() const {
            return content;
        }

        // Create ppm file from image
        void CreatePPM(std::string const &file_name) const;
    };

    template<size_t WIDTH, size_t HEIGHT>
    Image<WIDTH, HEIGHT>::Image()
            : content(ad::Vec3F(0.f)) {}

    template<size_t WIDTH, size_t HEIGHT>
    void Image<WIDTH, HEIGHT>::CreatePPM(std::string const &file_name) const {
        // Create output file
        std::fstream file = std::fstream(file_name, std::fstream::out);
        // Write header
        file << "P6\n";
        file << WIDTH << " " << HEIGHT << "\n";
        file << "255\n";
        // Write colors to file
        for (int i = static_cast<int>(HEIGHT - 1); i >= 0; i--) {
            for (size_t j = 0; j < WIDTH; j++) {
                ad::Vec3F c = (*this)(i, j);
                unsigned char r = static_cast<unsigned char> (c[0].Value() * 255.f);
                unsigned char g = static_cast<unsigned char> (c[1].Value() * 255.f);
                unsigned char b = static_cast<unsigned char> (c[2].Value() * 255.f);
                // Output color
                file << r << g << b;
            }
        }
        file.close();
    }

} // rt namespace

#endif //DRDEMO_IMAGE_HPP
