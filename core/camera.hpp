//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_CAMERA_HPP
#define DRDEMO_CAMERA_HPP

#include "geometry.hpp"

namespace drdemo {

    /**
     * Define CameraInterface, the base class that all other camera must implement
     */
    class CameraInterface {
    public:
        virtual ~CameraInterface() {}

        // Generate a ray given a pixel coordinates and a sample coordinates
        virtual Ray GenerateRay(uint32_t i, uint32_t j, float s_x, float s_y) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_CAMERA_HPP
