//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_PINHOLE_CAMERA_HPP
#define DRDEMO_PINHOLE_CAMERA_HPP

#include "camera.hpp"

namespace drdemo {

    /**
     * Define Pinhole Camera class
     */
    class PinholeCamera : public CameraInterface {
    private:
        // Position
        Vector3F eye_world;
        // Local base vectors
        Vector3F u, v, w;
        // Field of view
        float bottom, top, left, right;
        // Film size
        u_int32_t width, height;

    public:
        PinholeCamera(Vector3F const &e, Vector3F const &at, Vector3F const &up,
                      float fov, u_int32_t width, u_int32_t height);

        Ray GenerateRay(u_int32_t i, u_int32_t j, float s_x, float s_y) const override;
    };

} // drdemo namespace

#endif //DRDEMO_PINHOLE_CAMERA_HPP
