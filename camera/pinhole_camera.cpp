//
// Created by simon on 12.05.17.
//

#include "pinhole_camera.hpp"

namespace drdemo {

    PinholeCamera::PinholeCamera(Vector3F const &e, Vector3F const &at, Vector3F const &up,
                                 float fov, size_t width, size_t height)
            : eye_world(e), width(width), height(height) {
        // Compute camera view volume
        top = std::tan(DegToRad(fov / 2.f));
        bottom = -top;
        right = (width / static_cast<float>(height)) * top;
        left = -right;

        // Compute local base
        w = Normalize(e - at);
        u = Normalize(Cross(up, w));
        v = Cross(w, u);
    }

    Ray PinholeCamera::GenerateRay(size_t i, size_t j, float s_x, float s_y) const {
        // Compute point on view plane
        float view_plane_x = left + (right - left) * (i + s_x) / static_cast<float>(width);
        float view_plane_y = top - (top - bottom) * (j + s_y) / static_cast<float>(height);
        Vector3F s = view_plane_x * u + view_plane_y * v - w;

        return Ray(eye_world, Normalize(s));
    }

    Vector3F PinholeCamera::LookDir() const {
        return -w;
    }

}
