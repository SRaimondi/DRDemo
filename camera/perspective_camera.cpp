//
// Created by Simon on 01.10.2017.
//

#include "perspective_camera.hpp"

namespace drdemo {

    void PerspectiveCamera::MapImagePoint(const float *p_image, float *out_p) const {
        out_p[0] = p_i[0] * p_image[0] + p_i[1] * p_image[1] + p_i[2] * p_image[2];
        out_p[1] = p_i[3] * p_image[0] + p_i[4] * p_image[1] + p_i[5] * p_image[2];
        out_p[2] = p_i[6] * p_image[0] + p_i[7] * p_image[1] + p_i[8] * p_image[2];
        out_p[3] = p_i[9] * p_image[0] + p_i[10] * p_image[1] + p_i[11] * p_image[2];
        // Normalize
        for (int i = 0; i < 4; ++i) { out_p[i] /= out_p[3]; }
    }

    PerspectiveCamera::PerspectiveCamera(float *p, float *c, size_t width, size_t height)
            : width(width), height(height) {
        // Copy memory
        std::copy(p, p + 12, p_i);
        std::copy(c, c + 3, c_w);
        // Compute look at direction
        float i_center[3] = {width / 2.f + 0.5f, height / 2.f + 0.5f, 1.f};
        float i_center_w[4];
        MapImagePoint(i_center, i_center_w);
        look_dir.x = i_center_w[0] - c_w[0];
        look_dir.y = i_center_w[1] - c_w[1];
        look_dir.z = i_center_w[2] - c_w[2];
        // Normalize
        look_dir = Normalize(look_dir);
    }

    Ray PerspectiveCamera::GenerateRay(size_t i, size_t j, float s_x, float s_y) const {
        // Convert pixel position to homogeneous coordinates
        float p_h[3] = {i + s_x, j + s_y, 1.f};
        // Compute 3D point using pseudo-inverse
        float p_w[4];
        MapImagePoint(p_h, p_w);
        // Compute view direction
        Vector3F view_dir;
        view_dir.x = p_w[0] - c_w[0];
        view_dir.x = p_w[1] - c_w[1];
        view_dir.x = p_w[2] - c_w[2];

        return Ray(Vector3F(c_w[0], c_w[1], c_w[2]), Normalize(view_dir));
    }

    Vector3F PerspectiveCamera::LookDir() const {
        return ToFloat(look_dir);
    }

} // drdemo namespace