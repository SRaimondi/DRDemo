//
// Created by Simon on 01.10.2017.
//

#include "perspective_camera.hpp"

namespace drdemo {

//    void PerspectiveCamera::MapImagePoint(const float *p_image, float *out_p) const {
////        out_p[0] = p_i[0] * p_image[0] + p_i[1] * p_image[1] + p_i[2] * p_image[2];
////        out_p[1] = p_i[3] * p_image[0] + p_i[4] * p_image[1] + p_i[5] * p_image[2];
////        out_p[2] = p_i[6] * p_image[0] + p_i[7] * p_image[1] + p_i[8] * p_image[2];
////        out_p[3] = p_i[9] * p_image[0] + p_i[10] * p_image[1] + p_i[11] * p_image[2];
////        // Normalize
////        for (int i = 0; i < 4; ++i) { out_p[i] /= out_p[3]; }
//        out_p[0] = inv_m[0] * p_image[0] + inv_m[1] * p_image[1] + inv_m[2] * p_image[2];
//        out_p[1] = inv_m[3] * p_image[0] + inv_m[4] * p_image[1] + inv_m[5] * p_image[2];
//        out_p[2] = inv_m[6] * p_image[0] + inv_m[7] * p_image[1] + inv_m[8] * p_image[2];
//    }

    Vector3f PerspectiveCamera::MapImagePoint(const Vector3f &image_p) const {
        Vector3f world_p;
        world_p.x = inv_m[0] * image_p.x + inv_m[1] * image_p.y + inv_m[2] * image_p.z;
        world_p.y = inv_m[3] * image_p.x + inv_m[4] * image_p.y + inv_m[5] * image_p.z;
        world_p.z = inv_m[6] * image_p.x + inv_m[7] * image_p.y + inv_m[8] * image_p.z;

        return world_p;
    }

    PerspectiveCamera::PerspectiveCamera(float *p, const float *c, size_t width, size_t height)
            : width(width), height(height) {
        // Copy memory
        // std::copy(p, p + 12, p_i);
        std::copy(p, p + 9, inv_m);
        // std::copy(c, c + 3, c_w);
        c_w.x = c[0];
        c_w.y = c[1];
        c_w.z = c[2];
        // Compute look at direction
        Vector3f i_center((width - 1) / 2.f + 0.5f, (height - 1) / 2.f + 0.5f, 1.f);
        Vector3f i_center_w = MapImagePoint(i_center);
        look_dir = i_center_w - c_w;
        // Normalize
        look_dir = Normalize(look_dir);
    }

    Ray PerspectiveCamera::GenerateRay(size_t i, size_t j, float s_x, float s_y) const {
        // Convert pixel position to homogeneous coordinates
        // float p_h[3] = {i + s_x, j + s_y, 1.f};
        const Vector3f p_h(i + s_x, j + s_y, 1.f);
        // Compute 3D point using pseudo-inverse
        // float p_w[3];
        const Vector3f p_w = MapImagePoint(p_h);
        // Compute view direction
        Vector3f view_dir = p_w - c_w;
//        view_dir.x = p_w[0] - c_w[0];
//        view_dir.y = p_w[1] - c_w[1];
//        view_dir.z = p_w[2] - c_w[2];
        view_dir = Normalize(view_dir);

        return Ray(ToFloat(c_w), ToFloat(view_dir));
    }

    Vector3F PerspectiveCamera::LookDir() const {
        return ToFloat(look_dir);
    }

} // drdemo namespace