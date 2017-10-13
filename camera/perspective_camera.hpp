//
// Created by Simon on 01.10.2017.
//

#ifndef DRDEMO_PERSPECTIVE_CAMERA_HPP
#define DRDEMO_PERSPECTIVE_CAMERA_HPP

#include "camera.hpp"

namespace drdemo {

    /**
     * Define perspective camera class using projection pseudoinverse matrix
     */
    class PerspectiveCamera : public CameraInterface {
    private:
//        // Projection camera pseudoinverse
//        float p_i[12];

        // Inverse matrix inv(R) * inv(K)
        float inv_m[9];
        // Camera center
        Vector3f c_w;
        // float c_w[3];
        // Image high and width
        const size_t width, height;
        // Look direction
        Vector3f look_dir;

        // Transform 2D point on image plane to 3D points and normalize
        Vector3f MapImagePoint(const Vector3f &image_p) const;

    public:
        // c is a pointer to the pseudo-inverse matrix data in row major order
        PerspectiveCamera(float *p, const float *c, size_t width, size_t height);

        Ray GenerateRay(size_t i, size_t j, float s_x, float s_y) const override;

        Vector3F LookDir() const override; // TODO
    };

} // drdemo namespace

#endif //DRDEMO_PERSPECTIVE_CAMERA_HPP
