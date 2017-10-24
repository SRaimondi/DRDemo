//
// Created by Simon on 24.08.2017.
//

#ifndef DRDEMO_SDF_SPHERE_HPP
#define DRDEMO_SDF_SPHERE_HPP

#include <cstdlib>

namespace drdemo {

    /**
     * This two tests try to match the SDF grid to a sphere of radius different from 1 and smaller than 2 and to move the sphere
     */
//    void RadiusSphereTest(int grid_res, float start_radius, float target_radius);
//
//    void RadiusSphereTestMR(int start_grid_res, float start_radius, float target_radius, int ref_steps);
//
//    void MoveSphereTest(int grid_res);

    /**
     * Target render of .obj
     */
//    void OBJRenderTest(int grid_res, const std::string &file_name);

    enum class CAM_CONFIG {
        NORMAL,     // 4 views from y = 0 on the positive and negative sides of the x and z axis
        UP_DOWN,    // 4 views starting from 30 degree on z axis, with deltas of 90 degrees
        ALL         // All 8 vuews
    };

    void OBJRenderTestMR(int start_grid_res, const std::string &file_name, int ref_steps, const CAM_CONFIG &config);

}

#endif //DRDEMO_SDF_SPHERE_HPP
