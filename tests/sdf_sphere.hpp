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
    void RadiusSphereTest(int grid_res, float radius);

    void MoveSphereTest(int grid_res);

}

#endif //DRDEMO_SDF_SPHERE_HPP
