//
// Created by Simon on 12.09.2017.
//

#ifndef DRDEMO_ARMADILLO_RENDER_HPP
#define DRDEMO_ARMADILLO_RENDER_HPP

#include <string>

namespace drdemo {

    /**
     * Render target images for MVS of the armadillo model
     */
    void RenderArmadilloImages(const std::string &viewpoints_file, size_t w, size_t h);

} // drdemo namespace

#endif //DRDEMO_ARMADILLO_RENDER_HPP
