//
// Created by simon on 14.09.17.
//

#ifndef DRDEMO_SDF_LOADING_RENDER_TEST_HPP
#define DRDEMO_SDF_LOADING_RENDER_TEST_HPP

#include <string>

namespace drdemo {

    /**
     * This test tries to load a .sdf file and render it
     */
    void LoadAndTestSDF(const std::string &sdf_file_name, size_t w, size_t h);

} // drdemo namespace

#endif //DRDEMO_SDF_LOADING_RENDER_TEST_HPP
