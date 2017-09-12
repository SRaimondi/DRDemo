//
// Created by Simon on 10.09.2017.
//

#ifndef DRDEMO_TEST_COMMON_HPP
#define DRDEMO_TEST_COMMON_HPP

#include <string>
#include <geometry.hpp>

namespace drdemo {

    /**
     * Load from a file the list of camera viewpoints and return them in a list
     */
    std::vector<Vector3f>
    LoadViewPoints(const std::string &file_name, const char *format, int start_index = 0, int end_index = 0);

} // drdemo namespace

#endif //DRDEMO_TEST_COMMON_HPP
