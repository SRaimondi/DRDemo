//
// Created by Simon on 10.09.2017.
//

#include <iofile.hpp>
#include "test_common.hpp"

namespace drdemo {

    std::vector<Vector3f> LoadViewPoints(const std::string &file_name, const char *format,
                                         int start_index, int end_index) {
        std::vector<Vector3f> points;
        // Read file
        std::vector<std::string> file;
        if (ReadFile(file_name, file)) {
            // Parse lines and add points
            float p_x, p_y, p_z;
            // If start index and end index are both equal to zero, load all file
            if (start_index == 0 && end_index == 0) {
                for (auto const &line : file) {
                    if (sscanf(line.c_str(), format, &p_x, &p_y, &p_z) == 3) {
                        // Add point
                        points.emplace_back(Vector3f(p_x, p_y, p_z));
                    } else {
                        std::cerr << "Point has bad format, supposed to be " << format << std::endl;
                    }
                }
            } else {
                // Only load a range
                for (int i = start_index; i < end_index; i++) {
                    if (sscanf(file[i].c_str(), format, &p_x, &p_y, &p_z) == 3) {
                        // Add point
                        points.emplace_back(Vector3f(p_x, p_y, p_z));
                    } else {
                        std::cerr << "Point has bad format, supposed to be " << format << std::endl;
                    }
                }
            }
        } else {
            std::cerr << "Error reading input points file!" << std::endl;
            exit(EXIT_FAILURE);
        }

        return points;
    }

}
