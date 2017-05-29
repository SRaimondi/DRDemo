//
// Created by Simon on 29.05.2017.
//

#ifndef DRDEMO_IOFILE_HPP
#define DRDEMO_IOFILE_HPP

#include <string>
#include <vector>

namespace drdemo {

    // Utility function to read file and store string in a vector
    bool ReadFile(std::string const &file_name, std::vector<std::string> &lines);

} // drdemo namespace

#endif //DRDEMO_IOFILE_HPP