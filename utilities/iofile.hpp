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

    // Utility function to write a std::string to a file
    void WriteFile(std::string const &file_name, std::string const &content);

} // drdemo namespace

#endif //DRDEMO_IOFILE_HPP
