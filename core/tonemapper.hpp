//
// Created by simon on 11.05.17.
//

#ifndef DRDEMO_TONEMAPPER_HPP
#define DRDEMO_TONEMAPPER_HPP

#include "film.hpp"
#include <string>

namespace drdemo {

    /**
     * Define ToneMapper interface class
     */
    class ToneMapperInterface {
    public:
        // Process film and produce final .png image
        virtual void Process(std::string const &file_name, Film const &film) const = 0;
    };

} // drdemo namespace

#endif //DRDEMO_TONEMAPPER_HPP
