//
// Created by simon on 11.05.17.
//

#ifndef DRDEMO_CLAMP_TONEMAPPER_HPP
#define DRDEMO_CLAMP_TONEMAPPER_HPP

#include "tonemapper.hpp"

namespace drdemo {

    class ClampTonemapper : public ToneMapperInterface {
    public:
        ClampTonemapper() = default;

        void Process(std::string const &file_name, Film const &film) const override;
    };

} // drdemo namespace

#endif //DRDEMO_CLAMP_TONEMAPPER_HPP
