//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_COMMON_HPP
#define DRDEMO_COMMON_HPP

namespace drdemo {

    // Small tollerance value
#define EPS 10e-4f

    template<typename T>
    T Clamp(T val, T min, T max) {
        return (val < min ? min : (val > max ? max : val));
    }

} // drdemo namespace

#endif //DRDEMO_COMMON_HPP
