//
// Created by simon on 17.04.17.
//

#ifndef DRDEMO_COMMON_HPP
#define DRDEMO_COMMON_HPP

namespace utils {

    // Sign function
    template<typename T>
    inline int Sign(T val) {
        return (T(0) < val) - (val < T(0));
    }

} // utils namespace

#endif //DRDEMO_COMMON_HPP
