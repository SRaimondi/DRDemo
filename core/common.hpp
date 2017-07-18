//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_COMMON_HPP
#define DRDEMO_COMMON_HPP

namespace drdemo {

    // Small tolerance value used in the ray tracing part
#define EPS 10e-3f

    template<typename T>
    T Clamp(T val, T min, T max) {
        return (val < min ? min : (val > max ? max : val));
    }

    // Sign template
    template<typename T>
    inline int Sign(T v) noexcept {
        return static_cast<int>((static_cast<T>(0) < v) - (v < static_cast<T>(0)));
    }

    // Convert radians to degree
    template<typename T>
    inline T RadToDeg(T const &rad) { return (rad * (180.f * M_1_PI)); }

    // Convert degree to radians
    template<typename T>
    inline T DegToRad(T const &deg) { return (deg * (M_PI / 180.f)); }

} // drdemo namespace

#endif //DRDEMO_COMMON_HPP
