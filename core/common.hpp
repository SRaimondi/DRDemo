//
// Created by simon on 12.05.17.
//

#ifndef DRDEMO_COMMON_HPP
#define DRDEMO_COMMON_HPP

namespace drdemo {

    // Small tolerance value used in the ray tracing part
#define EPS 10e-3f

    // Define some constants used in the ray marching process
#define MAX_STEPS 1000
#define MIN_DIST 0.001f
#define MAX_DIST 100.f

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

    // Return squared number if positive, else return zero
    template<typename T>
    static inline T PositiveSQ(T x) {
        return x > (T) 0 ? x * x : (T) 0;
    }

    // Return squared number if negative, else return zero
    template<typename T>
    static inline T NegativeSQ(T x) {
        return x < (T) 0 ? x * x : (T) 0;
    }

    // Lerp
    template<typename T>
    static inline T Lerp(T t, T a, T b) {
        return ((T) 1 - t) * a + t * b;
    }

} // drdemo namespace

#endif //DRDEMO_COMMON_HPP
