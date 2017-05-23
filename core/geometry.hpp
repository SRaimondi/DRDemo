//
// Created by simon on 08.05.17.
//

#ifndef DRDEMO_GEOMETRY_HPP
#define DRDEMO_GEOMETRY_HPP

#include "rad.hpp"
#include "common.hpp"

namespace drdemo {

    /**
     * Define Vector3 class, components are public for practice and allocated on the stack frame
     */
    template<typename T>
    class Vector3 {
    public:
        // Vector components
        T x, y, z;

        // Constructor
        Vector3()
                : x(0), y(0), z(0) {}

        template<typename U>
        Vector3(U const &x, U const &y, U const &z)
                : x(x), y(y), z(z) {}

        // TODO Review this also
//        template<typename U>
//        Vector3(U const &v)
//                : x(v), y(v), z(v) {}

        Vector3(Vector3<T> const &other)
                : x(other.x), y(other.y), z(other.z) {}

        // Index element access
        T const &operator[](int i) const {
            if (i == 0) { return x; }
            if (i == 1) { return y; }
            return z;
        }

        T &operator[](int i) {
            if (i == 0) { return x; }
            if (i == 1) { return y; }
            return z;
        }

        // Sum
        Vector3<T> operator+(Vector3<T> const &v) const {
            return Vector3<T>(x + v.x, y + v.y, z + v.z);
        }

        Vector3<T> &operator+=(Vector3<T> const &v) {
            *this = *this + v;
            return *this;
        }

        // Subtraction
        Vector3<T> operator-(Vector3<T> const &v) const {
            return Vector3<T>(x - v.x, y - v.y, z - v.z);
        }

        Vector3<T> &operator-=(Vector3<T> const &v) {
            *this = *this - v;
            return *this;
        }

        // Scaling
        template<typename U>
        Vector3<T> operator*(U const &s) const {
            return Vector3<T>(s * x, s * y, s * z);
        }

        template<typename U>
        Vector3<T> &operator*=(U const &s) {
            *this = *this * s;
            return *this;
        }

        template<typename U>
        Vector3<T> operator/(U const &s) const {
            return Vector3<T>(x / s, y / s, z / s);
        }

        template<typename U>
        Vector3<T> operator/=(U const &s) {
            *this = *this / s;
            return *this;
        }

        // Negate vector
        Vector3<T> operator-() const {
            return Vector3<T>(-x, -y, -z);
        }
    };

    // Define common vector types
    using Vector3F = Vector3<Float>;
    using Vector3f = Vector3<float>;

    // Squared length of Vector3
    template<typename T>
    T LengthSquared(Vector3<T> const &v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    // Length of Vector3
    template<typename T>
    inline T Length(Vector3<T> const &v) {
        return std::sqrt(LengthSquared(v));
    }

    template<>
    inline Float Length<Float>(Vector3<Float> const &v) {
        return Sqrt(LengthSquared(v));
    }

    // Scaling operator
    template<typename U, typename T>
    inline Vector3<T> operator*(U const &s, Vector3<T> const &v) {
        return Vector3<T>(s * v.x, s * v.y, s * v.z);
    }

    // Dot product
    template<typename T>
    inline T Dot(Vector3<T> const &a, Vector3<T> const &b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    template<typename T>
    inline T AbsDot(Vector3<T> const &a, Vector3<T> const &b) {
        return std::abs(Dot(a, b));
    }

    template<>
    inline Float AbsDot<Float>(Vector3<Float> const &a, Vector3<Float> const &b) {
        return Abs(a.x * b.x + a.y * b.y + a.z * b.z);
    }

    // Cross product
    template<typename T>
    inline Vector3<T> Cross(Vector3<T> const &a, Vector3<T> const &b) {
        return Vector3<T>(a.y * b.z - a.z * b.y,
                          a.z * b.x - a.x * b.z,
                          a.x * b.y - a.y * b.x);
    }

    // Normalize vector
    template<typename T>
    inline Vector3<T> Normalize(Vector3<T> const &v) {
        return v / Length(v);
    }

    /**
     * Ray class, attribute are made public for practical reasons
     */
    class Ray {
    public:
        // Ray public data, all attributes are differentiable
        // Origin
        Vector3F o;
        // Direction, NOT forced to be normalized
        Vector3F d;
        // Minimum and maximum parameter of the ray
        Float t_min;
        mutable Float t_max;

        // Constructor
        Ray() : t_min(EPS), t_max(INFINITY) {}

        Ray(Vector3F const &o, Vector3F const &d, float t_min = EPS, float t_max = INFINITY)
                : o(o), d(d), t_min(t_min), t_max(t_max) {}

        template<typename T>
        inline Vector3F operator()(T const &t) const { return o + t * d; }
    };

} // drdemo namespace

#endif //DRDEMO_GEOMETRY_HPP
