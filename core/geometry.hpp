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

        // TODO Review this
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

        template<typename U>
        Vector3<T> operator/(const Vector3<U> &v) const {
            return Vector3<T>(x / v.x, y / v.y, z / v.z);
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
        const T length = Length(v);
        assert(length != 0.f);
        return v / length;
    }

    // Minimum vector between two
    template<typename T>
    inline Vector3<T> Min(Vector3<T> const &v1, Vector3<T> const &v2) {
        return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    template<>
    inline Vector3<Float> Min(Vector3<Float> const &v1, Vector3<Float> const &v2) {
        return Vector3<Float>(Min(v1.x, v2.x), Min(v1.y, v2.y), Min(v1.z, v2.z));
    }

    // Maximum vector between two
    template<typename T>
    inline Vector3<T> Max(Vector3<T> const &v1, Vector3<T> const &v2) {
        return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    template<>
    inline Vector3<Float> Max(Vector3<Float> const &v1, Vector3<Float> const &v2) {
        return Vector3<Float>(Max(v1.x, v2.x), Max(v1.y, v2.y), Max(v1.z, v2.z));
    }

    // Convert Float vector to float vector
    inline Vector3<float> Tofloat(Vector3<Float> const &v) {
        return Vector3<float>(v.x.GetValue(), v.y.GetValue(), v.z.GetValue());
    }

    // Convert float vector to Float vector
    inline Vector3<Float> ToFloat(Vector3<float> const &v) {
        return Vector3<Float>(v.x, v.y, v.z);
    }

    // Print vector
    template<typename T>
    std::ostream &operator<<(std::ostream &os, Vector3<T> const &v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";

        return os;
    }

    /**
     * Ray class, attribute are made public for practice
     */
    class Ray {
    public:
        // Origin
        Vector3F o;

        // Direction, NOT forced to be normalized
        Vector3F d;

        // Minimum and maximum parameter of the ray
        mutable Float t_min;
        mutable Float t_max;

        // Ray direction sign
        int sign[3];

        // Constructor
        Ray() : t_min(EPS), t_max(INFINITY) {}

        Ray(Vector3F const &o, Vector3F const &d, float t_min = EPS, float t_max = INFINITY)
                : o(o), d(d), t_min(t_min), t_max(t_max) {
            sign[0] = d.x < 0.f;
            sign[1] = d.y < 0.f;
            sign[2] = d.z < 0.f;
        }

        template<typename T>
        inline Vector3F operator()(T const &t) const { return o + t * d; }
    };

} // drdemo namespace

#endif //DRDEMO_GEOMETRY_HPP
