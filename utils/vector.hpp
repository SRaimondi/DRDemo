//
// Created by simon on 14.04.17.
//

#ifndef DRDEMO_VECTOR_HPP
#define DRDEMO_VECTOR_HPP

/**
 * Generic simple fixed size vector template that supports math operations
 */

#include <cstddef>
#include <cassert>
#include <algorithm>
#include "promote.hpp"

namespace utils {

    template<typename T, size_t SIZE>
    class Vector {
    protected:
        // Heap allocated vector content
        T *elements;

    public:
        Vector();

        template<typename E>
        explicit Vector(E const &v);

        // Copy constructor
        Vector(Vector<T, SIZE> const &other);

        // Assignment operator
        Vector<T, SIZE> &operator=(Vector<T, SIZE> const &other);

        // Move constructor
        Vector(Vector<T, SIZE> &&other);

        // Move operator
        Vector<T, SIZE> &operator=(Vector<T, SIZE> &&other);

        virtual ~Vector();

        // Get vector size
        constexpr inline size_t Size() const {
            return SIZE;
        }

        // Access vector element
        inline T const &operator[](size_t i) const {
#ifdef DEBUG
            assert(i < SIZE);
#endif
            return elements[i];
        }

        inline T &operator[](size_t i) {
#ifdef DEBUG
            assert(i < SIZE);
#endif
            return elements[i];
        }

        // Math operators on itself
        template<typename E>
        Vector<T, SIZE> &operator+=(Vector<E, SIZE> const &a) {
            for (size_t i = 0; i < SIZE; i++) {
                elements[i] += a.elements[i];
            }

            return *this;
        }

        template<typename E>
        Vector<T, SIZE> &operator-=(Vector<E, SIZE> const &a) {
            for (size_t i = 0; i < SIZE; i++) {
                elements[i] -= a.elements[i];
            }

            return *this;
        }

        template<typename E>
        Vector<T, SIZE> &operator*=(Vector<E, SIZE> const &a) {
            for (size_t i = 0; i < SIZE; i++) {
                elements[i] *= a.elements[i];
            }

            return *this;
        }

        template<typename E>
        Vector<T, SIZE> &operator*=(E const &a) {
            for (size_t i = 0; i < SIZE; i++) {
                elements[i] *= a;
            }

            return *this;
        }

        template<typename E>
        Vector<T, SIZE> &operator/=(Vector<E, SIZE> const &a) {
            for (size_t i = 0; i < SIZE; i++) {
#ifdef DEBUG
                assert(a.elements[i] != E(0));
#endif
                elements[i] /= a.elements[i];
            }

            return *this;
        }

        template<typename E>
        Vector<T, SIZE> &operator/=(E const &a) {
#ifdef DEBUG
            assert(a != E(0));
#endif
            for (size_t i = 0; i < SIZE; i++) {
                elements[i] /= a;
            }

            return *this;
        }
    };

    template<typename T, size_t SIZE>
    Vector<T, SIZE>::Vector()
            : elements(new T[SIZE]) {}

    template<typename T, size_t SIZE>
    template<typename E>
    Vector<T, SIZE>::Vector(E const &v)
            : elements(new T[SIZE]) {
        for (size_t i = 0; i < SIZE; i++) {
            elements[i] = v;
        }
    }

    template<typename T, size_t SIZE>
    Vector<T, SIZE>::Vector(Vector<T, SIZE> const &other)
            : elements(new T[SIZE]) {
        std::copy(other.elements, other.elements + SIZE, elements);
    }

    template<typename T, size_t SIZE>
    Vector<T, SIZE> &Vector<T, SIZE>::operator=(Vector<T, SIZE> const &other) {
        if (this != &other) {
            for (size_t i = 0; i < SIZE; ++i) {
                elements[i] = other.elements[i];
            }
        }

        return *this;
    }

    template<typename T, size_t SIZE>
    Vector<T, SIZE>::Vector(Vector<T, SIZE> &&other)
            : elements(other.elements) {
        other.elements = nullptr;
    }

    template<typename T, size_t SIZE>
    Vector<T, SIZE> &Vector<T, SIZE>::operator=(Vector<T, SIZE> &&other) {
        if (this != &other) {
            delete[] elements;
            elements = other.elements;
            other.elements = nullptr;
        }

        return *this;
    }

    template<typename T, size_t SIZE>
    Vector<T, SIZE>::~Vector() {
        if (elements != nullptr) {
            delete[] elements;
        }
    }

    /**
     * Declare the mathematical operators and function that allow to write expression using Vector<T> and Vector<E> and
     * Vector<T> / E
     */

    // Sum of Vector<T> and Vector<E>
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator+(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    // Subtraction of Vector<T> and Vector<E>
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator-(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a[i] - b[i];
        }
        return c;
    }

    // Multiplication of Vector<T> and Vector<E>
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator*(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a[i] * b[i];
        }
        return c;
    }

    // Multiplication of Vector<T> and E
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>, SIZE>
    operator*(Vector<T, SIZE> const &a, E const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a[i] * b;
        }
        return c;
    }

    // Multiplication of T and Vector<E>
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator*(T const &a, Vector<E, SIZE> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a * b[i];
        }
        return c;
    }

    // Division of Vector<T> and Vector<E>
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator/(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
#ifdef DEBUG
            assert(b[i] != E(0));
#endif
            c[i] = a[i] / b[i];
        }
        return c;
    }

    // Division of Vector<T> and E
    template<typename T, typename E, size_t SIZE>
    Vector<typename traits::Promote<T, E>::TResult, SIZE>
    operator/(Vector<T, SIZE> const &a, E const &b) {
#ifdef DEBUG
        assert(b != E(0));
#endif
        Vector<typename traits::Promote<T, E>::TResult, SIZE> c;
        for (size_t i = 0; i < SIZE; i++) {
            c[i] = a[i] / b;
        }
        return c;
    }

    /**
     * Define Vector mathematical functions
     */

    // Dot product
    template<typename T, typename E, size_t SIZE>
    typename traits::Promote<T, E>::TResult
    Dot(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        typename traits::Promote<T, E>::TResult dot = a[0] * b[0];
        for (size_t i = 1; i < SIZE; i++) {
            dot += a[i] * b[i];
        }
        return dot;
    }

    // Absolute value of dot function, the abs function we want to use must be passed as argument
    template<typename ABS_F, typename T, typename E, size_t SIZE>
    T AbsDot(Vector<T, SIZE> const &a, Vector<E, SIZE> const &b) {
        return ABS_F(Dot(a, b));
    }

    // Cross product for vector of SIZE = 3
    template<typename T, typename E>
    Vector<typename traits::Promote<T, E>::TResult, 3>
    Cross(Vector<T, 3> const &a, Vector<E, 3> const &b) {
        Vector<typename traits::Promote<T, E>::TResult, 3> cross;

        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[2] * b[0] - a[0] * b[2];
        cross[2] = a[0] * b[1] - a[1] * b[0];

        return cross;
    }

    // Length² of vector
    template<typename T, size_t SIZE>
    T Length2(Vector<T, SIZE> const &a) {
        return Dot(a, a);
    }

    // Length of vector
    template<typename SQRT_F, typename T, size_t SIZE>
    T Length(Vector<T, SIZE> const &a) {
        return SQRT_F(Length2(a));
    }

    /**
     * Create subclass of Vector<T, 3> as Vector3 for practice
     */
    template<typename T>
    class Vector3 : public Vector<T, 3> {
    public:
        Vector3();

        // Three value constructor
        template<typename E>
        Vector3(E const &x, E const &y, E const &z);

        // Single value constructor
        template<typename E>
        Vector3(E const &v);

        // Copy
        template<typename E>
        Vector3(Vector3<E> const &other);

        // Assignment operator
        Vector3<T> &operator=(Vector3<T> const &other);
    };

    template<typename T>
    Vector3<T>::Vector3()
            : Vector<T, 3>() {
        this->elements[0] = T(0);
        this->elements[1] = T(0);
        this->elements[2] = T(0);
    }

    template<typename T>
    template<typename E>
    Vector3<T>::Vector3(E const &x, E const &y, E const &z)
            : Vector<T, 3>() {
        this->elements[0] = x;
        this->elements[1] = y;
        this->elements[2] = z;
    }

    template<typename T>
    template<typename E>
    Vector3<T>::Vector3(E const &v)
            : Vector<T, 3>(v) {}

    template<typename T>
    template<typename E>
    Vector3<T>::Vector3(Vector3<E> const &other)
            : Vector<T, 3>(other) {}

    template<typename T>
    Vector3<T> &Vector3<T>::operator=(Vector3<T> const &other) {
        Vector<T, 3>::operator=(other);

        return *this;
    }

} // utils namespace

#endif //DRDEMO_VECTOR_HPP
