//
// Created by Simone Raimondi on 23.05.17.
//

#ifndef DRDEMO_VECTOR_HPP
#define DRDEMO_VECTOR_HPP

#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace drdemo {

    /**
     * Define generic size and type vector class, elements are allocated on the heap
     */
    template<typename T>
    class Vector {
    private:
        // Vector size
        size_t size;
        // Vector elements
        T *elements;

    public:
        // Default constructor
        explicit Vector(size_t size);

        // Copy constructor
        Vector(Vector<T> const &other);

        // Move constructor
        Vector(Vector<T> &&other);

        // Assignment operator
        Vector &operator=(Vector<T> const &other);

        // Move operator
        Vector &operator=(Vector<T> &&other);

        // Destructor
        ~Vector();

        // Element access operator
        inline T const &operator[](size_t i) const noexcept { return elements[i]; }

        inline T &operator[](size_t i) noexcept { return elements[i]; }

        inline T const &At(size_t i) const noexcept {
            assert(i < size);
            return elements[i];
        }

        inline T &At(size_t i) noexcept {
            assert(i < size);
            return elements[i];
        }

        // Size of the Vector
        inline size_t Size() const noexcept { return size; }
    };

    // Vector implementation
    template<typename T>
    Vector<T>::Vector(size_t size)
            : size(size), elements(new T[size]) {}

    template<typename T>
    Vector<T>::Vector(Vector<T> const &other)
            : size(other.size), elements(new T[size]) {
        std::copy(other.elements, other.elements + size, elements);
    }

    template<typename T>
    Vector<T>::Vector(Vector<T> &&other)
            : size(other.size) {
        // Acquire ownership of other vector memory
        elements = other.elements;
        // Set other vector pointer to nullptr
        other.size = 0;
        other.elements = nullptr;
    }

    template<typename T>
    Vector<T> &Vector<T>::operator=(Vector<T> const &other) {
        if (this != &other) {
            assert(size == other.size);
            // Copy memory
            std::copy(other.elements, other.elements + size, elements);
        }
        return *this;
    }

    template<typename T>
    Vector<T> &Vector<T>::operator=(Vector<T> &&other) {
        if (this != &other) {
            assert(size == other.size);
            // Delete current resources
            delete[] elements;
            // Transfer memory ownership
            elements = other.elements;
            // Set other vector to 0
            other.elements = nullptr;
            other.size = 0;
        }
        return *this;
    }

    template<typename T>
    Vector<T>::~Vector() {
        if (elements != nullptr) {
            delete[] elements;
        }
    }

    // Math operators
    template<typename T>
    Vector<T> operator+(Vector<T> const &v1, Vector<T> const &v2) {
        assert(v1.Size() == v2.Size());
        Vector<T> result(v1.Size());
        for (size_t i = 0; i < v1.Size(); ++i) {
            result[i] = v1[i] + v2[i];
        }
        return result;
    }

    template<typename T>
    Vector<T> operator-(Vector<T> const &v1, Vector<T> const &v2) {
        assert(v1.Size() == v2.Size());
        Vector<T> result(v1.Size());
        for (size_t i = 0; i < v1.Size(); ++i) {
            result[i] = v1[i] - v2[i];
        }
        return result;
    }

    template<typename T>
    Vector<T> operator*(Vector<T> const &v1, Vector<T> const &v2) {
        assert(v1.Size() == v2.Size());
        Vector<T> result(v1.Size());
        for (size_t i = 0; i < v1.Size(); ++i) {
            result[i] = v1[i] * v2[i];
        }
        return result;
    }

    template<typename T, typename U>
    Vector<T> operator*(Vector<T> const &v, U const &s) {
        Vector<T> result(v.Size());
        for (size_t i = 0; i < v.Size(); ++i) {
            result[i] = s * v[i];
        }
        return result;
    }

    template<typename U, typename T>
    Vector<T> operator*(U const &s, Vector<T> const &v) {
        Vector<T> result(v.Size());
        for (size_t i = 0; i < v.Size(); ++i) {
            result[i] = s * v[i];
        }
        return result;
    }

    template<typename T, typename U>
    Vector<T> operator/(Vector<T> const &v, U const &s) {
        Vector<T> result(v.Size());
        for (size_t i = 0; i < v.Size(); ++i) {
            result[i] = v[i] / s;
        }
        return result;
    }

    template<typename T>
    Vector<T> operator/(Vector<T> const &v1, Vector<T> const &v2) {
        assert(v1.Size() == v2.Size());
        Vector<T> result(v1.Size());
        for (size_t i = 0; i < v1.Size(); ++i) {
            result[i] = v1[i] / v2[i];
        }
        return result;
    }

    // Dot product
    template<typename T>
    T Dot(Vector<T> const &v1, Vector<T> const &v2) {
        assert(v1.Size() == v2.Size());
        T dot = v1[0] * v2[0];
        for (size_t i = 1; i < v1.Size(); ++i) {
            dot += v1[i] * v2[i];
        }
        return dot;
    }

    // Vector length squared
    template<typename T>
    T Length2(Vector<T> const &v) {
        return Dot(v, v);
    }

} // drdemo namespace

#endif //DRDEMO_VECTOR_HPP
