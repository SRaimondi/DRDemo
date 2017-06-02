//
// Created by Simone Raimondi on 23.05.17.
//

#ifndef DRDEMO_TAPE_STORAGE_HPP
#define DRDEMO_TAPE_STORAGE_HPP

#include <cstdlib>
#include <cassert>
#include <algorithm>

namespace drdemo {

    /**
     * Define custom Tape storage class, used in the RAD toolbox for th Tape class
     */
    template<typename T>
    class TapeStorage {
    private:
        // Maximum storage
        size_t max_size;
        // Used storage
        size_t size;
        // Allocated memory
        T *data;

    public:
        // Default constructor
        explicit TapeStorage(size_t start_size);

        // Copy constructor
        TapeStorage(TapeStorage<T> const &other);

        // Move constructor
        TapeStorage(TapeStorage<T> &&other);

        // Assignment operator
        TapeStorage &operator=(TapeStorage<T> const &other);

        // Move operator
        TapeStorage &operator=(TapeStorage<T> &&other);

        // Destructor
        ~TapeStorage();

        // Access element by index
        inline T const &operator[](size_t i) const noexcept { return data[i]; }

        inline T &operator[](size_t i) noexcept { return data[i]; }

        // Access element with boundaries check
        inline T const &At(size_t i) const {
            assert(i < size);
            return data[i];
        }

        inline T &At(size_t i) {
            assert(i < size);
            return data[i];
        }

        // Size of the Tape
        inline size_t Size() const noexcept { return size; }

        // Add element
        void Append(T const &element);

        // Cut a chunk of the TapeStorage starting from a given index included
        void Cut(size_t start_index);

        // Resize allocated memory to fit current content
        void Resize();
    };

    template<typename T>
    TapeStorage<T>::TapeStorage(size_t start_size)
            : max_size(start_size), size(0), data(new T[max_size]) {}

    template<typename T>
    TapeStorage<T>::TapeStorage(TapeStorage<T> const &other)
            : max_size(other.max_size), size(other.size), data(new T[max_size]) {
        std::copy(other.data, other.data + size, data);
    }

    template<typename T>
    TapeStorage<T>::TapeStorage(TapeStorage<T> &&other)
            : max_size(other.max_size), size(other.size) {
        // Acquire memory ownership
        data = other.data;
        // Set other tape to 0
        other.size = 0;
        other.max_size = 0;
        other.data = nullptr;
    }

    template<typename T>
    TapeStorage<T> &TapeStorage<T>::operator=(TapeStorage<T> const &other) {
        if (this != &other) {
            // TODO Review this
            assert(max_size == other.max_size);
            size = other.size;
            std::copy(other.data, other.data + size, data);
        }
        return *this;
    }

    template<typename T>
    TapeStorage<T> &TapeStorage<T>::operator=(TapeStorage<T> &&other) {
        if (this != &other) {
            assert(max_size == other.max_size);
            size = other.size;
            // Delete current content
            delete[] data;
            // Acquire memory ownership
            data = other.data;
            // Set other tape to 0
            other.data = nullptr;
            other.size = 0;
            other.max_size = 0;
        }
        return *this;
    }

    template<typename T>
    TapeStorage<T>::~TapeStorage() {
        if (data != nullptr) {
            delete[] data;
        }
    }

    template<typename T>
    void TapeStorage<T>::Append(T const &element) {
        if (size == max_size) {
            // Double the tape capacity
            max_size *= 2;
            T *new_data = new T[max_size];
            // Copy current data into new data
            std::copy(data, data + size, new_data);
            // Free old data
            delete[] data;
            // Set pointer to new allocated space
            data = new_data;
        }
        // Add new element
        data[size++] = element;
    }

    template<typename T>
    void TapeStorage<T>::Cut(size_t start_index) {
        assert(start_index < size);
        // Set size to new start
        size = start_index;
    }

    template<typename T>
    void TapeStorage<T>::Resize() {
        if (size < max_size) {
            // Set new max size to size
            max_size = size;
            // Allocate new space for Tape
            T *new_data = new T[max_size];
            // Copy data
            std::copy(data, data + size, new_data);
            // Free old data
            delete[] data;
            // Set pointer to new allocated space
            data = new_data;
        }
    }

} // drdemo namespace

#endif //DRDEMO_TAPE_STORAGE_HPP
