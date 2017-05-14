//
// Created by simon on 05.05.17.
//

#ifndef DRDEMO_RAD_HPP
#define DRDEMO_RAD_HPP

/**
 * Reverse automatic differentiation tools
 */

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>

namespace drdemo {

    /**
     * Define the TapeNode struct, used internally to build the tape that allows to compute
     * the derivatives in reverse mode
     */
    struct TapeNode {
        // Weight associated with the respective nodes
        float weights[2];
        // Parent node indices in the tape
        size_t parent_i[2];

        TapeNode(float w1, size_t p1, float w2, size_t p2) noexcept;
    };

    /**
     * Define the Tape class which holds the computation progress and allows then to compute the derivatives
     */
    class Tape {
    private:
        // List of Tape nodes
        std::vector<TapeNode> nodes;
        // Boolean flag to check if the Tape is enabled or not
        bool enabled;

    public:
        // Tape default constructor
        Tape();

        // Enable / Disable tape
        inline void Enable() { enabled = true; }

        inline void Disable() { enabled = false; }

        // Check if the Tape is enabled
        inline bool IsEnabled() const { return enabled; }

        // Access TapeNode at given index
        inline TapeNode const &At(size_t index) const {
            return nodes.at(index);
        }

        // Get size of the Tape
        inline size_t Size() const { return nodes.size(); }

        // Clear tape starting from a given index
        void Clear(size_t starting_index);

        // Push a Zero value node (leaf node), returns the index of the node on the Tape
        size_t PushLeaf();

        // Push a node that depends only on another single node given the value of the node and the parent index
        size_t PushSingleNode(float w, size_t p);

        // Push a node that depends on two children given both the values and the parents indices
        size_t PushTwoNode(float w1, size_t p1, float w2, size_t p2);
    };

    // Declare extern Tape variable
    extern Tape default_tape;

    /**
     * Define Float class that allows to do classical float computations while building the tape
     * structure for the reverse differentiation process
     */
    class Float {
    private:
        // Actual value
        float value;
        // Index to the node in the tape
        size_t node_index;

    public:
        // Value constructor, default value is 0
        explicit Float(float v = 0.f);

        // Construct Float given index on tape and value
        Float(size_t index, float v) noexcept;

        Float(Float const &other) noexcept;

        // Assignment operator, default is to create an alias of the variable, not pushing a new one on the stack
        Float &operator=(Float const &other);

        // Set Float value
        Float &operator=(float v) noexcept;

        // Convert to built-in float type
        inline float Value() const noexcept { return value; }

        // Get node index
        inline size_t NodeIndex() const noexcept { return node_index; }

        // Math operators

        // Negation
        Float operator-() const {
            return Float(default_tape.PushSingleNode(-1.f, node_index), -value);
        }

        // Addition
        Float operator+(Float const &v) const {
            return Float(default_tape.PushTwoNode(1.f, node_index, 1.f, v.NodeIndex()), value + v.Value());
        }

        Float operator+(float v) {
            return Float(default_tape.PushSingleNode(1.f, node_index), value + v);
        }

        Float &operator+=(Float const &v) {
            *this = *this + v;
            return *this;
        }

        Float &operator++() {
            value += 1.f;
            return *this;
        }

        Float operator++(int) {
            Float result(*this);
            ++(*this);
            return result;
        }

        // Subtraction
        Float operator-(Float const &v) const {
            return Float(default_tape.PushTwoNode(1.f, node_index, -1.f, v.NodeIndex()), value - v.Value());
        }

        Float operator-(float v) const {
            return Float(default_tape.PushSingleNode(1.f, node_index), value - v);
        }

        Float &operator-=(Float const &v) {
            *this = *this - v;
            return *this;
        }

        Float &operator--() {
            value -= 1.f;
            return *this;
        }

        Float operator--(int) {
            Float result(*this);
            --(*this);
            return result;
        }

        // Multiplication
        Float operator*(Float const &v) const {
            return Float(default_tape.PushTwoNode(v.Value(), node_index, value, v.NodeIndex()),
                         value * v.Value());
        }

        Float operator*(float v) const {
            return Float(default_tape.PushSingleNode(v, node_index), value * v);
        }

        // Division
        Float operator/(Float const &v) const {
            assert(v.Value() != 0.f);
            // If f(a,b) = a/b, then df/da = 1/b and df/db = -a/(b*b)
            return Float(default_tape.PushTwoNode(1.f / v.Value(), node_index,
                                                  -value / (v.Value() * v.Value()), v.NodeIndex()),
                         value / v.Value());
        }

        Float operator/(float v) const {
            assert(v != 0.f);
            return Float(default_tape.PushSingleNode(1.f / v, node_index), value / v);
        }
    };

    std::ostream &operator<<(std::ostream &os, Float const &var);

    /**
     * Declare all the Float mathematical operators
     */

    // Sum of float and Float
    inline Float operator+(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(1.f, b.NodeIndex()), a + b.Value());
    }

    // Subtraction of float and Float
    inline Float operator-(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(-1.f, b.NodeIndex()), a - b.Value());
    }

    // Multiplication of float and Float
    inline Float operator*(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(a, b.NodeIndex()), a * b.Value());
    }

    // Division of float and Float
    inline Float operator/(float a, Float const &b) {
        assert(b.Value() != 0.f);
        return Float(default_tape.PushSingleNode(-a / (b.Value() * b.Value()), b.NodeIndex()), a / b.Value());
    }

    // Comparison operator

    // <
    inline bool operator<(Float const &a, Float const &b) noexcept {
        return a.Value() < b.Value();
    }

    inline bool operator<(Float const &a, float b) noexcept {
        return a.Value() < b;
    }

    inline bool operator<(float a, Float const &b) noexcept {
        return a < b.Value();
    }

    // <=
    inline bool operator<=(Float const &a, Float const &b) noexcept {
        return a.Value() <= b.Value();
    }

    inline bool operator<=(Float const &a, float b) noexcept {
        return a.Value() <= b;
    }

    inline bool operator<=(float a, Float const &b) noexcept {
        return a <= b.Value();
    }

    // >
    inline bool operator>(Float const &a, Float const &b) noexcept {
        return a.Value() > b.Value();
    }

    inline bool operator>(Float const &a, float b) noexcept {
        return a.Value() > b;
    }

    inline bool operator>(float a, Float const &b) noexcept {
        return a > b.Value();
    }

    // >=
    inline bool operator>=(Float const &a, Float const &b) noexcept {
        return a.Value() >= b.Value();
    }

    inline bool operator>=(Float const &a, float b) noexcept {
        return a.Value() >= b;
    }

    inline bool operator>=(float a, Float const &b) noexcept {
        return a >= b.Value();
    }

    // ==
    inline bool operator==(Float const &a, Float const &b) noexcept {
        return a.Value() == b.Value();
    }

    inline bool operator==(Float const &a, float b) noexcept {
        return a.Value() == b;
    }

    inline bool operator==(float a, Float const &b) noexcept {
        return a == b.Value();
    }

    // !=
    inline bool operator!=(Float const &a, Float const &b) noexcept {
        return a.Value() != b.Value();
    }

    inline bool operator!=(Float const &a, float b) noexcept {
        return a.Value() != b;
    }

    inline bool operator!=(float a, Float const &b) noexcept {
        return a != b.Value();
    }

    // Sign of Float
    inline int Sign(Float const &v) noexcept {
        return static_cast<int>((0.f < v.Value()) - (v.Value() < 0.f));
    }

    // Sin of Float
    inline Float Sin(Float const &v) {
        return Float(default_tape.PushSingleNode(std::cos(v.Value()), v.NodeIndex()), std::sin(v.Value()));
    }

    // Cos of Float
    inline Float Cos(Float const &v) {
        return Float(default_tape.PushSingleNode(-std::sin(v.Value()), v.NodeIndex()), std::cos(v.Value()));
    }

    // Tan of Float
    inline Float Tan(Float const &v) {
        return Float(default_tape.PushSingleNode(2.f / (std::cos(2.f * v.Value()) + 1.f), v.NodeIndex()),
                     std::tan(v.Value()));
    }

    // Exp of Float
    inline Float Exp(Float const &v) {
        return Float(default_tape.PushSingleNode(std::exp(v.Value()), v.NodeIndex()), std::exp(v.Value()));
    }

    // Log of Float
    inline Float Log(Float const &v) {
        assert(v.Value() > 0.f);
        return Float(default_tape.PushSingleNode(1.f / v.Value(), v.NodeIndex()), std::log(v.Value()));
    }

    // Pow of Float
    inline Float Pow(Float const &v, float k) {
        return Float(default_tape.PushSingleNode(k * std::pow(v.Value(), k - 1.f), v.NodeIndex()),
                     std::pow(v.Value(), k));
    }

    // Sqrt of Float
    inline Float Sqrt(Float const &v) {
        return Float(default_tape.PushSingleNode(0.5f / std::sqrt(v.Value()), v.NodeIndex()), std::sqrt(v.Value()));
    }

    // Abs of Float
    inline Float Abs(Float const &v) {
        assert(v.Value() != 0.f);
        return Float(default_tape.PushSingleNode(Sign(v), v.NodeIndex()), std::abs(v.Value()));
    }

    // Max of two Float
    inline Float Max(Float const &a, Float const &b) {
        return (a.Value() > b.Value() ? a : b);
    }

    // Min of two Float
    inline Float Min(Float const &a, Float const &b) {
        return (a.Value() < b.Value() ? a : b);
    }

} // drdemo namespace

#endif //DRDEMO_RAD_HPP
