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
#include <tape_storage.hpp>

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

        TapeNode() = default;

        TapeNode(float w1, size_t p1, float w2, size_t p2) noexcept;
    };

    /**
     * Define the Tape class which holds the computation progress and allows then to compute the derivatives
     */
    class Tape {
    private:
        // List of Tape nodes
        // std::vector<TapeNode> nodes;
        TapeStorage<TapeNode> nodes;
        // List of checkpoints to reset the nodes to a certain point
        std::vector<size_t> clear_indices;
        // Boolean flag to check if the Tape is enabled or not
        // bool enabled;

    public:
        // Tape default constructor
        explicit Tape(size_t starting_size = 1000);

        // Enable / Disable tape
        // TODO This is not working, use Push / Pop
        // inline void Enable() { enabled = true; }

        // inline void Disable() { enabled = false; }

        // Check if the Tape is enabled
        // inline bool IsEnabled() const { return enabled; }

        // Access TapeNode at given index
        inline TapeNode const &At(size_t index) const {
            // return nodes.at(index);
            return nodes.At(index);
        }

        // Get size of the Tape
        inline size_t Size() const {
            // return nodes.size();
            return nodes.Size();
        }

        // Push current size of nodes, can be used to clear after
        inline void Push() {
            // clear_indices.push_back(nodes.size());
            clear_indices.push_back(nodes.Size());
        }

        // Pop current portion of stack, uses last checkpoint saved
        inline void Pop() {
            // nodes.erase(nodes.begin() + clear_indices[clear_indices.size() - 1], nodes.end());
            nodes.Cut(clear_indices[clear_indices.size() - 1]);
            // Delete index
            clear_indices.pop_back();
        }

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
        // GetValue constructor, default value is 0
        explicit Float(float v = 0.f);

        // Construct Float given index on tape and value
        Float(size_t index, float v) noexcept;

        Float(Float const &other) = default;

        // Assignment operator, default is to create an alias of the variable, not pushing a new one on the stack
        Float &operator=(Float const &other);

        // Set Float value
        Float &operator=(float v) noexcept;

        // Convert to built-in float type
        inline float GetValue() const noexcept { return value; }

        // Set value of Float
        inline void SetValue(float v) noexcept { value = v; }

        // Get node index
        inline size_t NodeIndex() const noexcept { return node_index; }

        // Math operators

        // Negation
        Float operator-() const {
            return Float(default_tape.PushSingleNode(-1.f, node_index), -value);
        }

        // Addition
        Float operator+(Float const &v) const {
            return Float(default_tape.PushTwoNode(1.f, node_index, 1.f, v.NodeIndex()), value + v.GetValue());
        }

        Float operator+(float v) {
            return Float(default_tape.PushSingleNode(1.f, node_index), value + v);
        }

        Float &operator+=(Float const &v) {
            *this = *this + v;
            return *this;
        }

        Float &operator+=(float v) {
            *this = *this + v;
            return *this;
        }

        Float &operator++() {
            value += 1.f;
            return *this;
        }

        // Subtraction
        Float operator-(Float const &v) const {
            return Float(default_tape.PushTwoNode(1.f, node_index, -1.f, v.NodeIndex()), value - v.GetValue());
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

        // Multiplication
        Float operator*(Float const &v) const {
            return Float(default_tape.PushTwoNode(v.GetValue(), node_index, value, v.NodeIndex()),
                         value * v.GetValue());
        }

        Float operator*(float v) const {
            return Float(default_tape.PushSingleNode(v, node_index), value * v);
        }

        // Division
        Float operator/(Float const &v) const {
            assert(v.GetValue() != 0.f);
            // If f(a,b) = a/b, then df/da = 1/b and df/db = -a/(b*b)
            return Float(default_tape.PushTwoNode(1.f / v.GetValue(), node_index,
                                                  -value / (v.GetValue() * v.GetValue()), v.NodeIndex()),
                         value / v.GetValue());
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
        return Float(default_tape.PushSingleNode(1.f, b.NodeIndex()), a + b.GetValue());
    }

    // Subtraction of float and Float
    inline Float operator-(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(-1.f, b.NodeIndex()), a - b.GetValue());
    }

    // Multiplication of float and Float
    inline Float operator*(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(a, b.NodeIndex()), a * b.GetValue());
    }

    // Division of float and Float
    inline Float operator/(float a, Float const &b) {
        assert(b.GetValue() != 0.f);
        return Float(default_tape.PushSingleNode(-a / (b.GetValue() * b.GetValue()), b.NodeIndex()), a / b.GetValue());
    }

    // Comparison operator

    // <
    inline bool operator<(Float const &a, Float const &b) noexcept {
        return a.GetValue() < b.GetValue();
    }

    inline bool operator<(Float const &a, float b) noexcept {
        return a.GetValue() < b;
    }

    inline bool operator<(float a, Float const &b) noexcept {
        return a < b.GetValue();
    }

    // <=
    inline bool operator<=(Float const &a, Float const &b) noexcept {
        return a.GetValue() <= b.GetValue();
    }

    inline bool operator<=(Float const &a, float b) noexcept {
        return a.GetValue() <= b;
    }

    inline bool operator<=(float a, Float const &b) noexcept {
        return a <= b.GetValue();
    }

    // >
    inline bool operator>(Float const &a, Float const &b) noexcept {
        return a.GetValue() > b.GetValue();
    }

    inline bool operator>(Float const &a, float b) noexcept {
        return a.GetValue() > b;
    }

    inline bool operator>(float a, Float const &b) noexcept {
        return a > b.GetValue();
    }

    // >=
    inline bool operator>=(Float const &a, Float const &b) noexcept {
        return a.GetValue() >= b.GetValue();
    }

    inline bool operator>=(Float const &a, float b) noexcept {
        return a.GetValue() >= b;
    }

    inline bool operator>=(float a, Float const &b) noexcept {
        return a >= b.GetValue();
    }

    // ==
    inline bool operator==(Float const &a, Float const &b) noexcept {
        return a.GetValue() == b.GetValue();
    }

    inline bool operator==(Float const &a, float b) noexcept {
        return a.GetValue() == b;
    }

    inline bool operator==(float a, Float const &b) noexcept {
        return a == b.GetValue();
    }

    // !=
    inline bool operator!=(Float const &a, Float const &b) noexcept {
        return a.GetValue() != b.GetValue();
    }

    inline bool operator!=(Float const &a, float b) noexcept {
        return a.GetValue() != b;
    }

    inline bool operator!=(float a, Float const &b) noexcept {
        return a != b.GetValue();
    }

    // Sign of Float
    inline int Sign(Float const &v) noexcept {
        return static_cast<int>((0.f < v.GetValue()) - (v.GetValue() < 0.f));
    }

    // Sin of Float
    inline Float Sin(Float const &v) {
        return Float(default_tape.PushSingleNode(std::cos(v.GetValue()), v.NodeIndex()), std::sin(v.GetValue()));
    }

    // Cos of Float
    inline Float Cos(Float const &v) {
        return Float(default_tape.PushSingleNode(-std::sin(v.GetValue()), v.NodeIndex()), std::cos(v.GetValue()));
    }

    // Tan of Float
    inline Float Tan(Float const &v) {
        return Float(default_tape.PushSingleNode(2.f / (std::cos(2.f * v.GetValue()) + 1.f), v.NodeIndex()),
                     std::tan(v.GetValue()));
    }

    // Exp of Float
    inline Float Exp(Float const &v) {
        return Float(default_tape.PushSingleNode(std::exp(v.GetValue()), v.NodeIndex()), std::exp(v.GetValue()));
    }

    // Log of Float
    inline Float Log(Float const &v) {
        assert(v > 0.f);
        return Float(default_tape.PushSingleNode(1.f / v.GetValue(), v.NodeIndex()), std::log(v.GetValue()));
    }

    // Pow of Float
    inline Float Pow(Float const &v, float k) {
        return Float(default_tape.PushSingleNode(k * std::pow(v.GetValue(), k - 1.f), v.NodeIndex()),
                     std::pow(v.GetValue(), k));
    }

    // Sqrt of Float
    inline Float Sqrt(Float const &v) {
        assert(v != 0.f);
        return Float(default_tape.PushSingleNode(0.5f / std::sqrt(v.GetValue()), v.NodeIndex()),
                     std::sqrt(v.GetValue()));
    }

    // Abs of Float
    inline Float Abs(Float const &v) {
        assert(v != 0.f);
        return Float(default_tape.PushSingleNode(Sign(v), v.NodeIndex()), std::abs(v.GetValue()));
    }

    // Max of two Float
    inline Float Max(Float const &a, Float const &b) {
        return (a.GetValue() > b.GetValue() ? a : b);
    }

    // Min of two Float
    inline Float Min(Float const &a, Float const &b) {
        return (a.GetValue() < b.GetValue() ? a : b);
    }

} // drdemo namespace

#endif //DRDEMO_RAD_HPP
