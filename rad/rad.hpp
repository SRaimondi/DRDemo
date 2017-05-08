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

    public:
        // Tape default constructor
        Tape();

        // Access TapeNode at given index
        inline TapeNode const &At(size_t index) const {
            return nodes.at(index);
        }

        // Get size of the Tape
        inline size_t Size() const { return nodes.size(); }

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
        explicit Float(float v = 0.f) noexcept;

        // Construct Float given index on tape and value
        Float(size_t index, float v) noexcept;

        Float(Float const &other) noexcept;

        // Assignment operator, default is to create an alias of the variable, not pushing a new one on the stack
        Float &operator=(Float const &other) = default;

        // Set Float value
        Float &operator=(float v) noexcept;

        // Convert to built-in float type
        inline float Value() const noexcept { return value; }

        // Get node index
        inline size_t NodeIndex() const noexcept { return node_index; }
    };

    /**
     * Declare all the Float mathematical operators
     */

    // Negation
    inline Float operator-(Float const &v) {
        return Float(default_tape.PushSingleNode(-1.f, v.NodeIndex()), -v.Value());
    }

    // Sum of Float and Float
    inline Float operator+(Float const &a, Float const &b) {
        return Float(default_tape.PushTwoNode(1.f, a.NodeIndex(), 1.f, b.NodeIndex()), a.Value() + b.Value());
    }

    // Sum of Float and float
    inline Float operator+(Float const &a, float b) {
        return Float(default_tape.PushSingleNode(1.f, a.NodeIndex()), a.Value() + b);
    }

    // Sum of float and Float
    inline Float operator+(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(1.f, b.NodeIndex()), a + b.Value());
    }

    // Subtraction of Float and Float
    inline Float operator-(Float const &a, Float const &b) {
        return Float(default_tape.PushTwoNode(1.f, a.NodeIndex(), -1.f, b.NodeIndex()), a.Value() - b.Value());
    }

    // Subtraction of Float and float
    inline Float operator-(Float const &a, float b) {
        return Float(default_tape.PushSingleNode(1.f, a.NodeIndex()), a.Value() - b);
    }

    // Subtraction of float and Float
    inline Float operator-(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(-1.f, b.NodeIndex()), a - b.Value());
    }

    // Multiplication of Float and Float
    inline Float operator*(Float const &a, Float const &b) {
        return Float(default_tape.PushTwoNode(b.Value(), a.NodeIndex(), a.Value(), b.NodeIndex()),
                     a.Value() * b.Value());
    }

    // Multiplication of Float and float
    inline Float operator*(Float const &a, float b) {
        return Float(default_tape.PushSingleNode(b, a.NodeIndex()), a.Value() * b);
    }

    // Multiplication of float and Float
    inline Float operator*(float a, Float const &b) {
        return Float(default_tape.PushSingleNode(a, b.NodeIndex()), a * b.Value());
    }

    // Division of Float and Float
    inline Float operator/(Float const &a, Float const &b) {
        // If f(a,b) = a/b, then df/da = 1/b and df/db = -a/(b*b)
        return Float(default_tape.PushTwoNode(1.f / b.Value(), a.NodeIndex(),
                                              -a.Value() / (b.Value() * b.Value()), b.NodeIndex()),
                     a.Value() / b.Value());
    }

    // Division of Float and float
    inline Float operator/(Float const &a, float b) {
        return Float(default_tape.PushSingleNode(1.f / b, a.NodeIndex()), a.Value() / b);
    }

    // Division of float and Float
    inline Float operator/(float a, Float const &b) {
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
    inline Float Sin(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(std::cos(v.Value()), v.NodeIndex()), std::sin(v.Value()));
    }

    // Cos of Float
    inline Float Cos(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(-std::sin(v.Value()), v.NodeIndex()), std::cos(v.Value()));
    }

    // Tan of Float
    inline Float Tan(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(2.f / (std::cos(2.f * v.Value()) + 1.f), v.NodeIndex()),
                     std::tan(v.Value()));
    }

    // Exp of Float
    inline Float Exp(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(std::exp(v.Value()), v.NodeIndex()), std::exp(v.Value()));
    }

    // Log of Float
    inline Float Log(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(1.f / v.Value(), v.NodeIndex()), std::log(v.Value()));
    }

    // Pow of Float
    inline Float Pow(Float const &v, float k) noexcept {
        return Float(default_tape.PushSingleNode(k * std::pow(v.Value(), k - 1.f), v.NodeIndex()),
                     std::pow(v.Value(), k));
    }

    // Sqrt of Float
    inline Float Sqrt(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(0.5f / std::sqrt(v.Value()), v.NodeIndex()), std::sqrt(v.Value()));
    }

    // Abs of Float
    inline Float Abs(Float const &v) noexcept {
        return Float(default_tape.PushSingleNode(Sign(v), v.NodeIndex()), std::abs(v.Value()));
    }

} // drdemo namespace

#endif //DRDEMO_RAD_HPP
