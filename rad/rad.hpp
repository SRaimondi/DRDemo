//
// Created by simon on 05.05.17.
//

#ifndef DRDEMO_RAD_HPP
#define DRDEMO_RAD_HPP

/**
 * Reverse automatic differentiation tools
 */

#include <cstdlib>
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

    // Declare static Tape variable
    static Tape default_tape = Tape();

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
        // List of different constructors
        Float(float v = 0.f) noexcept;

        // Construct Float given index on tape and value
        Float(size_t index, float v) noexcept;

        Float(Float const &other) = default;

        // Assignment operator, default is to create an alias of the variable, not pushing a new one on the stack
        Float &operator=(Float const &other) = default;

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


} // drdemo namespace

#endif //DRDEMO_RAD_HPP
