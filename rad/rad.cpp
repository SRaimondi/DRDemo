//
// Created by simon on 05.05.17.
//

#include "rad.hpp"

namespace drdemo {

    // Initialize default tape
    Tape default_tape = Tape();

    TapeNode::TapeNode(float w1, size_t p1, float w2, size_t p2) noexcept {
        // Set first child
        weights[0] = w1;
        parent_i[0] = p1;
        // Set second child
        weights[1] = w2;
        parent_i[1] = p2;
    }

    Tape::Tape() : enabled(true) {}

    void Tape::Clear(size_t starting_index) {
        nodes.erase(nodes.begin() + starting_index, nodes.end());
    }

    size_t Tape::PushLeaf() {
        if (IsEnabled()) {
            size_t const index = Size();
            nodes.push_back(TapeNode(0.f, index, 0.f, index));
            return index;
        }
        return 0;
    }

    size_t Tape::PushSingleNode(float w, size_t p) {
        if (IsEnabled()) {
            size_t const index = Size();
            nodes.push_back(TapeNode(w, p, 0.f, index));
            return index;
        }
        return 0;
    }

    size_t Tape::PushTwoNode(float w1, size_t p1, float w2, size_t p2) {
        if (IsEnabled()) {
            size_t const index = Size();
            nodes.push_back(TapeNode(w1, p1, w2, p2));
            return index;
        }
        return 0;
    }

    Float::Float(float v)
    // Set the value of the variable and push it on the default_tape
            : value(v), node_index(default_tape.PushLeaf()) {}

    Float::Float(size_t index, float v) noexcept
            : value(v), node_index(index) {}

    Float::Float(Float const &other) noexcept
            : value(other.value), node_index(other.node_index) {}

    Float &Float::operator=(float v) noexcept {
        value = v;
        return *this;
    }

    std::ostream &operator<<(std::ostream &os, Float const &var) {
        os << "Value: " << var.Value();
        return os;
    }

} // drdemo namespace