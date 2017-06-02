//
// Created by simon on 05.05.17.
//

#include "rad.hpp"

namespace drdemo {

    // Initialize default tape
    Tape default_tape = Tape();

    TapeNode::TapeNode(float w1, size_t p1, float w2, size_t p2) noexcept {
        assert(!std::isnan(w1));
        assert(!std::isnan(w2));
        assert(!std::isinf(w1));
        assert(!std::isinf(w2));
        // Set first child
        weights[0] = w1;
        parent_i[0] = p1;
        // Set second child
        weights[1] = w2;
        parent_i[1] = p2;
    }

    Tape::Tape(size_t starting_size) : nodes(starting_size) /* : enabled(true) */ {}

    void Tape::Clear(size_t starting_index) {
        nodes.Cut(starting_index);
    }

    size_t Tape::PushLeaf() {
        // if (IsEnabled()) {
        size_t const index = Size();
        nodes.Append(TapeNode(0.f, index, 0.f, index));
        return index;
        // }
        // return 0;
    }

    size_t Tape::PushSingleNode(float w, size_t p) {
        // if (IsEnabled()) {
        size_t const index = Size();
        nodes.Append(TapeNode(w, p, 0.f, index));
        return index;
        // }
        // return 0;
    }

    size_t Tape::PushTwoNode(float w1, size_t p1, float w2, size_t p2) {
        // if (IsEnabled()) {
        size_t const index = Size();
        nodes.Append(TapeNode(w1, p1, w2, p2));
        return index;
        // }
        // return 0;
    }

    Float::Float(float v)
    // Set the value of the variable and push it on the default_tape
            : value(v), node_index(default_tape.PushLeaf()) {}

    Float::Float(size_t index, float v) noexcept
            : value(v), node_index(index) {}

    Float::Float(Float const &other) noexcept
            : value(other.value), node_index(other.node_index) {}

    Float &Float::operator=(Float const &other) {
        if (this != &other) {
            value = other.value;
#ifdef FLOAT_NO_ALIAS
            node_index = default_tape.PushLeaf();
#else
            node_index = other.node_index;
#endif
        }
        return *this;
    }

    Float &Float::operator=(float v) noexcept {
        value = v;
        return *this;
    }

    std::ostream &operator<<(std::ostream &os, Float const &var) {
        os << var.GetValue();
        return os;
    }

} // drdemo namespace