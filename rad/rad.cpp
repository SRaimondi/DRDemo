//
// Created by simon on 05.05.17.
//

#include "rad.hpp"

namespace drdemo {

    TapeNode::TapeNode(float w1, size_t p1, float w2, size_t p2) noexcept {
        // Set first child
        weights[0] = w1;
        parent_i[0] = p1;
        // Set second child
        weights[1] = w2;
        parent_i[1] = p2;
    }

    Tape::Tape()
            : nodes() {}

    size_t Tape::PushLeaf() {
        size_t const index = Size();
        nodes.push_back(TapeNode(0.f, index, 0.f, index));
        return index;
    }

    size_t Tape::PushSingleNode(float w, size_t p) {
        size_t const index = Size();
        nodes.push_back(TapeNode(w, p, 0.f, index));
        return index;
    }

    size_t Tape::PushTwoNode(float w1, size_t p1, float w2, size_t p2) {
        size_t const index = Size();
        nodes.push_back(TapeNode(w1, p1, w2, p2));
        return index;
    }

    Float::Float(float v) noexcept
    // Set the value of the variable and push it on the default_tape
            : value(v), node_index(default_tape.PushLeaf()) {}

    Float::Float(size_t index, float v) noexcept
            : value(v), node_index(index) {}

} // drdemo namespace