//
// Created by simon on 13.04.17.
//

#ifndef DRDEMO_RAD_HPP
#define DRDEMO_RAD_HPP

/**
 * Reverse automatic differentiation tools file
 */

#include <vector>
#include <cstring>
#include <cassert>
#include <cmath>
#include <map>
#include "common.hpp"

namespace rad {

    /**
     * The TapeNode struct stores the information used in the construction of the Tape that is the used for the reverse
     * automatic differentiation process
     */
    template<typename T>
    struct TapeNode {
        // Weight associated with the respective nodes
        T weights[2];
        // Parent node indices in the tape
        size_t parent_i[2];

        TapeNode(T w1, size_t p1, T w2, size_t p2) noexcept;
    };

    template<typename T>
    TapeNode<T>::TapeNode(T w1, size_t p1, T w2, size_t p2) noexcept {
        // Set first child
        weights[0] = w1;
        parent_i[0] = p1;
        // Set second child
        weights[1] = w2;
        parent_i[1] = p2;
    }

    // Forward declare the Tape class, it's used in the Variable class
    template<typename T>
    class Tape;

    /**
     * The Variable class holds the information necessary for the evaluation of the derivatives of an algorithm using reverse
     * mode automatic differentiation.
     * It must be always associated with a Tape class in order to work
     */
    template<typename T>
    class Variable {
    private:
        // Pointer to the Tape associated with the Variable
        Tape<T> *tape;
        // TapeNode index in the Tape
        size_t node_index;
        // Value of the variable
        T value;

    public:
        Variable() noexcept;

        template<typename E>
        Variable(Tape<T> *const t, size_t n_i, E const &v) noexcept;

        // Assignment operator with another variable
        template<typename E>
        Variable<T> &operator=(Variable<E> const &var) noexcept;

        // Value assignment
        template<typename E>
        Variable<T> &operator=(E const &v) noexcept;

        // Access Variable fields
        inline Tape<T> *AssociatedTape() const noexcept {
            return tape;
        }

        inline size_t NodeIndex() const noexcept {
            return node_index;
        }

        inline T const &Value() const noexcept {
            return value;
        }

        // Math operators on itself
        inline Variable<T> &operator+=(Variable<T> const &a) {
            *this = *this + a;
            return *this;
        }

        inline Variable<T> &operator+=(T const &a) {
            *this = *this + a;
            return *this;
        }

        inline Variable<T> &operator-=(Variable<T> const &a) {
            *this = *this - a;
            return *this;
        }

        inline Variable<T> &operator-=(T const &a) {
            *this = *this - a;
            return *this;
        }

        inline Variable<T> &operator*=(Variable<T> const &a) {
            *this = *this * a;
            return *this;
        }

        inline Variable<T> &operator*=(T const &a) {
            *this = *this * a;
            return *this;
        }

        inline Variable<T> &operator/=(Variable<T> const &a) {
            *this = *this / a;
            return *this;
        }

        inline Variable<T> &operator/=(T const &a) {
            *this = *this / a;
            return *this;
        }
    };

    template<typename T>
    Variable<T>::Variable() noexcept
            : tape(nullptr), node_index(0), value(T(0)) {}

    template<typename T>
    template<typename E>
    Variable<T>::Variable(Tape<T> *const t, size_t n_i, E const &v) noexcept
            : tape(t), node_index(n_i), value(v) {}

    template<typename T>
    template<typename E>
    Variable<T> &Variable<T>::operator=(Variable<E> const &var) noexcept {
        // The assignment operator creates an alias by default
        if (this != &var) {
            this->tape = var.AssociatedTape();
            this->node_index = var.NodeIndex();
            this->value = var.Value();
        }

        return *this;
    }

    template<typename T>
    template<typename E>
    Variable<T> &Variable<T>::operator=(E const &v) noexcept {
        // Assign only a new value to the variable
        this->value = v;

        return *this;
    }

    /**
     * Declare the mathematical operators and function that allow to write expression using Variable<T>/Variable<T> and
     * Variable<T> / T
     */

    // Sum of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator+(Variable<T> const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr &&
               b.AssociatedTape() != nullptr &&
               a.AssociatedTape() == b.AssociatedTape());
#endif

        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushTwoNode(T(1), a.NodeIndex(), T(1), b.NodeIndex()),
                           a.Value() + b.Value());
    }

    // Sum of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator+(Variable<T> const &a, T const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr);
#endif
        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushSingleNode(T(1), a.NodeIndex()),
                           a.Value() + b);
    }

    template<typename T>
    inline Variable<T>
    operator+(T const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(b.AssociatedTape() != nullptr);
#endif
        return Variable<T>(b.AssociatedTape(),
                           b.AssociatedTape()->PushSingleNode(T(1), b.NodeIndex()),
                           a + b.Value());
    }

    // Subtraction of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator-(Variable<T> const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr &&
               b.AssociatedTape() != nullptr &&
               a.AssociatedTape() == b.AssociatedTape());
#endif

        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushTwoNode(T(1), a.NodeIndex(), T(-1), b.NodeIndex()),
                           a.Value() - b.Value());
    }

    // Subtraction of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator-(Variable<T> const &a, T const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr);
#endif
        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushSingleNode(T(1), a.NodeIndex()),
                           a.Value() - b);
    }

    template<typename T>
    inline Variable<T>
    operator-(T const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(b.AssociatedTape() != nullptr);
#endif
        return Variable<T>(b.AssociatedTape(),
                           b.AssociatedTape()->PushSingleNode(T(-1), b.NodeIndex()),
                           a - b.Value());
    }

    // Multiplication of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator*(Variable<T> const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr &&
               b.AssociatedTape() != nullptr &&
               a.AssociatedTape() == b.AssociatedTape());
#endif

        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushTwoNode(b.Value(), a.NodeIndex(), a.Value(), b.NodeIndex()),
                           a.Value() * b.Value());
    }

    // Multiplication of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator*(Variable<T> const &a, T const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr);
#endif
        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushSingleNode(b, a.NodeIndex()),
                           a.Value() * b);
    }

    template<typename T>
    inline Variable<T>
    operator*(T const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(b.AssociatedTape() != nullptr);
#endif
        return Variable<T>(b.AssociatedTape(),
                           b.AssociatedTape()->PushSingleNode(a, b.NodeIndex()),
                           a + b.Value());
    }

    // Division of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator/(Variable<T> const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr &&
               b.AssociatedTape() != nullptr &&
               a.AssociatedTape() == b.AssociatedTape());
#endif
        // Select tape to use
        auto tape = a.AssociatedTape() != nullptr ? a.AssociatedTape() : b.AssociatedTape();

        // If f(a,b) = a/b then df/da = 1/b and df/db = -a/(b*b)
        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushTwoNode(T(1) / b.Value(), a.NodeIndex(),
                                                           -a.Value() / (b.Value() * b.Value()), b.NodeIndex()),
                           a.Value() / b.Value());
    }

    // Division of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator/(Variable<T> const &a, T const &b) {
#ifdef DEBUG
        assert(a.AssociatedTape() != nullptr);
#endif
        return Variable<T>(a.AssociatedTape(),
                           a.AssociatedTape()->PushSingleNode(T(1) / b, a.NodeIndex()),
                           a.Value() / b);
    }

    template<typename T>
    inline Variable<T>
    operator/(T const &a, Variable<T> const &b) {
#ifdef DEBUG
        assert(b.AssociatedTape() != nullptr);
#endif
        return Variable<T>(b.AssociatedTape(),
                           b.AssociatedTape()->PushSingleNode(-a / (b.Value() * b.Value()), b.NodeIndex()),
                           a / b.Value());
    }

    // Comparison operators
    template<typename T>
    inline bool
    operator<(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() < b.Value();
    }

    template<typename T>
    inline bool
    operator<=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() <= b.Value();
    }

    template<typename T>
    inline bool
    operator>(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() > b.Value();
    }

    template<typename T>
    inline bool
    operator>=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() >= b.Value();
    }

    template<typename T>
    inline bool
    operator==(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() == b.Value();
    }

    template<typename T>
    inline bool
    operator!=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Value() != b.Value();
    }

    // Sign of Variable<T>
    template<typename T>
    inline int
    Sign(Variable<T> const &var) {
        return static_cast<int>((T(0) < var.Value()) - (var.Value() < T(0)));
    }

    // Sin of Variable<T>
    template<typename T>
    inline Variable<T>
    Sin(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(std::cos(var.Value()), var.NodeIndex()),
                           std::sin(var.Value()));
    }

    // Cos of Variable<T>
    template<typename T>
    inline Variable<T>
    Cos(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(std::sin(var.Value()), var.NodeIndex()),
                           std::cos(var.Value()));
    }

    // Exp of Variable<T>
    template<typename T>
    inline Variable<T>
    Exp(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(std::exp(var.Value()), var.NodeIndex()),
                           std::exp(var.Value()));
    }

    // Log of Variable<T>
    template<typename T>
    inline Variable<T>
    Log(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
        assert(var.Value() > T(0));
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(T(1) / var.Value(), var.NodeIndex()),
                           std::log(var.Value()));
    }

    // Pow of Variable<T>
    template<typename T>
    inline Variable<T>
    Pow(Variable<T> const &var, float k) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(k * std::pow(var.Value(), k - 1.f), var.NodeIndex()),
                           std::pow(var.Value(), k));
    }

    // Sqrt of Variable<T>
    template<typename T>
    inline Variable<T>
    Sqrt(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
        assert(var.Value() >= T(0));
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(T(0.5) / std::sqrt(var.Value()), var.NodeIndex()),
                           std::sqrt(var.Value()));
    }

    // Abs of Variable<T>
    template<typename T>
    inline Variable<T>
    Abs(Variable<T> const &var) {
#ifdef DEBUG
        assert(var != T(0));
#endif
        return Variable<T>(var.AssociatedTape(),
                           var.AssociatedTape()->PushSingleNode(utils::Sign(var.Value()), var.NodeIndex()),
                           std::abs(var.Value()));
    }

    /**
     * The Tape class is used to store information about the differentiation process during the execution of our program
     */
    template<typename T>
    class Tape {
    private:
        // Storage of the TapeNode's on the tape
        std::vector<TapeNode<T> > nodes;

    public:
        // Create Tape, optional argument is an already given size of the tape
        Tape(size_t size = 10) noexcept;

        // Request a new Variable from the Tape, optional argument is the value of it, otherwise it is set to zero
        Variable<T> NewVariable(T v = T(0)) noexcept;

        // Access tape node at given index
        inline TapeNode<T> const &At(size_t i) const {
#ifdef DEBUG
            return nodes.at(i);
#else
            return nodes[i];
#endif
        }

        // Access size of the Tape
        inline size_t Size() const {
            return nodes.size();
        }

        // Clear Tape, from the beginning or given a starting index
        void Clear(size_t start_index = 0);

        // Push a Zero value node (leaf node), returns the index of the node on the Tape
        size_t PushLeaf();

        // Push a node that depends only on another single node given the value of the node and the parent index
        size_t PushSingleNode(T w, size_t p);

        // Push a node that depends on two children given both the values and the parents indices
        size_t PushTwoNode(T w1, size_t p1, T w2, size_t p2);
    };

    template<typename T>
    Tape<T>::Tape(size_t size) noexcept
            : nodes() {
        nodes.reserve(size);
    }

    template<typename T>
    Variable<T> Tape<T>::NewVariable(T v) noexcept {
        return Variable<T>(this, PushLeaf(), v);
    }

    template<typename T>
    void Tape<T>::Clear(size_t start_index) {
        nodes.erase(nodes.begin() + start_index, nodes.end());
    }

    template<typename T>
    size_t Tape<T>::PushLeaf() {
        size_t const index = Size();
        nodes.push_back(TapeNode<T>(T(0), index, T(0), index));
        return index;
    }

    template<typename T>
    size_t Tape<T>::PushSingleNode(T w, size_t p) {
        size_t const index = Size();
        nodes.push_back(TapeNode<T>(w, p, T(0), index));
        return index;
    }

    template<typename T>
    size_t Tape<T>::PushTwoNode(T w1, size_t p1, T w2, size_t p2) {
        size_t const index = Size();
        nodes.push_back(TapeNode<T>(w1, p1, w2, p2));
        return index;
    }

    /**
     * The Derivatives class is used in conjugation with a variable to compute the derivatives
     */
    template<typename T>
    class Derivatives {
    private:
        // Associate variables with derivatives
        std::map<Variable<T>, std::vector<T> > var_derivatives_map;

    public:
        Derivatives();

        // Compute derivatives for a given variable
        void ComputeDerivatives(Variable<T> const &var);

        // Request derivative for a given outgoing variable with respect to a given input variable
        // Basically for df/dx, var_out is f and var_in x
        T D_Wrt(Variable<T> const &var_out, Variable<T> const &var_in) const;
    };

    template<typename T>
    Derivatives<T>::Derivatives()
            : var_derivatives_map() {}

    template<typename T>
    void Derivatives<T>::ComputeDerivatives(Variable<T> const &var) {
#ifdef DEBUG
        assert(var.AssociatedTape() != nullptr);
#endif
        size_t const size = var.AssociatedTape()->Size();
        // Check if element is already in the map
        if (var_derivatives_map.find(var) != var_derivatives_map.end()) {
            std::cout << "Variable derivatives already computed, exiting..." << std::endl;
            return;
        }
        // Create new entry
        var_derivatives_map[var] = std::vector<T>(size, T(0));
        // Get reference to make code more clear
        std::vector<T> &derivs = var_derivatives_map[var];

        // Seed at index of the input variable
        derivs[var.NodeIndex()] = T(1);

        // Traverse the tape in reverse
        for (size_t i = 0; i < size; i++) {
            size_t const rev_index = size - 1 - i;
            TapeNode<T> const &node = var.AssociatedTape()->At(rev_index);
            // Add children contribution
            derivs[node.parent_i[0]] += node.weights[0] * derivs[rev_index];
            derivs[node.parent_i[1]] += node.weights[1] * derivs[rev_index];
        }
    }

    template<typename T>
    T Derivatives<T>::D_Wrt(Variable<T> const &var_out, Variable<T> const &var_in) const {
        // Check if we computed the derivatives for the given out variable
        auto it = var_derivatives_map.find(var_out);
        if (it == var_derivatives_map.end()) {
            std::cout << "Could not find given output variable!" << std::endl;
            return T(0);
        } else {
            return it->second[var_in.NodeIndex()];
        }
    }

} // rad namespace

#endif //DRDEMO_RAD_HPP
