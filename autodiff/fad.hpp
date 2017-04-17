//
// Created by simon on 14.04.17.
//

#ifndef DRDEMO_FAD_HPP
#define DRDEMO_FAD_HPP

/**
 * Forward automatic differentiation tools file
 */

#include <cassert>
#include <cmath>
#include "common.hpp"

namespace FAD {

    /**
     * The Variable class implements DualNumber math to compute first order derivatives from expressions
     */
    template<typename T>
    class Variable {
    private:
        // Real and dual part
        T real;
        T dual;

    public:
        explicit Variable(T r = T(0), T d = T(0)) noexcept;

        // Access real and dual part
        inline T Real() const noexcept {
            return real;
        }

        inline T Dual() const noexcept {
            return dual;
        }

        // Set derivative with respect to this variable
        inline void SetDiff() noexcept {
            dual = T(1);
        }

        inline void UnsetDiff() noexcept {
            dual = T(0);
        }

        // Math operators on itself
        // Operators on self
        inline Variable<T> &operator+=(Variable const &a) noexcept {
            real += a.real;
            dual += a.dual;
            return *this;
        }

        inline Variable<T> &operator+=(T a) noexcept {
            real += a;
            return *this;
        }

        inline Variable<T> &operator-=(Variable<T> const &a) noexcept {
            real -= a.real;
            dual -= a.dual;
            return *this;
        }

        inline Variable<T> &operator-=(T a) noexcept {
            real -= a;
            return *this;
        }

        inline Variable<T> &operator*=(Variable<T> const &a) noexcept {
            real *= a.real;
            dual = dual * a.real + real * a.dual;
            return *this;
        }

        inline Variable<T> &operator*=(T a) noexcept {
            real *= a;
            dual *= a;
            return *this;
        }

        inline Variable<T> &operator/=(Variable<T> const &a) noexcept {
#ifdef DEBUG
            assert(a.real != T(0));
#endif
            real /= a.real;
            dual = (dual * a.real - real * a.dual) / (a.real * a.real);
            return *this;
        }

        inline Variable<T> &operator/=(T a) noexcept {
#ifdef DEBUG
            assert(a != T(0));
#endif
            real /= a;
            dual /= a;
            return *this;
        }
    };

    template<typename T>
    Variable<T>::Variable(T r, T d) noexcept
            : real(r), dual(d) {}

    /**
     * Declare the mathematical operators and function that allow to write expression using Variable<T>/Variable<T> and
     * Variable<T> / T
     */

    // Sum of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator+(Variable<T> const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a.Real() + b.Real(), a.Dual() + b.Dual());
    }

    // Sum of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator+(Variable<T> const &a, T const &b) noexcept {
        return Variable<T>(a.Real() + b, a.Dual());
    }

    template<typename T>
    inline Variable<T>
    operator+(T const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a + b.Real(), b.Dual());
    }

    // Subtraction of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator-(Variable<T> const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a.Real() - b.Real(), a.Dual() - b.Dual());
    }

    // Subtraction of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator-(Variable<T> const &a, T const &b) noexcept {
        return Variable<T>(a.Real() - b, a.Dual());
    }

    template<typename T>
    inline Variable<T>
    operator-(T const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a - b.Real(), -b.Dual());
    }

    // Multiplication of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator*(Variable<T> const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a.Real() * b.Real(), a.Dual() * b.Real() + a.Real() * b.Dual());
    }

    // Multiplication of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T> operator*(Variable<T> const &a, T const &b) noexcept {
        return Variable<T>(a.Real() * b, a.Dual() * b);
    }

    template<typename T>
    inline Variable<T> operator*(T const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a * b.Real(), a * b.Dual());
    }

    // Division of Variable<T> and Variable<T>
    template<typename T>
    inline Variable<T>
    operator/(Variable<T> const &a, Variable<T> const &b) noexcept {
#ifdef DEBUG
        assert(b.Real() != T(0));
#endif
        return Variable<T>(a.Real() / b.Real(), (a.Dual() * b.Real() - a.Real() * b.Dual()) / (b.Real() * b.Real()));
    }

    // Division of Variable<T> and T and vice-versa
    template<typename T>
    inline Variable<T>
    operator/(Variable<T> const &a, T const &b) noexcept {
        return Variable<T>(a.Real() / b, a.Dual() / b);
    }

    template<typename T>
    inline Variable<T>
    operator/(T const &a, Variable<T> const &b) noexcept {
        return Variable<T>(a / b.Real(), -a * b.Dual() / (b.Real() * b.Real()));
    }

    // Comparison operators
    template<typename T>
    inline bool
    operator<(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() < b.Real();
    }

    template<typename T>
    inline bool
    operator<=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() <= b.Real();
    }

    template<typename T>
    inline bool
    operator>(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() > b.Real();
    }

    template<typename T>
    inline bool
    operator>=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() >= b.Real();
    }

    template<typename T>
    inline bool
    operator==(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() == b.Real();
    }

    template<typename T>
    inline bool
    operator!=(Variable<T> const &a, Variable<T> const &b) noexcept {
        return a.Real() != b.Real();
    }

    // Sign of Variable<T>
    template<typename T>
    inline int
    Sign(Variable<T> const &var) {
        return static_cast<int>((T(0) < var.Real()) - (var.Real() < T(0)));
    }

    // Sin of Variable<T>
    template<typename T>
    inline Variable<T>
    Sin(Variable<T> const &a) noexcept {
        return Variable<T>(std::sin(a.Real()), a.Dual() * std::cos(a.Real()));
    }

    // Cos of Variable<T>
    template<typename T>
    inline Variable<T>
    Cos(Variable<T> const &a) noexcept {
        return Variable<T>(std::cos(a.Real()), -a.Dual() * std::sin(a.Real()));
    }

    // Tan of Variable<T>
    template<typename T>
    inline Variable<T>
    Tan(Variable<T> const &a) noexcept {
        return Variable<T>(std::tan(a.Real()), a.Dual() * T(2) / (std::cos(T(2) * a.Real()) + T(1)));
    }

    // Exp of Variable<T>
    template<typename T>
    inline Variable<T>
    Exp(Variable<T> const &a) noexcept {
        return Variable<T>(std::exp(a.Real()), a.Dual() * std::exp(a.Real()));
    }

    // Log of Variable<T>
    template<typename T>
    inline Variable<T>
    Log(Variable<T> const &a) noexcept {
#ifdef DEBUG
        assert(a.Real() > T(0));
#endif
        return Variable<T>(std::log(a.Real()), a.Dual() / a.Real());
    }

    // Power of Variable<T>
    template<typename T>
    inline Variable<T>
    Pow(Variable<T> const &a, float k) {
#ifdef DEBUG
        assert(a.Real() != T(0));
#endif
        return Variable<T>(std::pow(a.Real(), k), k * std::pow(a.Real(), k - T(1)) * a.Dual());
    }

    // Square of Variable<T>
    template<typename T>
    inline Variable<T>
    Sqrt(Variable<T> const &a) {
#ifdef DEBUG
        assert(a.Real() >= T(0));
#endif
        T sq = std::sqrt(a.Real());
        return Variable<T>(sq, a.Dual() * T(0.5) / sq);
    }

    // Absolute value of Variable<T>
    template<typename T>
    inline Variable<T>
    Abs(Variable<T> const &a) {
#ifdef DEBUG
        assert(a.Real() != T(0));
#endif
        return Variable<T>(std::abs(a.Real()), a.Dual() * utils::Sign(a.Real()));
    }

} // FAD namespace

#endif //DRDEMO_FAD_HPP
