//
// Created by simon on 17.04.17.
//

#ifndef DRDEMO_PROMOTE_HPP
#define DRDEMO_PROMOTE_HPP

/**
 * This file contains the definition of a template that is responsible to do the correct type promotion when multiple
 * types are involved in some mathematical operations for example
 */

#include "rad.hpp"
#include "fad.hpp"
#include "ifthenelse.hpp"

namespace traits {

    template<typename Ta, typename Tb>
    class Promote {
    public:
        // In general, return the "biggest" type
        typedef typename IfThenElse<(sizeof(Ta) > sizeof(Tb)), Ta, Tb>::TResult TResult;
    };

    // Specialization for the same type
    template<typename T>
    class Promote<T, T> {
    public:
        typedef T TResult;
    };

    // Specialization for some common types
    template<>
    class Promote<float, double> {
    public:
        typedef double TResult;
    };

    template<>
    class Promote<double, float> {
    public:
        typedef double TResult;
    };

    template<typename T>
    class Promote<rad::Variable<T>, T> {
    public:
        typedef typename rad::Variable<T> TResult;
    };

    template<typename T>
    class Promote<T, rad::Variable<T>> {
    public:
        typedef typename rad::Variable<T> TResult;
    };

    template<typename T>
    class Promote<fad::Variable<T>, T> {
    public:
        typedef typename fad::Variable<T> TResult;
    };

    template<typename T>
    class Promote<T, fad::Variable<T>> {
    public:
        typedef typename fad::Variable<T> TResult;
    };

} // traits namespace

#endif //DRDEMO_PROMOTE_HPP
