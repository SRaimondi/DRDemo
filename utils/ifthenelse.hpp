//
// Created by simon on 17.04.17.
//

#ifndef DRDEMO_IFTHENELSE_HPP
#define DRDEMO_IFTHENELSE_HPP

/**
 * Utility template to do compute time if then else conditional statements
 */

namespace traits {

    // Primary template declaration
    template<bool C, typename Ta, typename Tb>
    class IfThenElse;

    // Partial specialization: true yields first type
    template<typename Ta, typename Tb>
    class IfThenElse<true, Ta, Tb> {
    public:
        typedef Ta TResult;
    };

    // Partial specialization: false yields second type
    template<typename Ta, typename Tb>
    class IfThenElse<false, Ta, Tb> {
    public:
        typedef Tb Tresult;
    };

} // traits namespace

#endif //DRDEMO_IFTHENELSE_HPP
