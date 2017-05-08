//
// Created by simon on 08.05.17.
//

#include <iostream>
#include "derivative.hpp"

int main(void) {

    // Set namespace used
    using namespace drdemo;

    Float x(1.f), y(2.f);

    // Compute result
    Float z = Exp(x) * Sin(y);
    Float t(5.f);
    Float f = z * t * t;

    // Compute derivatives
    Derivatives derivatives;
    derivatives.ComputeDerivatives(z);
    std::cout << "dz/dx: " << derivatives.DfwrtDx(z, x) << std::endl;
    std::cout << "dz/dy: " << derivatives.DfwrtDx(z, y) << std::endl;
    std::cout << std::endl;

    derivatives.ComputeDerivatives(f);
    std::cout << "df/dx: " << derivatives.DfwrtDx(f, x) << std::endl;
    std::cout << "df/dy: " << derivatives.DfwrtDx(f, y) << std::endl;
    std::cout << "df/dt: " << derivatives.DfwrtDx(f, t) << std::endl;
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
