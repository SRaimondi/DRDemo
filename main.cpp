#include <iostream>
#include "rad.hpp"
#include "vector.hpp"

int main(int argc, char *argv[]) {
    // Define custom float type
    using Float = rad::Variable<float>;

    // Create tape
    rad::Tape<float> tape;

    // Create variables
    Float x = tape.NewVariable(1.f);
    Float y = tape.NewVariable(2.f);

    // Compute result
    Float z = rad::Exp(x) * rad::Sin(y);

    Float t = tape.NewVariable(5.f);

    Float f = z * t * t;

    // Compute derivatives
    rad::Derivatives<float> derivatives;
    derivatives.ComputeDerivatives(z);
    std::cout << "dz/dx: " << derivatives.D_Wrt(z, x) << std::endl;
    std::cout << "dz/dy: " << derivatives.D_Wrt(z, y) << std::endl;
    std::cout << std::endl;

    derivatives.ComputeDerivatives(f);
    std::cout << "df/dx: " << derivatives.D_Wrt(f, x) << std::endl;
    std::cout << "df/dy: " << derivatives.D_Wrt(f, y) << std::endl;
    std::cout << "df/dt: " << derivatives.D_Wrt(f, t) << std::endl;

    return 0;
}