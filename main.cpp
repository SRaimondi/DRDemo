#include <iostream>
#include "vector.hpp"

#define N 10

// Define spherical n-dimensional function
template<typename T, size_t SIZE>
T SphericalFunction(utils::Vector<T, SIZE> const &args) {
    return utils::Length2(args);
}

int main(/* int argc, char *argv[] */) {
//    // Define custom float type
//    using Float = rad::Variable<float>;
//    using FVector = utils::Vector<float, 10>;
//
//    FVector v1, v2;
//
//    for (size_t i = 0; i < 10; i++) {
//        v1[i] = i;
//        v2[i] = i * i;
//    }

    // Create tape
    rad::Tape<float> tape;

//    // Create variables
//    Float x = tape.NewVariable(1.f);
//    Float y = tape.NewVariable(2.f);
//
//    // Compute result
//    Float z = rad::Exp(x) * rad::Sin(y);
//    Float t = tape.NewVariable(5.f);
//    Float f = z * t * t;
//
//    // Compute derivatives
//    rad::Derivatives<float> derivatives;
//    derivatives.ComputeDerivatives(z);
//    std::cout << "dz/dx: " << derivatives.D_Wrt(z, x) << std::endl;
//    std::cout << "dz/dy: " << derivatives.D_Wrt(z, y) << std::endl;
//    std::cout << std::endl;
//
//    derivatives.ComputeDerivatives(f);
//    std::cout << "df/dx: " << derivatives.D_Wrt(f, x) << std::endl;
//    std::cout << "df/dy: " << derivatives.D_Wrt(f, y) << std::endl;
//    std::cout << "df/dt: " << derivatives.D_Wrt(f, t) << std::endl;
//    std::cout << std::endl;
//
//    FVector v3 = v1 * v1;
//    for (size_t i = 0; i < 10; i++) {
//        std::cout << v3[i] << std::endl;
//    }
//    std::cout << utils::Dot(v1, v1) << std::endl;

    // Create initial guess
    utils::Vector<rad::Variable<float>, N> x;
    for (size_t i = 0; i < N; i++) {
        x[i] = tape.NewVariable(i);
    }
    std::cout << "X:" << std::endl;
    for (size_t i = 0; i < N; i++) {
        std::cout << x[i].Value() << std::endl;
    }
    std::cout << std::endl;

    // Evaluate function
    rad::Variable<float> out = SphericalFunction(x);

    // Compute derivatives
    rad::Derivatives<float> derivatives;
    derivatives.ComputeDerivatives(out);

    // Compute gradient
    utils::Vector<float, N> gradient;
    for (size_t i = 0; i < N; i++) {
        gradient[i] = derivatives.D_Wrt(out, x[i]);
    }

    // Print gradient
    for (size_t i = 0; i < N; i++) {
        std::cout << gradient[i] << std::endl;
    }
    std::cout << std::endl;

    x = x - gradient;

    std::cout << "X after update" << std::endl;
    for (size_t i = 0; i < N; i++) {
        std::cout << x[i].Value() << std::endl;
    }

    return 0;
}