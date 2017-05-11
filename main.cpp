//
// Created by simon on 08.05.17.
//

//#include <iostream>
//#include "derivative.hpp"
//#include "geometry.hpp"

#include "clamp_tonemapper.hpp"
#include "box_film.hpp"

#define WIDTH 512
#define HEIGHT 512

int main(void) {

    // Set namespace used
    using namespace drdemo;

//    Float x(1.f);
//    Float y(2.f);
//
//    // Compute result
//    Float z = Exp(x) * Sin(y);
//    Float t(5.f);
//    Float f = t * t * z;
//
//    // Compute derivatives
//    Derivatives derivatives;
//    derivatives.ComputeDerivatives(z);
//    std::cout << "dz/dx: " << derivatives.Dwrt(z, x) << std::endl;
//    std::cout << "dz/dy: " << derivatives.Dwrt(z, y) << std::endl;
//    std::cout << std::endl;
//
//    derivatives.ComputeDerivatives(f);
//    std::cout << "df/dx: " << derivatives.Dwrt(f, x) << std::endl;
//    std::cout << "df/dy: " << derivatives.Dwrt(f, y) << std::endl;
//    std::cout << "df/dt: " << derivatives.Dwrt(f, t) << std::endl;
//    std::cout << std::endl;

    // Create film
    BoxFilterFilm film(WIDTH, HEIGHT);
    // Set all image to red
    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < HEIGHT; j++) {
            film.AddSample(Spectrum(1.f, 0.f, 0.5f), i, j, 0.5f, 0.5f);
        }
    }

    // Create image
    ClampTonemapper tonemapper;

    tonemapper.Process(std::string("test.ppm"), film);

    return EXIT_SUCCESS;
}
