cmake_minimum_required(VERSION 3.7)
project(DRDemo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set DEBUG build mode flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")

# Include directories
include_directories(
        #       old/autodiff old/utils old/raytracing
        rad core film)

set(SOURCE_FILES
        #        old/main.cpp
        #        autodiff/rad.hpp
        #        autodiff/fad.hpp
        #        utils/vector.hpp
        #        utils/common.hpp
        #        utils/common.hpp
        #        utils/ifthenelse.hpp
        #        utils/promote.hpp
        #        autodiff/diffobject.hpp
        #        raytracing/ray.hpp
        #        raytracing/image.hpp
        #        raytracing/shape.hpp
        #        raytracing/intersection.hpp
        #        raytracing/sphere.hpp
        #        raytracing/sphere.cpp
        main.cpp
        rad/rad.cpp
        rad/rad.hpp
        rad/derivative.cpp
        rad/derivative.hpp
        core/geometry.hpp
        core/spectrum.hpp
        core/camera.hpp
        core/film.hpp
        core/film.cpp
        core/interaction.hpp
        core/shape.hpp
        core/scene.hpp
        core/scene.cpp
        core/light.hpp film/box_film.cpp film/box_film.hpp)

add_executable(DRDemo ${SOURCE_FILES})