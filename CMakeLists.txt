cmake_minimum_required(VERSION 3.7)
project(DRDemo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set DEBUG build mode flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")

# Include directories
include_directories(rad core film tonemapper shapes camera)

set(SOURCE_FILES
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
        core/light.hpp
        film/box_film.cpp
        film/box_film.hpp
        core/tonemapper.hpp
        tonemapper/clamp_tonemapper.cpp
        tonemapper/clamp_tonemapper.hpp
        shapes/sphere.cpp
        shapes/sphere.hpp
        core/common.hpp
        core/renderer.hpp
        camera/pinhole_camera.cpp
        camera/pinhole_camera.hpp)

add_executable(DRDemo ${SOURCE_FILES})