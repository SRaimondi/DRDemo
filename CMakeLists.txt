cmake_minimum_required(VERSION 3.7)
project(DRDemo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set DEBUG build mode flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")

# Include directories
include_directories(old/autodiff old/utils old/raytracing)

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
        )
add_executable(DRDemo ${SOURCE_FILES})