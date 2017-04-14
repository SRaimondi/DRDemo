cmake_minimum_required(VERSION 3.7)
project(DRDemo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set DEBUG build mode flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")

# Include directories
include_directories(autodiff utils)

set(SOURCE_FILES main.cpp autodiff/rad.hpp autodiff/fad.hpp utils/vector.hpp)
add_executable(DRDemo ${SOURCE_FILES})