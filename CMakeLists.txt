cmake_minimum_required(VERSION 3.7)
project(DRDemo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set DEBUG build mode flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wpedantic -Wextra -DDEBUG")

# Include directories
include_directories(
        accelerators
        rad
        core
        film
        tonemapper
        shapes
        camera
        renderer
        integrator
        light
        utilities
        tests
        minimization)

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
        # shapes/sphere.cpp
        # shapes/sphere.hpp
        core/common.hpp
        core/renderer.hpp
        camera/pinhole_camera.cpp
        camera/pinhole_camera.hpp
        core/integrator.hpp
        renderer/simple_renderer.cpp
        renderer/simple_renderer.hpp
        # light/point_light.cpp
        # light/point_light.hpp
        integrator/direct_integrator.cpp
        integrator/direct_integrator.hpp
        rad/diff_object.hpp
        # light/directional_light.cpp
        # light/directional_light.hpp
        utilities/vector.hpp
        utilities/tape_storage.hpp
        core/bbox.hpp
        core/bbox.cpp
        utilities/iofile.cpp
        utilities/iofile.hpp
        shapes/triangle_mesh.cpp
        shapes/triangle_mesh.hpp
        accelerators/bvh.cpp
        accelerators/bvh.hpp
        shapes/grid.cpp
        shapes/grid.hpp
        core/lodepng.cpp
        core/lodepng.hpp
        minimization/scalar_function.hpp
        minimization/scalar_function.cpp
        # minimization/multi_view_energy.cpp
        # minimization/multi_view_energy.hpp
        minimization/gradient_descent.cpp
        minimization/gradient_descent.hpp
        minimization/reconstruction_energy.cpp
        minimization/reconstruction_energy.hpp
        # tests/geom_sphere.cpp
        # tests/geom_sphere.hpp
        tests/sdf_sphere.cpp
        tests/sdf_sphere.hpp
        light/ambient_light.cpp
        light/ambient_light.hpp
        tests/torus_test.cpp
        tests/torus_test.hpp
        minimization/reconstruction_energy_opt.cpp
        minimization/reconstruction_energy_opt.hpp
        tests/test_common.cpp tests/test_common.hpp
        tests/bunny_test.cpp
        tests/bunny_test.hpp
        tests/dragon_render.cpp
        tests/dragon_render.hpp
        tests/armadillo_render.cpp
        tests/armadillo_render.hpp
        tests/sdf_loading_render_test.cpp
        tests/sdf_loading_render_test.hpp
        tests/dragon_full_pipeline_test.cpp
        tests/dragon_full_pipeline_test.hpp
        camera/perspective_camera.cpp
        camera/perspective_camera.cpp
        camera/perspective_camera.hpp
        tests/dino_test.cpp
        tests/dino_test.hpp
        shapes/mac_grid.cpp
        shapes/mac_grid.hpp
        minimization/reconstruction_energy_light.cpp
        minimization/reconstruction_energy_light.hpp light/SH_light.cpp light/SH_light.hpp)

add_executable(DRDemo ${SOURCE_FILES})