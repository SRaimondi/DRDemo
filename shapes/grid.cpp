//
// Created by Simon on 19.06.2017.
//

#include "grid.hpp"

namespace drdemo {

    SignedDistanceGrid::SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b) {
        // Set number fo points along each dimension
        num_points[0] = n_x;
        num_points[1] = n_y;
        num_points[2] = n_z;
        total_points = n_x * n_y * n_z;
        // Set bounds
        bounds = b;
        // Compute voxel width
        Vector3f extent = bounds.Extent();
        for (int axis = 0; axis < 3; ++axis) {
            width[axis] = extent[axis] / static_cast<float>(num_points[axis]);
            inv_width[axis] = 1.f / width[axis];
        }
        // Allocate data
        data = new Float[total_points];
    }

    SignedDistanceGrid::SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b,
                                           float const *const raw_data) {
        // Set number fo points along each dimension
        num_points[0] = n_x;
        num_points[1] = n_y;
        num_points[2] = n_z;
        total_points = n_x * n_y * n_z;
        // Set bounds
        bounds = b;
        // Compute voxel width
        Vector3f extent = bounds.Extent();
        for (int axis = 0; axis < 3; ++axis) {
            width[axis] = extent[axis] / static_cast<float>(num_points[axis]);
            inv_width[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
        }
        // Allocate data
        data = new Float[total_points];
        // Copy values
        for (size_t i = 0; i < total_points; ++i) {
            data[i] = raw_data[i];
        }
    }

    SignedDistanceGrid::~SignedDistanceGrid() {
        if (data != nullptr) {
            delete[] data;
        }
    }

    bool SignedDistanceGrid::Intersect(Ray const &ray, Interaction *const interaction) const {
        return false;
    }

    bool SignedDistanceGrid::IntersectP(Ray const &ray) const {
        return false;
    }

    BBOX SignedDistanceGrid::BBox() const {
        return bounds;
    }

    Vector3f SignedDistanceGrid::Centroid() const {
        return Vector3f(0.f, 0.f, 0.f);
    }

    std::string SignedDistanceGrid::ToString() const {
        return "Signed Distance Grid with dimension: " + std::to_string(num_points[0]) + ", " +
               std::to_string(num_points[1]) + ", " + std::to_string(num_points[2]);
    }

    void SignedDistanceGrid::GetDiffVariables(std::vector<Float const *> &vars) const {
        for (size_t i = 0; i < total_points; ++i) {
            vars.push_back(&data[i]);
        }
    }

    size_t SignedDistanceGrid::GetNumVars() const noexcept {
        return total_points;
    }

    void SignedDistanceGrid::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {
        size_t used_vars = 0;
        for (size_t i = 0; i < total_points; ++i) {
            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
        }
    }

} // drdemo namespace
