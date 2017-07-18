//
// Created by Simon on 19.06.2017.
//

#include "grid.hpp"

namespace drdemo {

    void SignedDistanceGrid::PointsIndicesFromVoxel(size_t x, size_t y, size_t z, size_t *const indices) const {
        // Compute back face indices
        indices[0] = x + y * num_points[0] + z * num_points[0] * num_points[1];
        indices[1] = indices[0] + 1;
        indices[2] = x + (y + 1) * num_points[0] + z * num_points[0] * num_points[1];
        indices[3] = indices[2] + 1;
        // Compute front face indices
        indices[4] = x + y * num_points[0] + (z + 1) * num_points[0] * num_points[1];
        indices[5] = indices[4] + 1;
        indices[6] = x + (y + 1) * num_points[0] + (z + 1) * num_points[0] * num_points[1];
        indices[7] = indices[6] + 1;
    }

    Float SignedDistanceGrid::ValueAt(const Vector3F &p) const {
        // Convert position
        const Vector3f p_f = Tofloat(p);
        // Get voxel indices
        size_t voxel_i[3];
        for (int i = 0; i < 3; i++) { voxel_i[i] = PosToVoxel(p_f, i); }
        // Get point indices for the given voxel
        size_t indices[8];
        PointsIndicesFromVoxel(voxel_i[0], voxel_i[1], voxel_i[2], indices);

        // Compute minimum point of voxel
        const Vector3f voxel_min(bounds.MinPoint().x + voxel_i[0] * width.x,
                                 bounds.MinPoint().y + voxel_i[1] * width.y,
                                 bounds.MinPoint().z + voxel_i[2] * width.z);

        // Linear interpolate along x axis the eight values
        const Float tx = (p.x - voxel_min.x) * inv_width.x;
        const Float c01 = (1.f - tx) * data[indices[0]] + tx * data[indices[1]];
        const Float c23 = (1.f - tx) * data[indices[2]] + tx * data[indices[3]];
        const Float c45 = (1.f - tx) * data[indices[4]] + tx * data[indices[5]];
        const Float c67 = (1.f - tx) * data[indices[6]] + tx * data[indices[7]];

        // Linear interpolate along the y axis
        const Float ty = (p.y - voxel_min.y) * inv_width.y;
        const Float c0 = (1.f - ty) * c01 + ty * c23;
        const Float c1 = (1.f - ty) * c45 + ty * c67;

        // Return final value interpolated along z
        const Float tz = (p.z - voxel_min.z) * inv_width.z;

        return (1.f - tz) * c0 + tz * c1;
    }

    SignedDistanceGrid::SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b)
            : data(new Float[n_x * n_y * n_z]) {
        // Set number of points along each dimension
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
    }

    SignedDistanceGrid::SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b,
                                           float const *const raw_data)
            : data(new Float[n_x * n_y * n_z]) {
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
        // Copy values
        for (size_t i = 0; i < total_points; ++i) {
            data[i] = raw_data[i];
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
        return static_cast<size_t>(total_points);
    }

    void SignedDistanceGrid::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {
        // FIXME Initial simple attempt to propagate gradient directly into the data
        size_t used_vars = 0;
        for (size_t i = 0; i < total_points; ++i) {
            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
        }
    }

} // drdemo namespace