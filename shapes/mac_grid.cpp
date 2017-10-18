//
// Created by simon on 18.10.17.
//

#include "mac_grid.hpp"

namespace drdemo {

    MACGrid::MACGrid(int nx, int ny, int nz, const BBOX &b)
            : values(new Float[nx * ny * nz]) {
        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        total_voxels = nx * ny * nz;
        bounds = b;
        // Compute voxel size
        const Vector3f e = bounds.Extent();
        for (int axis = 0; axis < 3; ++axis) {
            v_width[axis] = e[axis] / (float) (dims[axis]);
            inv_v_width[axis] = (v_width[axis] == 0.f) ? 0.f : 1.f / v_width[axis];
        }
    }

    MACGrid::~MACGrid() {
        delete[] values;
    }

    Float MACGrid::SDFAt(const Vector3F &p) const {
        // Convert position to plain float for bbox
        const Vector3f p_float = Tofloat(p);
        // Check if we are outside of the grid
        if (!bounds.Inside(p_float)) {
            return Float(bounds.Distance(p_float) + 0.001f);
        }
        // Find voxel where the point is
        int v_i[3];
        for (int i = 0; i < 3; ++i) { v_i[i] = PosToVoxel(p_float, i); }
        // Compute the minimum and the maximum point of the voxel where we are in
        const Vector3F voxel_min(bounds.MinPoint().x + v_i[0] * v_width.x,
                                 bounds.MinPoint().y + v_i[1] * v_width.y,
                                 bounds.MinPoint().z + v_i[2] * v_width.z);
        // const Vector3f voxel_max = voxel_min + v_width;
        // Now look where we are in the voxel to see which other voxels we need to use to interpolate
        const Vector3F offset = (p - voxel_min) / v_width;
    }

    void MACGrid::Refine(const int *new_dims) {

    }

    bool MACGrid::Intersect(Ray const &ray, Interaction *interaction) const {
        return false;
    }

    bool MACGrid::IntersectP(Ray const &ray) const {
        return false;
    }

    BBOX MACGrid::BBox() const {
        return bounds;
    }

    Vector3f MACGrid::Centroid() const {
        return bounds.MinPoint() + 0.5f * bounds.Extent();
    }

    std::string MACGrid::ToString() const {
        return std::string("");
    }

    void MACGrid::GetDiffVariables(std::vector<Float const *> &vars) const {
        for (int i = 0; i < total_voxels; ++i) {
            vars.push_back(&values[i]);
        }
    }

    size_t MACGrid::GetNumVars() const noexcept {
        return static_cast<size_t>(total_voxels);
    }

    void MACGrid::UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) {
        size_t used_vars = 0;
        for (int i = 0; i < total_voxels; ++i) {
            values[i].SetValue(values[i].GetValue() + delta[starting_index + used_vars++]);
        }
    }

    void MACGrid::SetDiffVariables(const std::vector<float> &vals, size_t starting_index) {
        size_t used_vars = 0;
        for (int i = 0; i < total_voxels; ++i) {
            values[i].SetValue(vals[starting_index + used_vars++]);
        }
    }

} // drdemo namespace