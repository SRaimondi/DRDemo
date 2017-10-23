//
// Created by simon on 18.10.17.
//

#include <iofile.hpp>
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

    MACGrid::MACGrid(const std::string &sdf_file) {
        // Start by trying to reading the file
        std::vector<std::string> sdf_file_lines;
        if (ReadFile(sdf_file, sdf_file_lines)) {
            // No check on the format, we expect the file to be correct
            int n_x, n_y, n_z;
            // Get grid dimension from first line
            if (sscanf(sdf_file_lines[0].c_str(), "%d %d %d", &n_x, &n_y, &n_z) == 3) {
                // Set grid dimensions
                dims[0] = n_x;
                dims[1] = n_y;
                dims[2] = n_z;
                // Compute total number of points
                total_voxels = n_x * n_y * n_z;
                // Allocate memory for grid
                values = new Float[total_voxels];
            } else {
                std::cerr << "Error in first line of .sdf file" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Get grid size from third line
            float voxel_dim;
            if (sscanf(sdf_file_lines[2].c_str(), "%f", &voxel_dim) == 1) {
                // Set width and inv_width
                for (int i = 0; i < 3; i++) {
                    v_width[i] = voxel_dim;
                    inv_v_width[i] = 1.f / voxel_dim;
                }
            } else {
                std::cerr << "Error getting grid dimension" << std::endl;
                delete[] values;
                exit(EXIT_FAILURE);
            }

            // Get BBOX minimum point from second line
            float bbox_x, bbox_y, bbox_z;
            if (sscanf(sdf_file_lines[1].c_str(), "%f %f %f", &bbox_x, &bbox_y, &bbox_z) == 3) {
                const Vector3f bbox_min(bbox_x, bbox_y, bbox_z);
                const Vector3f bbox_dim(voxel_dim * dims[0], voxel_dim * dims[1],
                                        voxel_dim * dims[2]);
                bounds = BBOX(bbox_min, bbox_min + bbox_dim);
            } else {
                std::cerr << "Error reading bbox minimum point" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Read all the data and set the values
            float sdf_val;
            for (int i = 3; i < sdf_file_lines.size() - 1; ++i) {       // Minus 1 since last line is empty
                if (sscanf(sdf_file_lines[i].c_str(), "%f", &sdf_val) == 1) {
                    values[i - 3] = sdf_val;
                } else {
                    std::cerr << "Error reading sdf value" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            std::cerr << "Error trying to open SDF file!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    MACGrid::~MACGrid() {
        delete[] values;
    }

    Float MACGrid::SDFAt(const Vector3F &p) const {
        // Convert position to plain float for bbox
        const Vector3f p_float = Tofloat(p);
        // Compute bounds of "internal" grid
        const BBOX internal_bounds(bounds.MinPoint() + 0.5f * v_width, bounds.MaxPoint() - 0.5f * v_width);
        // Check if we are outside of the grid
        if (!internal_bounds.Inside(p_float)) {
            return Float(internal_bounds.Distance(p_float) + 0.001f);
        }
        // Find voxel where the point is
        // In this case we use the min point as the center of the min voxel, this makes computing the indices easier
        int v_i[3];
        for (int i = 0; i < 3; ++i) {
            v_i[i] = Clamp(static_cast<int>((p_float[i] - internal_bounds.MinPoint()[i]) * inv_v_width[i]), 0,
                           dims[i] - 2);
        }

        // Compute minimum point of voxel
        const Vector3f voxel_min(internal_bounds.MinPoint().x + v_i[0] * v_width.x,
                                 internal_bounds.MinPoint().y + v_i[1] * v_width.y,
                                 internal_bounds.MinPoint().z + v_i[2] * v_width.z);

        // Compute indices of the points of the voxel where we are in
        int indices[8];
        // Back face
        indices[0] = dims[1] * dims[0] * v_i[2] + dims[0] * v_i[1] + v_i[0];
        indices[1] = indices[0] + 1;
        indices[2] = dims[1] * dims[0] * v_i[2] + dims[0] * (v_i[1] + 1) + v_i[0];
        indices[3] = indices[2] + 1;
        // Front face
        indices[4] = dims[1] * dims[0] * (v_i[2] + 1) + dims[0] * v_i[1] + v_i[0];
        indices[5] = indices[4] + 1;
        indices[6] = dims[1] * dims[0] * (v_i[2] + 1) + dims[0] * (v_i[1] + 1) + v_i[0];
        indices[7] = indices[6] + 1;

        // Linear interpolate along x axis the eight values
        const Float tx = (p.x - voxel_min.x) * inv_v_width.x;
        const Float c01 = (1.f - tx) * values[indices[0]] + tx * values[indices[1]];
        const Float c23 = (1.f - tx) * values[indices[2]] + tx * values[indices[3]];
        const Float c45 = (1.f - tx) * values[indices[4]] + tx * values[indices[5]];
        const Float c67 = (1.f - tx) * values[indices[6]] + tx * values[indices[7]];

        // Linear interpolate along the y axis
        const Float ty = (p.y - voxel_min.y) * inv_v_width.y;
        const Float c0 = (1.f - ty) * c01 + ty * c23;
        const Float c1 = (1.f - ty) * c45 + ty * c67;

        // Return final value interpolated along z
        const Float tz = (p.z - voxel_min.z) * inv_v_width.z;

        return (1.f - tz) * c0 + tz * c1;
    }

    Vector3F MACGrid::NormalAt(const Vector3F &p) const {
        // Convert position to plain float
        const Vector3f p_float = Tofloat(p);
        // Find the voxel of the whole grid where we are
        int v_i[3];
        for (int i = 0; i < 3; ++i) { v_i[i] = PosToVoxel(p_float, i); }
        // Compute minimum point of voxel
        const Vector3f voxel_min(bounds.MinPoint().x + v_i[0] * v_width.x,
                                 bounds.MinPoint().y + v_i[1] * v_width.y,
                                 bounds.MinPoint().z + v_i[2] * v_width.z);
        // Derivatives
        Float dx, dy, dz;

        // Check if we are at boundaries
        if (v_i[0] == 0) {
            dx = (values[OffsetVoxel(v_i[0] + 1, v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                 inv_v_width.x;
        } else if (v_i[0] == dims[0] - 1) {
            dx = (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0] - 1, v_i[1], v_i[2])]) *
                 inv_v_width.x;
        } else {
            // Compute halfway derivatives and interpolate
            const Float dx_i_plus_12 =
                    (values[OffsetVoxel(v_i[0] + 1, v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                    inv_v_width.x;
            const Float dx_i_minus_12 =
                    (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0] - 1, v_i[1], v_i[2])]) *
                    inv_v_width.x;
            const Float tx = (p.x - voxel_min.x) * inv_v_width.x;
            dx = (1.f - tx) * dx_i_minus_12 + tx * dx_i_plus_12;
        }

        // Check if we are at boundaries
        if (v_i[1] == 0) {
            dy = (values[OffsetVoxel(v_i[0], v_i[1] + 1, v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                 inv_v_width.y;
        } else if (v_i[1] == dims[1] - 1) {
            dy = (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1] - 1, v_i[2])]) *
                 inv_v_width.y;
        } else {
            // Compute halfway derivatives and interpolate
            const Float dy_j_plus_12 =
                    (values[OffsetVoxel(v_i[0], v_i[1] + 1, v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                    inv_v_width.y;
            const Float dy_j_minus_12 =
                    (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1] - 1, v_i[2])]) *
                    inv_v_width.y;
            const Float ty = (p.y - voxel_min.y) * inv_v_width.y;
            dy = (1.f - ty) * dy_j_minus_12 + ty * dy_j_plus_12;
        }

        // Check if we are at boundaries
        if (v_i[2] == 0) {
            dz = (values[OffsetVoxel(v_i[0], v_i[1], v_i[2] + 1)] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                 inv_v_width.z;
        } else if (v_i[2] == dims[2] - 1) {
            dz = (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2] - 1)]) *
                 inv_v_width.z;
        } else {
            // Compute halfway derivatives and interpolate
            const Float dz_k_plus_12 =
                    (values[OffsetVoxel(v_i[0], v_i[1], v_i[2] + 1)] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2])]) *
                    inv_v_width.z;
            const Float dz_k_minus_12 =
                    (values[OffsetVoxel(v_i[0], v_i[1], v_i[2])] - values[OffsetVoxel(v_i[0], v_i[1], v_i[2] - 1)]) *
                    inv_v_width.z;
            const Float tz = (p.z - voxel_min.z) * inv_v_width.z;
            dz = (1.f - tz) * dz_k_minus_12 + tz * dz_k_plus_12;
        }

        return Vector3F(dx, dy, dz);
    }

    Vector3F MACGrid::NormalAt(int i, int j, int k) const {
        // TODO
    }

    void MACGrid::Refine(const int *new_dims) {
        // TODO
    }

    bool MACGrid::Intersect(Ray const &ray, Interaction *interaction) const {
        // The intersection procedure uses ray marching to check if we have an interaction with the stored surface

        // Current depth
        Float depth(0.f);

        // Find tollerance
        const float min_dist = std::min(MIN_DIST, std::min(std::min(v_width.x, v_width.y), v_width.z));

        for (int steps = 0; steps < MAX_STEPS; steps++) {
            // Compute distance from surface
            const Float distance = SDFAt(ray(depth));
            // Check if we are close enough to the surface
            if (distance < min_dist) {
                // Fill interaction
                interaction->p = ray(depth);

                // Estimate normal
                interaction->n = NormalAt(interaction->p);

                // Interaction parameter
                interaction->t = depth;

                // Outgoing direction
                interaction->wo = -Normalize(ray.d);
                // Set albedo to 1
                interaction->albedo = Spectrum(1.f);

                return true;
            }
            // Increase distance
            depth += distance;
            // Check for end
            if (depth > MAX_DIST) { return false; }
        }
        return false;
    }

    bool MACGrid::IntersectP(Ray const &ray) const {
        // The intersection procedure uses ray marching to check if we have a hit with the surface

        // Current depth
        Float depth(0.f);

        for (int steps = 0; steps < MAX_STEPS; steps++) {
            // Compute distance from surface
            const Float distance = SDFAt(ray(depth));
            // Check if we are close enough to the surface
            if (distance < MIN_DIST) { return true; }
            // Increase distance
            depth += distance;
            // Check for end
            if (depth > MAX_DIST) { return false; }
        }
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