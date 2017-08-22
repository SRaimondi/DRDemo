//
// Created by Simon on 19.06.2017.
//

#include "grid.hpp"

namespace drdemo {

    void SignedDistanceGrid::PointsIndicesFromVoxel(int x, int y, int z, int *const indices) const {
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
        // Check if we are outside the BBOX
        if (!bounds.Inside(p_f)) {
            return Float(bounds.Distance(p_f) + 0.001f);
            // FIXME The 0.001 is there to make the next point go inside the Grid if the ray direction is perpendicular to
            // the normal of the grid intersected face
        }

        // Get voxel indices
        int voxel_i[3];
        for (int i = 0; i < 3; i++) {
            voxel_i[i] = PosToVoxel(p_f, i);
        }
        // Get point indices for the given voxel
        int indices[8];
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

    Vector3F SignedDistanceGrid::NormalAt(const Vector3F &p) const {
        // Convert position
        const Vector3f p_f = Tofloat(p);
        // Get voxel indices
        int voxel_i[3];
        for (int i = 0; i < 3; i++) {
            voxel_i[i] = PosToVoxel(p_f, i);
        }
        // Get point indices for the given voxel
        int indices[8];
        PointsIndicesFromVoxel(voxel_i[0], voxel_i[1], voxel_i[2], indices);

        // Compute minimum point of voxel
        const Vector3f voxel_min(bounds.MinPoint().x + voxel_i[0] * width.x,
                                 bounds.MinPoint().y + voxel_i[1] * width.y,
                                 bounds.MinPoint().z + voxel_i[2] * width.z);

        int x, y, z;
        // Compute normal at the 8 vertices of the voxel
        IndicesFromLinear(indices[0], x, y, z);
        const Vector3F n0 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[1], x, y, z);
        const Vector3F n1 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[2], x, y, z);
        const Vector3F n2 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[3], x, y, z);
        const Vector3F n3 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[4], x, y, z);
        const Vector3F n4 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[5], x, y, z);
        const Vector3F n5 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[6], x, y, z);
        const Vector3F n6 = NormalAtPoint(x, y, z);
        IndicesFromLinear(indices[7], x, y, z);
        const Vector3F n7 = NormalAtPoint(x, y, z);

        // Interpolate normals
        const Float tx = (p.x - voxel_min.x) * inv_width.x;
        const Vector3F n01 = (1.f - tx) * n0 + tx * n1;
        const Vector3F n23 = (1.f - tx) * n2 + tx * n3;
        const Vector3F n45 = (1.f - tx) * n4 + tx * n5;
        const Vector3F n67 = (1.f - tx) * n6 + tx * n7;

        // Interpolate along y
        const Float ty = (p.y - voxel_min.y) * inv_width.y;
        const Vector3F nf1 = (1.f - ty) * n01 + ty * n23;
        const Vector3F nf2 = (1.f - ty) * n45 + ty * n67;

        // Return final value interpolated along z
        const Float tz = (p.z - voxel_min.z) * inv_width.z;

        return (1.f - tz) * nf1 + tz * nf2;
    }

    Vector3F SignedDistanceGrid::NormalAtPoint(int x, int y, int z) const {
        // TODO Maybe we could use some higher order scheme?
        // FIXME This needs to be VERY CAREFULLY designed now
        // Compute normal at given grid point using forward difference method

        return Vector3F(
                (data[LinearIndex(x + 1, y, z)] - data[LinearIndex(x, y, z)]) * inv_width.x,
                (data[LinearIndex(x, y + 1, z)] - data[LinearIndex(x, y, z)]) * inv_width.y,
                (data[LinearIndex(x, y, z + 1)] - data[LinearIndex(x, y, z)]) * inv_width.z
        );
    }

    Vector3F SignedDistanceGrid::EstimateNormal(const Vector3F &p, float eps) const {
        return Normalize(
                Vector3F(ValueAt(p + Vector3F(eps, 0.f, 0.f)) - ValueAt(p - Vector3F(eps, 0.f, 0.f)),
                         ValueAt(p + Vector3F(0.f, eps, 0.f)) - ValueAt(p - Vector3F(0.f, eps, 0.f)),
                         ValueAt(p + Vector3F(0.f, 0.f, eps)) - ValueAt(p - Vector3F(0.f, 0.f, eps)))
        );
    }

    SignedDistanceGrid::SignedDistanceGrid(int n_x, int n_y, int n_z, BBOX const &b)
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
            width[axis] = extent[axis] / static_cast<float>(num_points[axis] - 1);
            inv_width[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
        }
    }

    SignedDistanceGrid::SignedDistanceGrid(int n_x, int n_y, int n_z, BBOX const &b,
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
        const Vector3f extent = bounds.Extent();
        for (int axis = 0; axis < 3; ++axis) {
            width[axis] = extent[axis] / static_cast<float>(num_points[axis] - 1);
            inv_width[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
        }
        // Copy values
        for (int i = 0; i < total_points; ++i) {
            data[i] = raw_data[i];
        }
    }

    bool SignedDistanceGrid::Intersect(Ray const &ray, Interaction *const interaction) const {
        // The intersection procedure uses ray marching to check if we have an interaction with the stored surface

        // Current depth
        Float depth(0.f);

        for (int steps = 0; steps < MAX_STEPS; steps++) {
            // Compute distance from surface
            const Float distance = ValueAt(ray(depth));
            // Check if we are close enough to the surface
            if (distance < MIN_DIST) {
                // Fill interaction
                interaction->p = ray(depth);

                // Estimate normal with finite difference
                // interaction->n = EstimateNormal(interaction->p);
                interaction->n = Normalize(NormalAt(interaction->p));
                // interaction->n = NormalAt(ray(depth));
//                if (Length(interaction->n).GetValue() - 1.f > 0.001f) {
//                    std::cout << "Normal length: " << Length(interaction->n).GetValue() << std::endl;
//                }

                // Interaction parameter
                interaction->t = depth;

                // Outgoing direction
                interaction->wo = -Normalize(ray.d);

                return true;
            }
            // Increase distance
            depth += distance;
            // Check for end
            if (depth > MAX_DIST) { return false; }
        }
        return false;
    }

    bool SignedDistanceGrid::IntersectP(Ray const &ray) const {
        // The intersection procedure uses ray marching to check if we have a hit with the surface

        // Current depth
        Float depth(0.f);

        for (int steps = 0; steps < MAX_STEPS; steps++) {
            // Compute distance from surface
            const Float distance = ValueAt(ray(depth));
            // Check if we are close enough to the surface
            if (distance < MIN_DIST) { return true; }
            // Increase distance
            depth += distance;
            // Check for end
            if (depth > MAX_DIST) { return false; }
        }
        return false;
    }

    BBOX SignedDistanceGrid::BBox() const {
        return bounds;
    }

    Vector3f SignedDistanceGrid::Centroid() const {
        return Vector3f(0.f, 0.f, 0.f);
    }

    std::string SignedDistanceGrid::ToString() const {
        std::string content("(");
        for (int i = 0; i < total_points; i++) {
            content += std::to_string(data[i].GetValue());
            if (i != total_points - 1) {
                content += ", ";
            }
        }
        content += ")";

        return content;
    }

    void SignedDistanceGrid::GetDiffVariables(std::vector<Float const *> &vars) const {
        for (int i = 0; i < total_points; ++i) {
            vars.push_back(&data[i]);
        }
    }

    size_t SignedDistanceGrid::GetNumVars() const noexcept {
        return static_cast<size_t>(total_points);
    }

//    void SignedDistanceGrid::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {
//        // FIXME Initial simple attempt to propagate gradient directly into the data
//        // THIS DOES NOT WORK AT ALL
//        size_t used_vars = 0;
//        for (int i = 0; i < total_points; ++i) {
//            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
//        }
//    }

    void SignedDistanceGrid::UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) {
        // Second attempt, trying to correct the SDF after propagating directly the gradient into it


        // Propagate gradient inside the grid
        size_t used_vars = 0;
        for (int i = 0; i < total_points; ++i) {
            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
        }
        // Reinitialise the grid to be a SDF after gradient update
        // ReinitializeSDF(*this, 0.0001f, 10000, 0.005f, 100.f);
    }

    void SignedDistanceGrid::SetDiffVariables(const std::vector<float> &vals, size_t starting_index) {
        size_t used_vars = 0;
        for (int i = 0; i < total_points; ++i) {
            data[i].SetValue(vals[starting_index + used_vars++]);
        }
    }

    float GradNorm2(const SignedDistanceGrid &grid, int x, int y, int z) {
        // Compute derivatives using finite difference
        const float grad_x = (grid(x + 1, y, z).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().x;
        const float grad_y = (grid(x, y + 1, z).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().y;
        const float grad_z = (grid(x, y, z + 1).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().z;

        return grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
    }

    void ComputeSignGrid(const SignedDistanceGrid &grid, float *sign_grid) {
        for (int z = 0; z < grid.Size(2); z++) {
            for (int y = 0; y < grid.Size(1); y++) {
                for (int x = 0; x < grid.Size(0); x++) {
                    // Get SDF value
                    const float sdf_v = grid(x, y, z).GetValue();
                    // Compute gradient squared norm
                    const float grad_2 = GradNorm2(grid, x, y, z);
                    // Update sign
                    sign_grid[grid.LinearIndex(x, y, z)] = sdf_v /
                                                           std::sqrt(sdf_v * sdf_v + grad_2 * EPS * EPS);
                }
            }
        }
    }

    float RightTermReinitializePDE(const SignedDistanceGrid &grid, const float *sign_grid, int x, int y, int z) {
        // Compute gradient norm term using Godunov method
        const float dphi_dx_plus = (grid(x + 1, y, z).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().x;
        const float dphi_dx_minus = (grid(x, y, z).GetValue() - grid(x - 1, y, z).GetValue()) * grid.InvVoxelSize().x;

        const float dphi_dy_plus = (grid(x, y + 1, z).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().y;
        const float dphi_dy_minus = (grid(x, y, z).GetValue() - grid(x, y - 1, z).GetValue()) * grid.InvVoxelSize().y;

        const float dphi_dz_plus = (grid(x, y, z + 1).GetValue() - grid(x, y, z).GetValue()) * grid.InvVoxelSize().z;
        const float dphi_dz_minus = (grid(x, y, z).GetValue() - grid(x, y, z - 1).GetValue()) * grid.InvVoxelSize().z;

        // Check sign and compute gradient norm
        float grad_x_sq, grad_y_sq, grad_z_sq;

        if (sign_grid[grid.LinearIndex(x, y, z)] > 0.f) {
            grad_x_sq = std::max(PositiveSQ(dphi_dx_minus), NegativeSQ(dphi_dx_plus));
            grad_y_sq = std::max(PositiveSQ(dphi_dy_minus), NegativeSQ(dphi_dy_plus));
            grad_z_sq = std::max(PositiveSQ(dphi_dz_minus), NegativeSQ(dphi_dz_plus));
        } else {
            grad_x_sq = std::max(PositiveSQ(dphi_dx_plus), NegativeSQ(dphi_dx_minus));
            grad_y_sq = std::max(PositiveSQ(dphi_dy_plus), NegativeSQ(dphi_dy_minus));
            grad_z_sq = std::max(PositiveSQ(dphi_dz_plus), NegativeSQ(dphi_dz_minus));
        }

        const float grad_norm = std::sqrt(grad_x_sq + grad_y_sq + grad_z_sq);

        return sign_grid[grid.LinearIndex(x, y, z)] * (1.f - grad_norm);
    }

    void ReinitializeSDF(SignedDistanceGrid &grid, float dt, size_t max_iters, float tolerance, float band) {
        std::cout << "Reinitializing SDF..." << std::endl;
        // Allocate space for storing the sign grid
        auto sign_grid = new float[grid.GetNumVars()];

        // Compute initial sign grid
        ComputeSignGrid(grid, sign_grid);

        // Allocate space to store dphi_dt
        auto sdf_dt = new float[grid.GetNumVars()];
        float max_dt = tolerance + 1.f;
        size_t iters = 0;

        // Construction loop
        while (iters < max_iters && max_dt > tolerance) {
            max_dt = 0.f;

            // Loop over all grid and compute the right term of our PDE
            for (int z = 0; z < grid.Size(2); z++) {
                for (int y = 0; y < grid.Size(1); y++) {
                    for (int x = 0; x < grid.Size(0); x++) {
                        const int linear_index = grid.LinearIndex(x, y, z);
                        // Check if value is inside the band
                        if (std::abs(grid(x, y, z).GetValue()) >= band) {
                            sdf_dt[linear_index] = 0.f;
                        } else {
                            // Compute right part of PDE
                            sdf_dt[linear_index] = RightTermReinitializePDE(grid, sign_grid, x, y, z);
                            // Check for new max_dt
                            if (std::abs(sdf_dt[linear_index]) > max_dt) {
                                max_dt = std::abs(sdf_dt[linear_index]);
                            }
                        }
                    }
                }
            }

            // Update our SDF
            for (int z = 0; z < grid.Size(2); z++) {
                for (int y = 0; y < grid.Size(1); y++) {
                    for (int x = 0; x < grid.Size(0); x++) {
                        grid(x, y, z).SetValue(grid(x, y, z).GetValue() + dt * sdf_dt[grid.LinearIndex(x, y, z)]);
                    }
                }
            }

            // Update sign grid
            ComputeSignGrid(grid, sign_grid);
//            if (iters != 0 && iters % 100 == 0) {
//                std::cout << iters << " iterations..." << std::endl;
//            }
            iters++;
        }

        std::cout << "Reinitialised SDF in " << iters << " iterations." << std::endl;

        // Free sign grid and sdf_dt
        delete[] sign_grid;
        delete[] sdf_dt;
    }

} // drdemo namespace
