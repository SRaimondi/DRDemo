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
        // Check if we are outside the BBOX
        if (!bounds.Inside(p_f)) {
            return Float(bounds.Distance(p_f) + 0.001f);
            // TODO The 0.001 is there to make the next point go inside the Grid if the ray direction is perpendicular to
            // the normal of the grid intersected face
        }

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

    Vector3F SignedDistanceGrid::EstimateNormal(const Vector3F &p, float eps) const {
        return Normalize(
                Vector3F(ValueAt(p + Vector3F(eps, 0.f, 0.f)) - ValueAt(p - Vector3F(eps, 0.f, 0.f)),
                         ValueAt(p + Vector3F(0.f, eps, 0.f)) - ValueAt(p - Vector3F(0.f, eps, 0.f)),
                         ValueAt(p + Vector3F(0.f, 0.f, eps)) - ValueAt(p - Vector3F(0.f, 0.f, eps)))
        );
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
            width[axis] = extent[axis] / static_cast<float>(num_points[axis] - 1);
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
            width[axis] = extent[axis] / static_cast<float>(num_points[axis] - 1);
            inv_width[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
        }
        // Copy values
        for (size_t i = 0; i < total_points; ++i) {
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
                interaction->n = EstimateNormal(interaction->p);
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
        // The intersection procedure uses ray marching to check if we have an interaction with the stored surface

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
        for (size_t i = 0; i < total_points; i++) {
            content += std::to_string(data[i].GetValue());
            if (i != total_points - 1) {
                content += ", ";
            }
        }
        content += ")";

        return content;
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
        // THIS DOES NOT WORK AT ALL
        size_t used_vars = 0;
        for (size_t i = 0; i < total_points; ++i) {
            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
        }
    }

//    void SignedDistanceGrid::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {
//        // FIXME Initial simple attempt to propagate gradient directly into the data
//        size_t used_vars = 0;
//        for (size_t i = 0; i < total_points; ++i) {
//            data[i].SetValue(data[i].GetValue() + delta[starting_index + used_vars++]);
//        }
//    }

} // drdemo namespace
