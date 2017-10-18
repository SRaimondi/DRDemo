//
// Created by Simon on 19.06.2017.
//

#ifndef DRDEMO_GRID_HPP
#define DRDEMO_GRID_HPP

#include "shape.hpp"
#include <memory>

namespace drdemo {

    /**
     * This file defines an implementation of a signed distance filed mesh representation using a 3D grid
     * We store the value of the signed distance function at each vertex
     *
     * The storage order is along x, y, z.
     */
    class SignedDistanceGrid : public Shape, public DiffObjectInterface {
    private:
        // Number of points for each axis (number of voxels is the value for a given axis minus 1)
        int num_points[3];
        // Total points
        int total_points;
        // Bounds of the Grid
        BBOX bounds;
        // Size of the voxels and inverse of width
        Vector3f width, inv_width;
        // Signed distance values, stored at each voxel intersection in x, y, z order
        // std::unique_ptr<Float[]> data;
        Float *data;

        // Private utility methods
        inline int OffsetPoint(int x, int y, int z) const {
            // Check for boundaries condition
            if (x >= num_points[0]) { x -= num_points[0]; }
            else if (x < 0) { x += num_points[0]; }

            if (y >= num_points[1]) { y -= num_points[1]; }
            else if (y < 0) { y += num_points[1]; }

            if (z >= num_points[2]) { z -= num_points[2]; }
            else if (z < 0) { z += num_points[2]; }

            return z * num_points[0] * num_points[1] + y * num_points[0] + x;
        }

        inline int OffsetVoxel(int x, int y, int z) const {
            return z * (num_points[0] - 1) * (num_points[1] - 1) + y * (num_points[0] - 1) + x;
        }

        // Convert 3d point to VOXEL coordinate given an axis (0: x, 1: y, 2: z)
        inline int PosToVoxel(const Vector3f &p, int axis) const {
            auto v_i = static_cast<int>((p[axis] - bounds.MinPoint()[axis]) * inv_width[axis]);

            // return v_i;
            return Clamp(v_i, 0, num_points[axis] - 2);
        }

        // Get the indices of the eight vertices of the voxel given his coordinates
        // The indices start from the corner with the smallest (x,y,z) coordinates and then follow the order of the grid
        void PointsIndicesFromVoxel(int x, int y, int z, int *indices) const;

        // Compute normal at given point inside the SDF using tri-linear interpolation
        Vector3F NormalAt(const Vector3F &p) const;

        // Estimate SDF normal at given point
//        Vector3F EstimateNormal(const Vector3F &p, float eps = EPS) const;

    public:
        // Default constructor, initialises an empty grid
        SignedDistanceGrid(int n_x, int n_y, int n_z, BBOX const &b);

        // Construct grid from values
        SignedDistanceGrid(int n_x, int n_y, int n_z, BBOX const &b, float const *raw_data);

        // Construct grid from file, input file from https://github.com/christopherbatty/SDFGen
        explicit SignedDistanceGrid(const std::string &sdf_file);

        // Destructor
        ~SignedDistanceGrid();

        // Access Grid point at given indices
        inline Float const &operator()(int x, int y, int z) const {
            return data[OffsetPoint(x, y, z)];
        }

        inline Float &operator()(int x, int y, int z) {
            return data[OffsetPoint(x, y, z)];
        }

        // Compute the value of the Signed Distance Function sampled by the grid using trilinear interpolation given
        // a point in the grid
        Float ValueAt(const Vector3F &p) const;

        // Compute the normal at a given grid point, last argument is to use backward difference
        Vector3F NormalAtPoint(int x, int y, int z /* , bool bd = false */) const;

        // Access size of the voxels
        inline Vector3f const &VoxelSize() const { return width; }

        inline Vector3f const &InvVoxelSize() const { return inv_width; }

        // Compute coordinates of points at given indices
        inline Vector3f CoordsAt(int x, int y, int z) const {
            return Vector3f(bounds.MinPoint().x + x * width.x,
                            bounds.MinPoint().y + y * width.y,
                            bounds.MinPoint().z + z * width.z);
        }

        // Get grid dimensions
        inline int Size(int axis) const {
            return num_points[axis];
        }

        // Convert 3 indices to linear
        inline int LinearIndex(int x, int y, int z) const {
            return OffsetPoint(x, y, z);
        }

        // Convert linear index to 3 indices
        // TODO Check this again, should be correct
        inline void IndicesFromLinear(int linear_index, int &x, int &y, int &z) const {
            // Compute z index
            z = linear_index / (num_points[0] * num_points[1]);
            linear_index -= z * (num_points[0] * num_points[1]);
            // Compute y index
            y = linear_index / num_points[0];
            linear_index -= y * num_points[0];
            // Last index
            x = linear_index;
        }

        // Refine grid to new higher resolution resolution
        void Refine(int const *new_dims);

        // Shape methods
        bool Intersect(Ray const &ray, Interaction *interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        BBOX BBox() const override;

        Vector3f Centroid() const override;

        std::string ToString() const override;

        // Differentiable object methods
        void GetDiffVariables(std::vector<Float const *> &vars) const override;

        size_t GetNumVars() const noexcept override;

        void UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) override;

        void SetDiffVariables(const std::vector<float> &vals, size_t starting_index) override;
    };


//    /**
//     * Reinitialization equation utilities
//     */
//
//    // Compute squared norm of the gradient using forward finite difference method
//    float GradNorm2(const SignedDistanceGrid &grid, int x, int y, int z);
//
//    // Compute sign grid
//    void ComputeSignGrid(const SignedDistanceGrid &grid, float *sign_grid);
//
//    // Compute right part of PDF
//    float RightTermReinitializePDE(const SignedDistanceGrid &grid, const float *sign_grid, int x, int y, int z);
//
//    // Reinitialize SDF
//    void ReinitializeSDF(SignedDistanceGrid &grid, float dt, size_t max_iters, float tolerance, float band);

} // drdemo namespace

#endif //DRDEMO_GRID_HPP
