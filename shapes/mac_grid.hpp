//
// Created by simon on 18.10.17.
//

#ifndef DRDEMO_MAC_GRID_HPP
#define DRDEMO_MAC_GRID_HPP

#include "shape.hpp"

namespace drdemo {

    /**
     * This file defines an implementation of a MAC grid to represent a SDF
     * We store the value at the center of the voxels so we can compute the normal
     * at the boundaries between the voxels and interpolate on that
     *
     * The grid boundaries are represented with a bounding box
     * The storage order is x, y, z
     */
    class MACGrid : public Shape, public DiffObjectInterface {
    private:
        // Number of voxels along each axis
        int dims[3];
        // Total number of points
        int total_voxels;
        // Grid bounds
        BBOX bounds;
        // Size and inverse size of the voxel
        Vector3f v_width, inv_v_width;
        // Signed distance function values stored at voxels centers
        Float *values;

        // Compute offset of voxel given the 3 indices
        inline int OffsetVoxel(int i, int j, int k) const {
            return k * (dims[1] * dims[0]) + j * dims[0] + i;
        }

        // Convert 3d point to voxel coordinate given an axis (x:0, y:1, z:2)
        inline int PosToVoxel(const Vector3f &p, int axis) const {
            auto index = static_cast<int>((p[axis] - bounds.MinPoint()[axis]) * inv_v_width[axis]);
            return Clamp(index, 0, dims[axis] - 1);
        }

    public:
        // Create and empty grid
        MACGrid(int nx, int ny, int nz, const BBOX &b);

        // Destructor
        ~MACGrid();

        // Access voxel value at given indices
        inline Float const &operator()(int i, int j, int k) const {
            return values[OffsetVoxel(i, j, k)];
        }

        inline Float &operator()(int i, int j, int k) {
            return values[OffsetVoxel(i, j, k)];
        }

        // Compute SDF value at given point
        Float SDFAt(const Vector3F &p) const;

        // Compute coordinates of points at given indices
        inline Vector3f PointAt(int i, int j, int k) const {
            return Vector3f(bounds.MinPoint().x + i * v_width.x,
                            bounds.MinPoint().y + j * v_width.y,
                            bounds.MinPoint().z + k * v_width.z);
        }

        // Get grid dimension along axis
        inline int Dim(int axis) const { return dims[axis]; }

        // Refine grid
        void Refine(const int new_dims[3]);

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

} // drdemo namespace

#endif //DRDEMO_MAC_GRID_HPP
