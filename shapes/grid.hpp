//
// Created by Simon on 19.06.2017.
//

#ifndef DRDEMO_GRID_HPP
#define DRDEMO_GRID_HPP

#include "shape.hpp"

namespace drdemo {

    /**
     * This file defines an implementation of a signed distance filed mesh representation using a 3D grid
     * We store the value of the signed distance function at each vertex
     *
     * The grid is centered at the origin for simplicity
     * The storage order is along x, y, z.
     */
    class SignedDistanceGrid : public Shape {
    private:
        // Number of points for each axis (number of voxels is the value for a given axis minus 1)
        size_t num_points[3];
        // Total points
        size_t total_points;
        // Bounds of the Grid
        BBOX bounds;
        // Size of the voxels and inverse of width
        Vector3f width, inv_width;
        // Signed distance values, stored at each voxel intersection
        Float *data;

        // Private utility methods
        inline size_t OffsetPoint(int x, int y, int z) const {
            return z * num_points[0] * num_points[1] + y * num_points[0] + x;
        }

        inline size_t OffsetVoxel(int x, int y, int z) const {
            return z * (num_points[0] - 1) * (num_points[1] - 1) + y * (num_points[0] - 1) + x;
        }

    public:
        // Default constructor, initialises an empty grid
        SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b);

        // Construct grid from values
        SignedDistanceGrid(size_t n_x, size_t n_y, size_t n_z, BBOX const &b, float const *const raw_data);

        // Destructor
        ~SignedDistanceGrid();

        // Access Grid point at given indices
        inline Float const &operator()(int x, int y, int z) const {
            return data[OffsetPoint(x, y, z)];
        }

        inline Float &operator()(int x, int y, int z) {
            return data[OffsetPoint(x, y, z)];
        }

        // Shape methods
        bool Intersect(Ray const &ray, Interaction *const interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        BBOX BBox() const override;

        Vector3f Centroid() const override;

        std::string ToString() const override;

        // Differentiable object methods
        void GetDiffVariables(std::vector<Float const *> &vars) const override;

        size_t GetNumVars() const noexcept override;

        void UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index = 0) override;
    };


} // drdemo namespace

#endif //DRDEMO_GRID_HPP
