//
// Created by simon on 02.06.17.
//

#ifndef DRDEMO_TRIANGLE_MESH_HPP
#define DRDEMO_TRIANGLE_MESH_HPP

#include <cstdint>
#include <memory>
#include "bvh.hpp"
#include "shape.hpp"

namespace drdemo {

    /**
     * This file defines the classes to store a triangle mesh loaded from an .obj file
     */

    // Forward declare TriangleMesh class
    class TriangleMesh;

    /**
     * Triangle storage information for the mesh
     * This struct hold the indices of each triangle (vertices and normals indices)
     */
    struct TriangleIndices {
        // Vertices indices
        uint32_t v[3];
        // Normal indices
        uint32_t n[3];

        // Default constructor
        TriangleIndices() = default;

        TriangleIndices(uint32_t v0, uint32_t v1, uint32_t v2,
                        uint32_t n0 = 0, uint32_t n1 = 0, uint32_t n2 = 0);
    };

    /**
     * Triangle class
     * 10.9.2017: Refactored to not use differentiable variables
     */
    class Triangle : public Shape {
    private:
        // Reference to the TriangleMesh holding the data
        TriangleMesh &mesh;
        // Index of the triangle data
        uint32_t const triangle_index;

    public:
        Triangle(TriangleMesh &mesh, uint32_t t_i);

        // Shape methods
        bool Intersect(Ray const &ray, Interaction *interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        BBOX BBox() const override;

        Vector3f Centroid() const override;

        std::string ToString() const override;

//        // Differentiable object methods
//        void GetDiffVariables(std::vector<Float const *> &vars) const override;
//
//        size_t GetNumVars() const noexcept override;
//
//        void UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) override;
//
//        void SetDiffVariables(std::vector<float> const &vals, size_t starting_index) override;
    };

    /**
     * TriangleMesh class, hold the information data for a triangle mesh
     *
     * This is now used to load a mesh and render it as target
     */
    class TriangleMesh : public Shape {
    private:
        // Friend class Triangle for practice
        friend class Triangle;

        // Loaded triangles
        std::vector<TriangleIndices> triangles;

        // Loaded vertices
        // std::vector<Vector3F> vertices;
        std::vector<Vector3f> vertices;

        // Loaded normals
        // std::vector<Vector3F> normals;
        std::vector<Vector3f> normals;

        // List of triangles
        std::vector<std::shared_ptr<const Shape> > shapes;

        // Intersection acceleration structure
        BVH bvh;

        // Adds the triangles in the mesh to a vector of shape
        void CreateTriangles();

    public:
        explicit TriangleMesh(std::string const &file_name);

        // Access mesh information
        inline size_t NumTriangles() const { return triangles.size(); }

        inline size_t NumVertices() const { return vertices.size(); }

        // Shape methods
        bool Intersect(Ray const &ray, Interaction *interaction) const override;

        bool IntersectP(Ray const &ray) const override;

        BBOX BBox() const override;

        Vector3f Centroid() const override;

        std::string ToString() const override;

//        // Differentiable object methods
//        void GetDiffVariables(std::vector<Float const *> &vars) const override;
//
//        size_t GetNumVars() const noexcept override;
//
//        void UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) override;
//
//        void SetDiffVariables(const std::vector<float> &vals, size_t starting_index) override;
    };

} // drdemo namespace

#endif //DRDEMO_TRIANGLE_MESH_HPP
