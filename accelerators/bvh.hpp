//
// Created by simon on 02.06.17.
//

#ifndef DRDEMO_BVH_HPP
#define DRDEMO_BVH_HPP

#include <memory>
#include "shape.hpp"

namespace drdemo {

    /**
     * Define BVH flat node structure
     */
    struct BVHFlatNode {
        // Node BBOX
        BBOX bbox;
        // Starting index of Shape
        uint32_t start;
        // Number of primitives in the node
        uint32_t num_prims;
        // Right node offset
        uint32_t right_offset;
    };

    /**
     * Define BVH acceleration structure for fast triangle mesh ray intersection
     * TODO: Can we just use float for BVH?
     */
    class BVH /* : public Shape */ {
    private:
        // Number of nodes, leafs and leaf size
        uint32_t num_nodes;
        uint32_t num_leafs;
        uint32_t leaf_size;
        // Shapes in the BVH
        std::vector<std::shared_ptr<Shape> > &shapes;
        // Flat tree data
        std::vector<BVHFlatNode> flat_tree;

    public:
        BVH(std::vector<std::shared_ptr<Shape> > &s, uint32_t leaf_size = 4);

        // Build the tree
        void Build();

        // Rebuild BVH data structure, needed when we update the triangles positions
        void Rebuild();

        // Shape methods
        // bool Intersect(Ray const &ray, Interaction *const interaction) const override;

        bool Intersect(Ray const &ray, Interaction *const interaction) const;

        // bool IntersectP(Ray const &ray) const override;

        bool IntersectP(Ray const &ray) const;

        // BBOX BBox() const override;

        // Vector3F Centroid() const override;

        // std::string ToString() const override;

        std::string ToString() const;

        // Differentiable object methods
        // void GetDiffVariables(std::vector<Float const *> &vars) const override;

        // size_t GetNumVars() const noexcept override;

        // void UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index = 0) override;
    };

} // drdemo namespace

#endif //DRDEMO_BVH_HPP
