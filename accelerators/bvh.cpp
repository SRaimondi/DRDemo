//
// Created by simon on 02.06.17.
//

#include "bvh.hpp"

namespace drdemo {

    // Build node data structure
    struct BVHBuildEntry {
        BVHBuildEntry() = default;

        BVHBuildEntry(uint32_t s, uint32_t e, uint32_t p)
                : parent(p), start(s), end(e) {}

        // If non-zero then this is the index of the parent (used in offsets)
        uint32_t parent;
        // The range of objects in the object list covered by this entry
        uint32_t start, end;
    };

    void BVH::Build() {
        // Create entryy stack
        BVHBuildEntry todo[128];
        uint32_t stack_ptr = 0;

        // Constants for construction
        const uint32_t untouched = 0xffffffff;
        const uint32_t touched_twice = 0xfffffffd;

        // Push the root
        todo[stack_ptr++] = BVHBuildEntry(0, (uint32_t) objects.size(), 0xfffffffc);

        // Declare local variables for construction
        BVHFlatNode node;
        // Reserve space for build nodes (num_objects * 2)
        flat_tree.reserve(objects.size() * 2);

        // Loop until we have object on the stack to do
        while (stack_ptr > 0) {
            // Pop item from the stack
            BVHBuildEntry &build_node = todo[--stack_ptr];
            uint32_t start = build_node.start;
            uint32_t end = build_node.end;
            uint32_t num_prims = end - start;

            num_nodes++;
            // Set node values
            node.start = start;
            node.num_prims = num_prims;
            node.right_offset = untouched;

            // Compute BBox for the node
            BBOX bb = objects[start]->BBox();
            BBOX bc = BBOX(objects[start]->Centroid());
            for (uint32_t o = start + 1; o < end; o++) {
                bb.ExpandTo(objects[o]->BBox());
                bc.ExpandTo(objects[o]->Centroid());
            }
            // Set node BBox
            node.bbox = bb;

            // Check how many primitives we have in this node
            // Set as leaf if <= leaf_size
            if (num_prims <= leaf_size) {
                node.right_offset = 0;
                num_leafs++;
            }

            // Add node
            flat_tree.push_back(node);

            // Check if this is the root node
            if (build_node.parent != 0xfffffffc) {
                flat_tree[build_node.parent].right_offset--;
                // If this is the second touch, set offset
                if (flat_tree[build_node.parent].right_offset == touched_twice) {
                    flat_tree[build_node.parent].right_offset = num_nodes - 1 - build_node.parent;
                }
            }

            // Check if this is a leaf
            if (node.right_offset == 0) { continue; }

            // Find split dimensions
            uint32_t split_axis = bc.MaxDimension();
            // Split on the center of the longest axis
            float split_coord = 0.5f * (bc.MinPoint()[split_axis] + bc.MaxPoint()[split_axis]).GetValue();

            // Partition the list of objects
            uint32_t mid = start;
            for (uint32_t i = start; i < end; i++) {
                // Check if the object centroid is below the split coordinate
                if (objects[i]->Centroid()[split_axis] < split_coord) {
                    // Swap primitives position
                    std::swap(objects[i], objects[mid]);
                    ++mid;
                }
            }

            // If we get a bad split, change center
            if (mid == start || mid == end) {
                mid = start + (end - start) / 2;
            }

            // Push right and left child
            todo[stack_ptr++] = BVHBuildEntry(mid, end, num_nodes - 1);
            todo[stack_ptr++] = BVHBuildEntry(start, mid, num_nodes - 1);
        }
    }

    BVH::BVH(std::vector<std::shared_ptr<Shape const> > &objs, uint32_t leaf_size)
            : num_nodes(0), num_leafs(0), leaf_size(leaf_size), objects(objs), flat_tree() {
        if (!objects.empty()) {
            // Build tree
            Build();
        }

        // Print tree data
        std::cout << "Built BVH with " << num_nodes << " nodes, "
                  << num_leafs << " leafs and " << objects.size() << " objects" << std::endl;
    }

    void BVH::Rebuild() {

    }

    bool BVH::Intersect(Ray const &ray, Interaction *const interaction) const {
        return false;
    }

    bool BVH::IntersectP(Ray const &ray) const {
        return false;
    }

    BBOX BVH::BBox() const {
        return BBOX();
    }

    Vector3F BVH::Centroid() const {
        return drdemo::Vector3F();
    }

    std::string BVH::ToString() const {
        return std::string();
    }

    void BVH::GetDiffVariables(std::vector<Float const *> &vars) const {

    }

    size_t BVH::GetNumVars() const noexcept {
        return 0;
    }

    void BVH::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {

    }

} // drdemo namespace