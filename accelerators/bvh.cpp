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
        // The range of shapes in the object list covered by this entry
        uint32_t start, end;
    };

    void BVH::Build() {
        // Create entry stack
        BVHBuildEntry todo[128];
        uint32_t stack_ptr = 0;

        // Constants for construction
        const uint32_t untouched = 0xffffffff;
        const uint32_t touched_twice = 0xfffffffd;

        // Push the root
        todo[stack_ptr++] = BVHBuildEntry(0, (uint32_t) shapes.size(), 0xfffffffc);

        // Declare local variables for construction
        BVHFlatNode node;
        // Reserve space for build nodes (num_objects * 2)
        flat_tree.reserve(shapes.size() * 2);

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
            BBOX bb = shapes[start]->BBox();
            BBOX bc = BBOX(shapes[start]->Centroid());
            for (uint32_t o = start + 1; o < end; o++) {
                bb.ExpandTo(shapes[o]->BBox());
                bc.ExpandTo(shapes[o]->Centroid());
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

            // Partition the list of shapes
            uint32_t mid = start;
            for (uint32_t i = start; i < end; i++) {
                // Check if the object centroid is below the split coordinate
                if (shapes[i]->Centroid()[split_axis] < split_coord) {
                    // Swap primitives position
                    std::swap(shapes[i], shapes[mid]);
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

    BVH::BVH(std::vector<std::shared_ptr<Shape> > &s, uint32_t leaf_size)
            : num_nodes(0), num_leafs(0), leaf_size(leaf_size), shapes(s), flat_tree() {
        if (!shapes.empty()) {
            // Build tree
            Build();
        }

        // Print tree data
        std::cout << ToString() << std::endl;
    }

    void BVH::Rebuild() {
        std::cout << "Rebuilding BVH..." << std::endl;
        // Clear flat tree data and rebuild
        flat_tree.clear();
        Build();
    }

    // Define tree traversal struct
    struct BVHTraversal {
        BVHTraversal() = default;

        BVHTraversal(uint32_t i, float min_t)
                : i(i), min_t(min_t) {}

        // Node
        uint32_t i;
        // Minimum hit for this node
        float min_t;
    };

    bool BVH::Intersect(Ray const &ray, Interaction *const interaction) const {
        if (shapes.empty()) { return false; }
        // Local used interval
        Float child_0_min, child_0_max;
        Float child_1_min, child_1_max;

        // Working set stack
        BVHTraversal todo[64];
        int32_t stack_ptr = 0;

        // Flag fot hit
        bool hit = false;

        // Push the root node on the stack
        todo[stack_ptr] = BVHTraversal(0, ray.t_min.GetValue());

        while (stack_ptr >= 0) {
            // Pop node to work on
            uint32_t ni = todo[stack_ptr].i;
            float near = todo[stack_ptr].min_t;
            stack_ptr--;
            const BVHFlatNode &node = flat_tree[ni];

            if (near > ray.t_max) { continue; }

            // Check if node is a leag
            if (node.right_offset == 0) {
                for (uint32_t o = 0; o < node.num_prims; o++) {
                    if (shapes[node.start + o]->Intersect(ray, interaction)) {
                        hit = true;
                    }
                }
            } else {
                // Not a leaf, check intersection with child BBox
                bool hit_c0 = flat_tree[ni + 1].bbox.Intersect(ray, &child_0_min, &child_0_max);
                bool hit_c1 = flat_tree[ni + node.right_offset].bbox.Intersect(ray, &child_1_min, &child_1_max);

                // Check if we hit both
                if (hit_c0 && hit_c1) {
                    // Check which child was closer
                    if (child_0_max <= child_1_min) {
                        // Add farther child first
                        todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                        todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                    } else {
                        todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                        todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                    }
                } else if (hit_c0) {
                    todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                } else if (hit_c1) {
                    todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                }
            }
        }

        return hit;
    }

    bool BVH::IntersectP(Ray const &ray) const {
        if (shapes.empty()) { return false; }
        // Local used interval
        Float child_0_min, child_0_max;
        Float child_1_min, child_1_max;

        // Working set stack
        BVHTraversal todo[64];
        int32_t stack_ptr = 0;

        // Push the root node on the stack
        todo[stack_ptr] = BVHTraversal(0, ray.t_min.GetValue());

        while (stack_ptr >= 0) {
            // Pop node to work on
            uint32_t ni = todo[stack_ptr].i;
            float near = todo[stack_ptr].min_t;
            stack_ptr--;
            const BVHFlatNode &node = flat_tree[ni];

            if (near > ray.t_max) { continue; }

            // Check if node is a leag
            if (node.right_offset == 0) {
                for (uint32_t o = 0; o < node.num_prims; o++) {
                    if (shapes[node.start + o]->IntersectP(ray)) {
                        return true;
                    }
                }
            } else {
                // Not a leaf, check intersection with child BBox
                bool hit_c0 = flat_tree[ni + 1].bbox.Intersect(ray, &child_0_min, &child_0_max);
                bool hit_c1 = flat_tree[ni + node.right_offset].bbox.Intersect(ray, &child_1_min, &child_1_max);

                // Check if we hit both
                if (hit_c0 && hit_c1) {
                    // Check which child was closer
                    if (child_0_max <= child_1_min) {
                        // Add farther child first
                        todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                        todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                    } else {
                        todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                        todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                    }
                } else if (hit_c0) {
                    todo[++stack_ptr] = BVHTraversal(ni + 1, child_0_min.GetValue());
                } else if (hit_c1) {
                    todo[++stack_ptr] = BVHTraversal(ni + node.right_offset, child_1_min.GetValue());
                }
            }
        }

        return false;
    }

    BBOX BVH::BBox() const {
        return flat_tree[0].bbox;
    }

    Vector3F BVH::Centroid() const {
        Vector3F extent = flat_tree[0].bbox.Extent();

        return flat_tree[0].bbox.MinPoint() + 0.5f * extent;
    }

    std::string BVH::ToString() const {
        return "BVH acceleration structure with " + std::to_string(num_nodes) + " nodes, "
               + std::to_string(num_leafs) + " and " + std::to_string(shapes.size()) + " shapes.";
    }

    void BVH::GetDiffVariables(std::vector<Float const *> &vars) const { // FIXME
        // Loop over the list of all Shapes and request variables
        for (auto const &shape : shapes) {
            shape->GetDiffVariables(vars);
        }
    }

    size_t BVH::GetNumVars() const noexcept { // FIXME
        size_t total_vars = 0;
        for (auto const &shape : shapes) {
            total_vars += shape->GetNumVars();
        }

        return total_vars;
    }

    void BVH::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) { // FIXME
        size_t used_vars = 0;
        for (auto &shape : shapes) {
            shape->UpdateDiffVariables(delta, starting_index + used_vars);
            used_vars += shape->GetNumVars();
        }
    }

} // drdemo namespace