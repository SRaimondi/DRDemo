//
// Created by simon on 02.06.17.
//

#include "triangle_mesh.hpp"

namespace drdemo {

    TriangleIndices::TriangleIndices(uint32_t v0, uint32_t v1, uint32_t v2,
                                     uint32_t n0, uint32_t n1, uint32_t n2) {
        // Set vertices indices
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        // Set normal indices
        n[0] = n0;
        n[1] = n1;
        n[2] = n2;
    }

    Triangle::Triangle(TriangleMesh const &mesh, uint32_t t_i)
            : mesh(mesh), triangle_index(t_i) {}

    bool Triangle::Intersect(Ray const &ray, Interaction *const interaction) const {
        return false;
    }

    bool Triangle::IntersectP(Ray const &ray) const {
        return false;
    }

    BBOX Triangle::BBox() const {
        return BBOX();
    }

    Vector3F Triangle::Centroid() const {
        return drdemo::Vector3F();
    }

    std::string Triangle::ToString() const {
        return std::__cxx11::string();
    }

    void Triangle::GetDiffVariables(std::vector<Float const *> &vars) const {

    }

    size_t Triangle::GetNumVars() const {
        return 0;
    }

    void Triangle::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) {

    }

} // drdemo namespace
