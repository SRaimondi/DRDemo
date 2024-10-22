//
// Created by simon on 02.06.17.
//

#include "triangle_mesh.hpp"
#include "iofile.hpp"

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

    Triangle::Triangle(TriangleMesh &mesh, uint32_t t_i)
            : mesh(mesh), triangle_index(t_i) {}

    bool Triangle::Intersect(Ray const &ray, Interaction *const interaction) const {
        // Get vertices from mesh
        Vector3f const &v0 = mesh.vertices[mesh.triangles[triangle_index].v[0]];
        Vector3f const &v1 = mesh.vertices[mesh.triangles[triangle_index].v[1]];
        Vector3f const &v2 = mesh.vertices[mesh.triangles[triangle_index].v[2]];

        // Compute e1 and e2
        Vector3f const e1 = v1 - v0;
        Vector3f const e2 = v2 - v0;

        // Compute variables for triangle intersection
        Vector3f s1 = Cross(Tofloat(ray.d), e2);
        float divisor = Dot(s1, e1);
        if (divisor == 0.f) { return false; }
        float inv_divisor = 1.f / divisor;

        // Compute first barycentric coordinate
        Vector3f s = Tofloat(ray.o) - v0;
        float b1 = Dot(s, s1) * inv_divisor;
        if (b1 < 0.f || b1 > 1.f) { return false; }

        // Compute second barycentric coordinate
        Vector3f s2 = Cross(s, e1);
        float b2 = Dot(Tofloat(ray.d), s2) * inv_divisor;
        if (b2 < 0.f || b1 + b2 > 1.f) { return false; }

        // Compute t of the intersection
        float t = Dot(e2, s2) * inv_divisor;
        // Check we are inside ray boundaries
        if (t < ray.t_min || t > ray.t_max) { return false; }

        // Update ray max
        ray.t_max = t;

        // Fill interaction
        interaction->p = ray(t);

        if (mesh.normals.empty()) {
            interaction->n = ToFloat(Normalize(Cross(e1, e2)));
        } else {
            interaction->n = ToFloat(Normalize((1.f - b1 - b2) * mesh.normals[mesh.triangles[triangle_index].n[0]] +
                                               b1 * mesh.normals[mesh.triangles[triangle_index].n[1]] +
                                               b2 * mesh.normals[mesh.triangles[triangle_index].n[2]]));
        }
        interaction->t = t;
        interaction->wo = Normalize(-ray.d);
        // FIXME For the moment we fix the albedo of triangle mesh to be 1
        interaction->albedo = Spectrum(1.f);

        return true;
    }

    bool Triangle::IntersectP(Ray const &ray) const {
        // Get vertices from mesh
        Vector3f const &v0 = mesh.vertices[mesh.triangles[triangle_index].v[0]];
        Vector3f const &v1 = mesh.vertices[mesh.triangles[triangle_index].v[1]];
        Vector3f const &v2 = mesh.vertices[mesh.triangles[triangle_index].v[2]];

        // Compute e1 and e2
        Vector3f const e1 = v1 - v0;
        Vector3f const e2 = v2 - v0;

        // Compute variables for triangle intersection
        Vector3f s1 = Cross(Tofloat(ray.d), e2);
        float divisor = Dot(s1, e1);
        if (divisor == 0.f) { return false; }
        float inv_divisor = 1.f / divisor;

        // Compute first barycentric coordinate
        Vector3f s = Tofloat(ray.o) - v0;
        float b1 = Dot(s, s1) * inv_divisor;
        if (b1 < 0.f || b1 > 1.f) { return false; }

        // Compute second barycentric coordinate
        Vector3f s2 = Cross(s, e1);
        float b2 = Dot(Tofloat(ray.d), s2) * inv_divisor;
        if (b2 < 0.f || b1 + b2 > 1.f) { return false; }

        // Compute t of the intersection
        float t = Dot(e2, s2) * inv_divisor;
        // Check we are inside ray boundaries
        return (t > ray.t_min && t < ray.t_max);
    }

    BBOX Triangle::BBox() const {
        // Get vertices from mesh
        Vector3f const &v0 = mesh.vertices[mesh.triangles[triangle_index].v[0]];
        Vector3f const &v1 = mesh.vertices[mesh.triangles[triangle_index].v[1]];
        Vector3f const &v2 = mesh.vertices[mesh.triangles[triangle_index].v[2]];

        BBOX bbox(v0, v1);
        bbox.ExpandTo(v2);

        return bbox;
    }

    Vector3f Triangle::Centroid() const {
        // Get vertices from mesh
        Vector3f const &v0 = mesh.vertices[mesh.triangles[triangle_index].v[0]];
        Vector3f const &v1 = mesh.vertices[mesh.triangles[triangle_index].v[1]];
        Vector3f const &v2 = mesh.vertices[mesh.triangles[triangle_index].v[2]];

        return (1.f / 3.f) * (v0 + v1 + v2);
    }

    std::string Triangle::ToString() const {
        // Get vertices from mesh
        Vector3f const &v0 = mesh.vertices[mesh.triangles[triangle_index].v[0]];
        Vector3f const &v1 = mesh.vertices[mesh.triangles[triangle_index].v[1]];
        Vector3f const &v2 = mesh.vertices[mesh.triangles[triangle_index].v[2]];

        return "Vertices: \n(" +
               std::to_string(v0.x) + ", " +
               std::to_string(v0.y) + ", " +
               std::to_string(v0.z) + ")\n" +
               "(" +
               std::to_string(v1.x) + ", " +
               std::to_string(v1.y) + ", " +
               std::to_string(v1.z) + ")\n" +
               "(" +
               std::to_string(v2.x) + ", " +
               std::to_string(v2.y) + ", " +
               std::to_string(v2.z) + ")\n";

    }

//    void Triangle::GetDiffVariables(std::vector<Float const *> &vars) const { // FIXME
//        // Loop over all three vertices and get variables
//        for (unsigned i = 0; i < 3; ++i) {
//            // Get reference to vertex
//            Vector3F const &v = mesh.vertices[mesh.triangles[triangle_index].v[i]];
//            // Push pointers to variables
//            vars.push_back(&(v.x));
//            vars.push_back(&(v.y));
//            vars.push_back(&(v.z));
//        }
//    }
//
//    size_t Triangle::GetNumVars() const noexcept { // FIXME This is wrong!!!
//        return 3 * 3; // 3 Floats per 3 vertex
//    }
//
//    void
//    Triangle::UpdateDiffVariables(std::vector<float> const &delta, size_t starting_index) { // FIXME This is wrong!!!
//        for (unsigned i = 0; i < 3; ++i) {
//            // Get reference to vertex
//            Vector3F &v = mesh.vertices[mesh.triangles[triangle_index].v[i]];
//            // Update vertex values
//            v.x.SetValue(v.x.GetValue() + delta[starting_index + 3 * i]);
//            v.y.SetValue(v.y.GetValue() + delta[starting_index + 3 * i + 1]);
//            v.z.SetValue(v.z.GetValue() + delta[starting_index + 3 * i + 2]);
//        }
//    }
//
//    void Triangle::SetDiffVariables(std::vector<float> const &vals, size_t starting_index) {
//        for (unsigned i = 0; i < 3; ++i) {
//            // Get reference to vertex
//            Vector3F &v = mesh.vertices[mesh.triangles[triangle_index].v[i]];
//            // Update vertex values
//            v.x.SetValue(vals[starting_index + 3 * i]);
//            v.y.SetValue(vals[starting_index + 3 * i + 1]);
//            v.z.SetValue(vals[starting_index + 3 * i + 2]);
//        }
//    }

    void TriangleMesh::CreateTriangles() {
        for (size_t t = 0; t < triangles.size(); ++t) {
            shapes.push_back(std::make_shared<const Triangle>(*this, t));
        }
    }

    TriangleMesh::TriangleMesh(std::string const &file_name)
            : bvh(shapes) {
        // Loaded file lines
        std::vector<std::string> loaded_file;
        if (ReadFile(file_name, loaded_file)) {
            // Loop over all lines loaded and read the input data
            for (auto &line_it : loaded_file) {
                // Check if line is a comment
                if (line_it.compare(0, 1, "#") == 0) {
                    continue;
                }
                if (line_it.compare(0, 2, "v ") == 0) {
                    // Read vertex data
                    float x, y, z;
                    sscanf(line_it.c_str(), "v %f %f %f", &x, &y, &z);
                    vertices.push_back(Vector3f(x, y, z));
                } /* else if (line_it.compare(0, 3, "vt ") == 0) {
                    // Read uv data
                    float u, v;
                    sscanf(line_it.c_str(), "vt %f %f", &u, &v);
                    uvs.push_back(u);
                    uvs.push_back(v);
                } */ else if (line_it.compare(0, 3, "vn ") == 0) {
                    // Read normal data
                    float x, y, z;
                    sscanf(line_it.c_str(), "vn %f %f %f", &x, &y, &z);
                    normals.push_back(Normalize(Vector3f(x, y, z)));
                } else if (line_it.compare(0, 1, "f") == 0) {
                    // Read face data
                    uint32_t v0, v1, v2;
                    uint32_t n0, n1, n2;
                    uint32_t uv0, uv1, uv2;

                    // Check if face is just composed of vertices
                    if (sscanf(line_it.c_str(), "f %u %u %u", &v0, &v1, &v2) == 3) {
                        // triangles.push_back(MeshTriangle(v0 - 1, v1 - 1, v2 - 1));
                        triangles.push_back(TriangleIndices(v0 - 1, v1 - 1, v2 - 1));
                        // Check if face contains vertices and textures coordinates
                    } else if (sscanf(line_it.c_str(), "f %u/%u %u/%u %u/%u", &v0, &uv0, &v1, &uv1, &v2, &uv2) == 6) {
                        triangles.push_back(TriangleIndices(v0 - 1, v1 - 1, v2 - 1));
                        // Check if face contains textures coordinate, vertices and normals
                    } else if (sscanf(line_it.c_str(), "f %u/%u/%u %u/%u/%u %u/%u/%u", &v0, &uv0, &n0, &v1, &uv1, &n1,
                                      &v2, &uv2, &n2) == 9) {
                        triangles.push_back(TriangleIndices(v0 - 1, v1 - 1, v2 - 1, n0 - 1, n1 - 1, n2 - 1));
                        // Check if face contains normals and vertices only
                    } else if (sscanf(line_it.c_str(), "f %u//%u %u//%u %u//%u", &v0, &n0, &v1, &n1, &v2, &n2) == 6) {
                        triangles.push_back(TriangleIndices(v0 - 1, v1 - 1, v2 - 1, n0 - 1, n1 - 1, n2 - 1));
                    } else {
                        std::cerr << "Can not read face format!" << std::endl;
                    }
                } else {
                    std::cerr << "Unrecognized token during .obj parsing" << std::endl;
                }
            }

            // Clear loaded file
            loaded_file.clear();

            // Fit vectors
            triangles.shrink_to_fit();
            vertices.shrink_to_fit();
            normals.shrink_to_fit();

            // Print out some data about the loaded mesh
            std::cout << "Loaded triangle mesh with " << triangles.size() << " triangles and "
                      << vertices.size() << " vertices." << std::endl;

            // Create triangles
            CreateTriangles();

            // Build BVH
            bvh.Build();
        }
    }

    bool TriangleMesh::Intersect(Ray const &ray, Interaction *const interaction) const {
        return bvh.Intersect(ray, interaction);
    }

    bool TriangleMesh::IntersectP(Ray const &ray) const {
        return bvh.IntersectP(ray);
    }

    BBOX TriangleMesh::BBox() const {
        // Return BBOX from BVH acceleration structure
        return bvh.BBox();
    }

    Vector3f TriangleMesh::Centroid() const {
        return drdemo::Vector3f(0.f, 0.f, 0.f);
    }

    std::string TriangleMesh::ToString() const {
        return std::string(); // TODO
    }
//
//    void TriangleMesh::GetDiffVariables(std::vector<Float const *> &vars) const {
//        // Add pointer to all components of each vertex
//        for (auto const &v : vertices) {
//            vars.push_back(&(v.x));
//            vars.push_back(&(v.y));
//            vars.push_back(&(v.z));
//        }
//    }
//
//    size_t TriangleMesh::GetNumVars() const noexcept {
//        return vertices.size() * 3;
//    }
//
//    void TriangleMesh::UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) {
//        size_t used_vars = 0;
//        for (auto &v : vertices) {
//            // Update vertex values
//            v.x.SetValue(v.x.GetValue() + delta[starting_index + used_vars++]);
//            v.y.SetValue(v.y.GetValue() + delta[starting_index + used_vars++]);
//            v.z.SetValue(v.z.GetValue() + delta[starting_index + used_vars++]);
//        }
//
//        // Rebuild BVH acceleration structure, not a good solution but it's the best at the moment
//        bvh.Rebuild();
//    }
//
//    void TriangleMesh::SetDiffVariables(const std::vector<float> &vals, size_t starting_index) {
//        size_t used_vars = 0;
//        for (auto &v : vertices) {
//            // Update vertex values
//            v.x.SetValue(vals[starting_index + used_vars++]);
//            v.y.SetValue(vals[starting_index + used_vars++]);
//            v.z.SetValue(vals[starting_index + used_vars++]);
//        }
//
//        // Rebuild BVH acceleration structure, not a good solution but it's the best at the moment
//        bvh.Rebuild();
//    }

} // drdemo namespace
