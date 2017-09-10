//
// Created by Simon on 29.05.2017.
//

#include "bbox.hpp"

namespace drdemo {

    void BBOX::ExpandTo(const Vector3f &p) {
        bounds[0] = Min(bounds[0], p);
        bounds[1] = Max(bounds[1], p);
    }

    void BBOX::ExpandTo(const BBOX &b) {
        bounds[0] = Min(bounds[0], b.bounds[0]);
        bounds[1] = Max(bounds[1], b.bounds[1]);
    }

    void BBOX::Expand(float s) {
        bounds[0] *= s;
        bounds[1] *= s;
    }

    float BBOX::Distance(const Vector3f &p) const {
//        // Compute extent
//        const Vector3f extent = Extent();
//        const Vector3f d = Vector3f(std::abs(p.x) - extent.x,
//                                    std::abs(p.y) - extent.y,
//                                    std::abs(p.z) - extent.z);
//
//        return Length(Vector3f(std::max(d.x, 0.f), std::max(d.y, 0.f), std::max(d.z, 0.f))) +
//               std::min(std::max(d.x, std::max(d.y, d.z)), 0.f);
        float dist = 0.f;

        for (int i = 0; i < 3; i++) {
            // For each axis count any excess distance outside box extents
            float v = p[i];
            if (v < bounds[0][i]) { dist += (bounds[0][i] - v) * (bounds[0][i] - v); }
            if (v > bounds[1][i]) { dist += (v - bounds[1][i]) * (v - bounds[1][i]); }
        }

        return std::sqrt(dist);
    }

    unsigned BBOX::MaxDimension() const {
        Vector3f extent = Extent();
        unsigned dim = 0;
        if (extent.y > extent.x) { dim = 1; }
        if (extent.z > extent.y) { dim = 2; }

        return dim;
    }

    float BBOX::Surface() const {
        const Vector3f extent = Extent();

        return 2.f * (extent.x * extent.y + extent.x * extent.z + extent.y * extent.z);
    }

    bool BBOX::Intersect(const Ray &ray, float *const t_min, float *const t_max) const {
        // Find minimum intersection
        float tx_min = (bounds[ray.sign[0]].x - ray.o.x.GetValue()) / ray.d.x.GetValue();
        float ty_min = (bounds[ray.sign[1]].y - ray.o.y.GetValue()) / ray.d.y.GetValue();
        float tz_min = (bounds[ray.sign[2]].z - ray.o.z.GetValue()) / ray.d.z.GetValue();
        // Find maximum intersection
        float tx_max = (bounds[1 - ray.sign[0]].x - ray.o.x.GetValue()) / ray.d.x.GetValue();
        float ty_max = (bounds[1 - ray.sign[1]].y - ray.o.y.GetValue()) / ray.d.y.GetValue();
        float tz_max = (bounds[1 - ray.sign[2]].z - ray.o.z.GetValue()) / ray.d.z.GetValue();

        *t_min = std::max(tz_min, std::max(ty_min, std::max(tx_min, ray.t_min.GetValue())));
        *t_max = std::min(tz_max, std::min(ty_max, std::min(tx_max, ray.t_max.GetValue())));

        return *t_min <= *t_max;
    }

} // drdemo namespace