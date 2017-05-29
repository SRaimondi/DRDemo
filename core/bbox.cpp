//
// Created by Simon on 29.05.2017.
//

#include "bbox.hpp"

namespace drdemo {

    void BBOX::ExpandTo(Vector3F const &p) {
        bounds[0] = Min(bounds[0], p);
        bounds[1] = Max(bounds[1], p);
    }

    void BBOX::ExpandTo(BBOX const &b) {
        bounds[0] = Min(bounds[0], b.bounds[0]);
        bounds[1] = Min(bounds[1], b.bounds[1]);
    }

    unsigned BBOX::MaxDimension() const {
        Vector3F extent = Extent();
        unsigned dim = 0;
        if (extent.y > extent.x) { dim = 1; }
        if (extent.z > extent.y) { dim = 2; }

        return dim;
    }

    Float BBOX::Surface() const {
        Vector3F extent = Extent();

        return 2.f * (extent.x * extent.y + extent.x * extent.z + extent.y * extent.z);
    }

    bool BBOX::Intersect(Ray const &ray, Float *const t_min, Float *const t_max) const {
        // Find minimum intersection
        Float tx_min = (bounds[ray.sign[0]].x - ray.o.x) / ray.d.x;
        Float ty_min = (bounds[ray.sign[1]].y - ray.o.y) / ray.d.y;
        Float tz_min = (bounds[ray.sign[2]].z - ray.o.z) / ray.d.z;
        // Find maximum intersection
        Float tx_max = (bounds[1 - ray.sign[0]].x - ray.o.x) / ray.d.x;
        Float ty_max = (bounds[1 - ray.sign[1]].y - ray.o.y) / ray.d.y;
        Float tz_max = (bounds[1 - ray.sign[2]].z - ray.o.z) / ray.d.z;

        *t_min = Max(tz_min, Max(ty_min, Max(tx_min, ray.t_min)));
        *t_max = Min(tz_max, Min(ty_max, Min(tx_max, ray.t_max)));

        return *t_min <= *t_max;
    }

} // drdemo namespace