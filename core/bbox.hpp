//
// Created by Simon on 29.05.2017.
//

#ifndef DRDEMO_BBOX_HPP
#define DRDEMO_BBOX_HPP

#include "geometry.hpp"

namespace drdemo {

    /**
     * Define BBOX class
     */
    class BBOX {
    private:
        // BBOX bounds
        Vector3f bounds[2];

    public:
        // Constructor
        inline BBOX() {
            bounds[0] = Vector3f(INFINITY, INFINITY, INFINITY);
            bounds[1] = Vector3f(-INFINITY, -INFINITY, -INFINITY);
        }

        inline explicit BBOX(Vector3f const &p) {
            bounds[0] = p;
            bounds[1] = p;
        }

        inline BBOX(Vector3f const &p1, Vector3f const &p2) {
            bounds[0] = Min(p1, p2);
            bounds[1] = Max(p1, p2);
        }

        // Access minimum and maximum
        inline Vector3f const &MinPoint() const { return bounds[0]; }

        inline Vector3f const &MaxPoint() const { return bounds[1]; }

        // BBOX utility methods
        void ExpandTo(const Vector3f &p);

        void ExpandTo(const BBOX &b);

        // Check if a point is inside
        inline bool Inside(const Vector3f &p) const {
            return p.x >= bounds[0].x && p.x <= bounds[1].x &&
                   p.y >= bounds[0].y && p.y <= bounds[1].y &&
                   p.z >= bounds[0].z && p.z <= bounds[1].z;
        }

        // Compute the distance of a point to the BBOX
        float Distance(const Vector3f &p) const;

        // Compute BBOX extent
        inline Vector3f Extent() const {
            return bounds[1] - bounds[0];
        }

        // Return axis with maximum dimension
        unsigned MaxDimension() const;

        // Surface area
        float Surface() const;

        // Compute intersection with ray
        bool Intersect(const Ray &ray, float *t_min, float *t_max) const;
    };

} // drdemo namespace

#endif //DRDEMO_BBOX_HPP
