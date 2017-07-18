//
// Created by Simon on 29.05.2017.
//

#ifndef DRDEMO_BBOX_HPP
#define DRDEMO_BBOX_HPP

#include "geometry.hpp"

namespace drdemo {

    /**
     * Define BBOX class
     * TODO: this class maybe could work with simple float and not need the Float class
     */
    class BBOX {
    private:
        // Vector3F bounds[2];
        Vector3f bounds[2];

    public:
        // Constructor
        inline BBOX() {
            bounds[0] = Vector3f(INFINITY, INFINITY, INFINITY);
            bounds[1] = Vector3f(-INFINITY, -INFINITY, -INFINITY);
        }

        inline BBOX(Vector3f const &p) {
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
        void ExpandTo(Vector3f const &p);

        void ExpandTo(BBOX const &b);

        // Compute BBOX extent
        inline Vector3f Extent() const {
            return bounds[1] - bounds[0];
        }

        // Return axis with maximum dimension
        unsigned MaxDimension() const;

        // Surface area
        float Surface() const;

        // Compute intersection with ray
        bool Intersect(Ray const &ray, float *const t_min, float *const t_max) const;
    };

} // drdemo namespace

#endif //DRDEMO_BBOX_HPP
