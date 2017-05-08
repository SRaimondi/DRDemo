//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_SCENE_HPP
#define DRDEMO_SCENE_HPP

#include <bits/shared_ptr.h>
#include "shape.hpp"

namespace drdemo {

    /**
     * Define simple Scene class
     */
    class Scene {
    protected:
        // List of Shape in the scene
        std::vector<std::shared_ptr<Shape> > shapes;

    public:
        Scene();

        // Add shapeto scene
        void AddShape(std::shared_ptr<Shape> const &shape);

        // Intersect Ray with Scene
        bool Intersect(Ray const &ray, Interaction *const interaction) const;

        bool IntersectP(Ray const &ray) const;
    };

} // drdemo namespace

#endif //DRDEMO_SCENE_HPP
