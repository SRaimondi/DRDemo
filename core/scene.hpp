//
// Created by simon on 09.05.17.
//

#ifndef DRDEMO_SCENE_HPP
#define DRDEMO_SCENE_HPP

#include <memory>
#include "shape.hpp"
#include "light.hpp"

namespace drdemo {

    /**
     * Define simple Scene class
     */
    class Scene {
    protected:
        // List of Shape in the scene
        std::vector<std::shared_ptr<Shape> > shapes;
        // List of Lights in the scene
        std::vector<std::shared_ptr<LightInterface> > lights;

    public:
        Scene();

        // Add Shape to scene
        void AddShape(std::shared_ptr<Shape> const &shape);

        // Add Light to scene
        void AddLight(std::shared_ptr<LightInterface> const &light);

        // Access list of shapes
        std::vector<std::shared_ptr<Shape> > const &GetShapes() const { return shapes; }

        std::vector<std::shared_ptr<Shape> > &GetShapes() { return shapes; }

        // Clear list of shapes
        void ClearShapes() { shapes.clear(); }

        // Access list of lights
        std::vector<std::shared_ptr<LightInterface> > const &GetLights() const { return lights; }

        std::vector<std::shared_ptr<LightInterface> > &GetLights() { return lights; }

        // Clear list of lights
        void ClearLights() { lights.clear(); }

        // Intersect Ray with Scene
        bool Intersect(Ray const &ray, Interaction *const interaction) const;

        bool IntersectP(Ray const &ray) const;
    };

} // drdemo namespace

#endif //DRDEMO_SCENE_HPP
