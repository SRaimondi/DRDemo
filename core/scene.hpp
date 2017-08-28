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
        Scene() = default;

        // Add Shape to scene
        void AddShape(std::shared_ptr<Shape> const &shape);

        // Add Light to scene
        void AddLight(std::shared_ptr<LightInterface> const &light);

        // Access list of shapes
        inline std::vector<std::shared_ptr<Shape> > const &GetShapes() const { return shapes; }

        // Clear list of shapes
        void ClearShapes() { shapes.clear(); }

        // Access list of lights
        inline std::vector<std::shared_ptr<LightInterface> > const &GetLights() const { return lights; }

        // Clear list of lights
        void ClearLights() { lights.clear(); }

        // Enable / disable all lights
        void EnableLights();

        void DisableLights();

        // Enable only light at given index
        void EnableLight(size_t index);

        // Intersect Ray with Scene
        bool Intersect(Ray const &ray, Interaction *interaction) const;

        bool IntersectP(Ray const &ray) const;
    };

} // drdemo namespace

#endif //DRDEMO_SCENE_HPP
