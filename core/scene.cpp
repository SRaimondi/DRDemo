//
// Created by simon on 09.05.17.
//

#include "scene.hpp"

namespace drdemo {

    Scene::Scene() {}

    void Scene::AddShape(std::shared_ptr<Shape> const &shape) {
        shapes.push_back(shape);
    }

    void Scene::AddLight(std::shared_ptr<LightInterface> const &light) {
        lights.push_back(light);
    }

    bool Scene::Intersect(Ray const &ray, Interaction *const interaction) const {
        // Hit flag
        bool hit = false;
        // Loop over all shapes and look for closest interaction
        for (auto shape : shapes) {
            if (shape->Intersect(ray, interaction)) {
                hit = true;
            }
        }

        return hit;
    }

    bool Scene::IntersectP(Ray const &ray) const {
        for (auto shape : shapes) {
            if (shape->IntersectP(ray)) { return true; }
        }

        return false;
    }

} // drdemo namespace