//
// Created by Simon on 28.08.2017.
//

#ifndef DRDEMO_AMBIENT_LIGHT_HPP
#define DRDEMO_AMBIENT_LIGHT_HPP

#include "light.hpp"
#include "diff_object.hpp"

namespace drdemo {

    /**
     * Define ambient light class
     */
    class AmbientLight : public LightInterface, public DiffObjectInterface {
        // private:
    public:
        // Ambient intensity
        Spectrum intensity;

        // public:
        explicit AmbientLight(const Spectrum &i);

        Spectrum SampleLi(const Interaction &interaction, float u0, float u1,
                          Vector3F *wi, Float *pdf) const override;

        void GetDiffVariables(std::vector<Float const *> &vars) const override;

        size_t GetNumVars() const noexcept override;

        void UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) override;

        void SetDiffVariables(const std::vector<float> &vals, size_t starting_index) override;
    };

} // drdemo namespace

#endif //DRDEMO_AMBIENT_LIGHT_HPP
