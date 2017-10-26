//
// Created by Simon on 28.08.2017.
//

#include "ambient_light.hpp"

namespace drdemo {

    AmbientLight::AmbientLight(const Spectrum &i) : intensity(i) {}

    Spectrum
    AmbientLight::SampleLi(Interaction const &interaction, float u0, float u1, Vector3F *wi, Float *pdf) const {
        // Set as interaction normal
        wi->x = interaction.n.x.GetValue();
        wi->y = interaction.n.y.GetValue();
        wi->z = interaction.n.z.GetValue();

        // Set pdf
        *pdf = 1.f;

        return intensity;
    }

    void AmbientLight::GetDiffVariables(std::vector<Float const *> &vars) const {
        // Add RGB values as differentiable variables
        vars.push_back(&(intensity.r));
        vars.push_back(&(intensity.g));
        vars.push_back(&(intensity.b));
    }

    size_t AmbientLight::GetNumVars() const noexcept {
        return 3;
    }

    void AmbientLight::UpdateDiffVariables(const std::vector<float> &delta, size_t starting_index) {
        intensity.r.SetValue(intensity.r.GetValue() + delta[starting_index]);
        intensity.g.SetValue(intensity.g.GetValue() + delta[starting_index + 1]);
        intensity.b.SetValue(intensity.b.GetValue() + delta[starting_index + 2]);
    }

    void AmbientLight::SetDiffVariables(const std::vector<float> &vals, size_t starting_index) {
        intensity.r.SetValue(vals[starting_index]);
        intensity.g.SetValue(vals[starting_index + 1]);
        intensity.b.SetValue(vals[starting_index + 2]);
    }

} // drdemo namespace