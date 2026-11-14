#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>

namespace nano_geo_matrix::saturation {

// Overrides select an authority, rather than silently accepting inconsistent units.
enum class Tau1Policy { require_match, energy_inputs_authoritative, supplied_tau1_authoritative };

inline const char* policy_name(Tau1Policy policy) {
    switch (policy) {
        case Tau1Policy::require_match: return "require_match";
        case Tau1Policy::energy_inputs_authoritative: return "energy_inputs_authoritative";
        case Tau1Policy::supplied_tau1_authoritative: return "supplied_tau1_authoritative";
    }
    throw std::invalid_argument("Unknown saturation tau1 policy");
}

struct Inputs {
    double plasma_energy_eV;
    double linewidth_eV; // hbar*Delta, with Delta=2/tau2_SI
    double gain;        // G=Gmax for prescribed inversion one; magnitude is used
    double tau1_over_tau2 = 5.;
    std::optional<double> supplied_tau1_normalized = std::nullopt;
    Tau1Policy tau1_policy = Tau1Policy::require_match;
};

struct Scale {
    static constexpr const char* convention = "inversion_one_plasma_scaled_saturation_v1";
    static constexpr double tau1_match_relative_tolerance = 1.e-12;
    Inputs inputs; // Original values and explicit policy retained for provenance.
    double energy_derived_tau1_normalized;
    double tau1_normalized; // Effective value selected by the recorded policy.
    double tau2_normalized;
    double effective_tau1_over_tau2;
    double Esat_normalized;
    double Isat_normalized; // Squared electric field, not SI irradiance.
    bool supplied_tau1_matches; // True when no supplied value needs checking.
    bool mismatch_overridden;
};

namespace detail {
inline void positive(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.)
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
}
inline double representable(long double value, const char* name) {
    const double result = static_cast<double>(value);
    if (!std::isfinite(result) || result <= 0.)
        throw std::overflow_error(std::string(name) + " is outside the positive finite double range");
    return result;
}
} // namespace detail

// Pure scale calculation: does not configure a solver or infer input units.
inline Scale calculate(const Inputs& in) {
    detail::positive(in.plasma_energy_eV, "plasma_energy_eV");
    detail::positive(in.linewidth_eV, "linewidth_eV");
    detail::positive(std::abs(in.gain), "abs(gain)");
    detail::positive(in.tau1_over_tau2, "tau1_over_tau2");
    policy_name(in.tau1_policy); // Validate even if no comparison is necessary.
    if (in.supplied_tau1_normalized)
        detail::positive(*in.supplied_tau1_normalized, "supplied_tau1_normalized");
    if (in.tau1_policy != Tau1Policy::require_match && !in.supplied_tau1_normalized)
        throw std::invalid_argument("An explicit tau1 override requires supplied_tau1_normalized");

    const long double tau2 = 2.L * in.plasma_energy_eV / in.linewidth_eV;
    const long double derived_tau1 = in.tau1_over_tau2 * tau2;
    const double derived = detail::representable(derived_tau1, "energy-derived tau1_normalized");
    const double tau2_result = detail::representable(tau2, "tau2_normalized");
    bool matches = true;
    if (in.supplied_tau1_normalized) {
        const long double supplied = *in.supplied_tau1_normalized;
        matches = std::abs(supplied-derived_tau1) <= Scale::tau1_match_relative_tolerance *
            std::max(supplied, derived_tau1);
        if (!matches && in.tau1_policy == Tau1Policy::require_match)
            throw std::invalid_argument("Supplied normalized tau1 disagrees with energy inputs and tau1_over_tau2; select an explicit authority to override");
    }
    const long double effective_tau1 = in.tau1_policy == Tau1Policy::supplied_tau1_authoritative
        ? static_cast<long double>(*in.supplied_tau1_normalized) : derived_tau1;
    const long double intensity = 1.L / (std::abs(static_cast<long double>(in.gain)) * effective_tau1);
    return {in, derived, detail::representable(effective_tau1, "tau1_normalized"), tau2_result,
        detail::representable(effective_tau1/tau2, "effective_tau1_over_tau2"),
        detail::representable(std::sqrt(intensity), "Esat_normalized"),
        detail::representable(intensity, "Isat_normalized"), matches, !matches};
}

} // namespace nano_geo_matrix::saturation
