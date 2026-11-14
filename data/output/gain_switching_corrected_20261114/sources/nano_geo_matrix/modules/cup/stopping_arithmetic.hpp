#pragma once

#include <cmath>

namespace cup_stopping_detail {

inline bool finite_positive(double value)
{
    return std::isfinite(value) && value > 0.;
}

// Legacy Fourier sample gating only; not a sustained convergence criterion.
inline bool relative_change_below(double current, double previous, double tolerance)
{
    return finite_positive(current) && std::isfinite(previous) &&
           std::fabs((current - previous)/current) < tolerance;
}

inline bool amplitude_above(double intensity, double field, double threshold)
{
    return finite_positive(field) && intensity/field > threshold;
}

} // namespace cup_stopping_detail
