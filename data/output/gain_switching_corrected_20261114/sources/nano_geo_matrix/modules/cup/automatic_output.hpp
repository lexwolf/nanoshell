#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace cup_output_detail {

struct Cadence {
    int steps = 1;
    double actual_ps = 0.;
};

// Round upward: the saved interval is never shorter than requested, and is
// at least one integration step. frequency is Ep/hbar, in radians/second.
inline Cadence automatic_cadence(double requested_ps, double frequency, double dt)
{
    if (!std::isfinite(requested_ps) || requested_ps <= 0. ||
        !std::isfinite(frequency) || frequency <= 0. ||
        !std::isfinite(dt) || dt <= 0.) {
        throw std::invalid_argument("Automatic output interval, frequency and step must be finite and positive");
    }
    const long double step_ps = static_cast<long double>(dt)*1.e12L/frequency;
    const long double limit = std::numeric_limits<int>::max();
    // Check before dividing or converting to int, including enormous requests.
    if (!std::isfinite(step_ps) || step_ps <= 0. || requested_ps/limit > step_ps) {
        throw std::overflow_error("Automatic output interval exceeds the step-count range");
    }
    const long double rounded = std::max(1.L, std::ceil(requested_ps/step_ps));
    if (rounded > limit || rounded*step_ps > std::numeric_limits<double>::max()) {
        throw std::overflow_error("Automatic output interval is not representable");
    }
    return {static_cast<int>(rounded), static_cast<double>(rounded*step_ps)};
}

// completed_step is the timestamp index (i-1), rather than the loop's i after
// increment. The first scheduled automatic row is at steps*dt, not at zero.
inline bool scheduled(int completed_step, int steps)
{
    return completed_step > 0 && completed_step % steps == 0;
}

// Preserve distinguishable automatic timestamps without changing value columns.
inline std::ostream& write_timestamp(std::ostream& file, double time_ps, bool automatic)
{
    const auto precision = file.precision();
    if (automatic) file.precision(std::numeric_limits<double>::max_digits10);
    file << "  " << std::setw(8) << std::setiosflags(std::ios::left) << time_ps;
    file.precision(precision);
    return file;
}

inline void write_metadata(std::initializer_list<std::ostream*> files,
                           double requested_ps, const Cadence& cadence)
{
    std::ostringstream text;
    text << std::setprecision(17)
         << "# requested_output_interval_ps = " << requested_ps << '\n'
         << "# actual_output_interval_ps = " << cadence.actual_ps << '\n'
         << "# output_interval_steps = " << cadence.steps << '\n'
         << "# output_rounding = ceil_to_steps_minimum_one\n";
    for (auto* file : files) *file << text.str();
}

inline void write_termination(std::initializer_list<std::ostream*> files, const char* reason)
{
    for (auto* file : files) *file << "# termination_reason = " << reason << '\n';
}

} // namespace cup_output_detail
