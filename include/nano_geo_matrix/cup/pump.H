#pragma once

#include <cmath>
#include <cstring>

enum class PumpMode { Step, Sin };

inline PumpMode pump_mode(const char* s) {
    if (!s) return PumpMode::Step;
    if (std::strcmp(s, "step") == 0) return PumpMode::Step;
    if (std::strcmp(s, "sin")  == 0) return PumpMode::Sin;
    return PumpMode::Step; // safe fallback
}

/*
 * Returns target inversion tildeN(t)
 *
 * t        : current time (same units as solver loop)
 * tpump    : pump start time
 * mode     : Step or Sin
 * tNmin    : minimum inversion
 * tNmax    : maximum inversion
 * omega    : angular frequency for Sin mode
 * phase    : phase for Sin mode
 */
inline double gimme_tildeN(
    double t,
    double tpump,
    PumpMode mode,
    double tNmin = 0.0,
    double tNmax = 1.0,
    double omega = 0.0,
    double phase = 0.0
) {

    if (tNmax < tNmin) {
        double tmp = tNmax;
        tNmax = tNmin;
        tNmin = tmp;
    }

    switch (mode) {

        case PumpMode::Step:
            return (t >= tpump) ? tNmax : tNmin;

        case PumpMode::Sin: {
            if (t < tpump)
                return tNmin;

            double mid  = 0.5 * (tNmax + tNmin);
            double half = 0.5 * (tNmax - tNmin);

            double val = mid + half * std::sin(omega * (t - tpump) + phase);

            // strict clamping (numerical safety)
            if (val > tNmax) val = tNmax;
            if (val < tNmin) val = tNmin;

            return val;
        }
    }

    return tNmin; // unreachable but keeps compiler calm
}
