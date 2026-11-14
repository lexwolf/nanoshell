#pragma once

#include <functional>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <optional>
#include <stdexcept>

namespace cup_convergence {

// Times are ps; absolute_flatness is in the solver's observable intensity units.
struct Options {
    double flat_window_ps = 10.;
    int flat_windows = 3;
    int flat_confirmation_windows = 1;
    int periodic_cycles = 5;
    double periodic_min_span_ps = 30.;
    int periodic_confirmation_cycles = 2;
    double relative_flatness = 1.e-4;
    double absolute_flatness = 0.;
    double relative_period = .01;
    double relative_amplitude = .01;
    double relative_cycle_mean = .01;
    // Hard storage bound; periods requiring more cycles to span the minimum
    // observation remain unclassified. Independent of saved-output cadence.
    int max_periodic_cycles = 4096;

    void validate() const {
        const auto positive = [](double x) { return std::isfinite(x) && x > 0.; };
        const auto nonnegative = [](double x) { return std::isfinite(x) && x >= 0.; };
        if (!positive(flat_window_ps) || !positive(periodic_min_span_ps) ||
            !nonnegative(relative_flatness) || !nonnegative(absolute_flatness) ||
            (relative_flatness == 0. && absolute_flatness == 0.) ||
            !positive(relative_period) || !positive(relative_amplitude) ||
            !positive(relative_cycle_mean) || flat_windows < 1 ||
            flat_confirmation_windows < 0 || flat_windows > 4096 ||
            flat_confirmation_windows > 4096-flat_windows || periodic_cycles < 2 ||
            periodic_confirmation_cycles < 0 || max_periodic_cycles < periodic_cycles ||
            max_periodic_cycles > 65536 ||
            periodic_confirmation_cycles > max_periodic_cycles-periodic_cycles)
            throw std::invalid_argument("Invalid sustained-convergence options");
    }
};

struct CycleObservation { double start_ps, end_ps, mean, minimum, maximum; };

struct Result {
    const char* classification = "unclassified";
    double start_ps = 0., qualifying_end_ps = 0., confirmation_end_ps = 0., detected_at_ps = 0.;
    double mean_intensity = 0., peak_to_peak = 0., mean_spread = 0.;
    double mean_period_ps = 0., period_spread_ps = 0., mean_amplitude = 0., amplitude_spread = 0.;
    double maxima_spread = 0., minima_spread = 0.;
    int windows_or_cycles = 0;
};

class Monitor {
    std::function<void(const CycleObservation&)> cycle_observer_;
    using Real = long double;
    struct Point { double t = 0., y = 0.; Real area = 0.; };
    struct Window { double start, end; Real low, high, mean; };
    struct Cycle { double start, end; Real low, high, mean, amplitude, start_peak; };
    Options opt_;
    bool periodic_, have_sample_ = false, flat_candidate_ = false, cycle_candidate_ = false;
    Point previous_, peak_, valley_, last_peak_;
    bool have_peak_ = false;
    int direction_ = 0, flat_confirmed_ = 0, cycle_confirmed_ = 0;
    Real integral_ = 0., window_area_ = 0., window_low_ = 0., window_high_ = 0.;
    double window_start_ = 0., flat_start_ = 0., flat_end_ = 0., cycle_start_ = 0., cycle_end_ = 0.;
    std::deque<Window> windows_;
    std::deque<Cycle> cycles_;
    Result result_;
    bool accepted_ = false;

    Real resolution(Real mean) const { return opt_.absolute_flatness + Real(opt_.relative_flatness)*mean; }
    static Real range(Real low, Real high) { return high-low; }
    void start(const Point& p) {
        have_sample_ = true; previous_ = peak_ = valley_ = p;
        window_start_ = p.t; window_low_ = window_high_ = p.y;
    }
    bool flat_statistics(Result& r) const {
        if (windows_.empty()) return false;
        Real low = windows_.front().low, high = windows_.front().high;
        Real ml = windows_.front().mean, mh = ml, area = 0., span = 0.;
        for (const auto& w : windows_) {
            low = std::min(low, w.low); high = std::max(high, w.high);
            ml = std::min(ml, w.mean); mh = std::max(mh, w.mean);
            area += w.mean*(w.end-w.start); span += w.end-w.start;
        }
        if (!(span > 0.)) return false;
        const Real mean = area/span, tolerance = resolution(mean);
        if (!std::isfinite(mean) || range(low, high) > tolerance || range(ml, mh) > tolerance) return false;
        for (const auto& w : windows_) if (w.high-w.low > tolerance) return false;
        r.mean_intensity = double(mean); r.peak_to_peak = double(high-low);
        r.mean_spread = double(mh-ml); r.windows_or_cycles = int(windows_.size());
        return true;
    }
    void complete_window(double end) {
        windows_.push_back({window_start_, end, window_low_, window_high_, window_area_/(end-window_start_)});
        if (!flat_candidate_) while (windows_.size() > std::size_t(opt_.flat_windows)) windows_.pop_front();
        Result stats;
        if (windows_.size() >= std::size_t(opt_.flat_windows) && flat_statistics(stats)) {
            if (!flat_candidate_) {
                flat_candidate_ = true; flat_confirmed_ = 0;
                flat_start_ = windows_.front().start; flat_end_ = end;
            } else ++flat_confirmed_;
            if (flat_confirmed_ >= opt_.flat_confirmation_windows) {
                result_ = stats; result_.classification = "sustained_flat";
                result_.start_ps = flat_start_; result_.qualifying_end_ps = flat_end_;
                result_.confirmation_end_ps = end; accepted_ = true;
            }
        } else {
            flat_candidate_ = false;
            while (windows_.size() > std::size_t(opt_.flat_windows)) windows_.pop_front();
        }
    }
    void flat_sample(const Point& p) {
        double left = previous_.t;
        Real value = previous_.y;
        while (left < p.t && !accepted_) {
            const double boundary = window_start_+opt_.flat_window_ps;
            if (!std::isfinite(boundary) || boundary <= window_start_) { reset(); return; }
            const double right = std::min(boundary, p.t);
            const Real next = right == p.t ? Real(p.y) : Real(previous_.y) +
                (Real(p.y)-previous_.y)*(Real(right)-previous_.t)/(Real(p.t)-previous_.t);
            window_area_ += (value+next)*.5L*(right-left);
            window_low_ = std::min({window_low_, value, next});
            window_high_ = std::max({window_high_, value, next});
            if (right == boundary) {
                complete_window(boundary);
                window_start_ = boundary; window_area_ = 0.;
                window_low_ = window_high_ = next;
            }
            left = right; value = next;
        }
    }
    bool periodic_statistics(Result& r) const {
        if (cycles_.size() < std::size_t(opt_.periodic_cycles)) return false;
        const double span = cycles_.back().end-cycles_.front().start;
        if (!std::isfinite(span) || span <= 0. || span < opt_.periodic_min_span_ps) return false;
        Real pl = std::numeric_limits<Real>::max(), ph = 0., al = pl, ah = 0.;
        Real xl = pl, xh = 0., nl = pl, nh = 0., ml = pl, mh = 0., area = 0., sum_a = 0.;
        for (const auto& c : cycles_) {
            const Real period = c.end-c.start;
            pl = std::min(pl, period); ph = std::max(ph, period);
            al = std::min(al, c.amplitude); ah = std::max(ah, c.amplitude);
            xl = std::min({xl, c.high, c.start_peak}); xh = std::max(xh, c.high);
            nl = std::min(nl, c.low); nh = std::max(nh, c.low);
            ml = std::min(ml, c.mean); mh = std::max(mh, c.mean);
            area += c.mean*period; sum_a += c.amplitude;
        }
        const Real mean = area/span, amplitude = sum_a/cycles_.size(), period = Real(span)/cycles_.size();
        if (!std::isfinite(mean) || !std::isfinite(amplitude) || amplitude <= resolution(mean) ||
            ph-pl > opt_.relative_period*period || ah-al > opt_.relative_amplitude*amplitude ||
            xh-xl > opt_.relative_amplitude*amplitude || nh-nl > opt_.relative_amplitude*amplitude ||
            mh-ml > opt_.absolute_flatness+Real(opt_.relative_cycle_mean)*mean) return false;
        r.mean_intensity = double(mean); r.peak_to_peak = double(xh-nl); r.mean_spread = double(mh-ml);
        r.mean_period_ps = double(period); r.period_spread_ps = double(ph-pl);
        r.mean_amplitude = double(amplitude); r.amplitude_spread = double(ah-al);
        r.maxima_spread = double(xh-xl); r.minima_spread = double(nh-nl);
        r.windows_or_cycles = int(cycles_.size());
        return true;
    }
    void maximum(const Point& p) {
        if (have_peak_) {
            const Real period = p.t-last_peak_.t;
            if (!(period > 0.)) return;
            const Real high = std::max(p.y, last_peak_.y), low = valley_.y;
            if (cycle_observer_) {
                cycle_observer_({last_peak_.t, p.t, double((p.area-last_peak_.area)/period), double(low), double(high)});
                last_peak_ = p; have_peak_ = true;
                return;
            }
            cycles_.push_back({last_peak_.t, p.t, low, high, (p.area-last_peak_.area)/period, high-low, last_peak_.y});
            if (cycles_.size() > std::size_t(opt_.max_periodic_cycles)) {
                cycles_.pop_front(); cycle_candidate_ = false;
            }
            if (!cycle_candidate_) {
                while (cycles_.size() > std::size_t(opt_.periodic_cycles) &&
                       cycles_.back().end-cycles_[1].start >= opt_.periodic_min_span_ps) cycles_.pop_front();
            }
            Result stats;
            if (periodic_statistics(stats)) {
                if (!cycle_candidate_) {
                    cycle_candidate_ = true; cycle_confirmed_ = 0;
                    cycle_start_ = cycles_.front().start; cycle_end_ = p.t;
                } else ++cycle_confirmed_;
                if (cycle_confirmed_ >= opt_.periodic_confirmation_cycles) {
                    result_ = stats; result_.classification = "sustained_periodic";
                    result_.start_ps = cycle_start_; result_.qualifying_end_ps = cycle_end_;
                    result_.confirmation_end_ps = p.t; accepted_ = true;
                }
            } else cycle_candidate_ = false;
        }
        last_peak_ = p; have_peak_ = true;
    }
    void periodic_sample(Point p) {
        // Hysteresis uses the configured flatness resolution at the local level.
        // Strict updates choose the FIRST sample of an exactly equal plateau.
        if (direction_ == 0) {
            if (p.y > peak_.y) peak_ = p;
            if (p.y < valley_.y) valley_ = p;
            if (Real(p.y)-valley_.y > resolution((Real(p.y)+valley_.y)*.5L)) { direction_ = 1; peak_ = p; }
            else if (Real(peak_.y)-p.y > resolution((Real(peak_.y)+p.y)*.5L)) { direction_ = -1; valley_ = p; }
        } else if (direction_ > 0) {
            if (p.y > peak_.y) peak_ = p;
            if (Real(peak_.y)-p.y > resolution((Real(peak_.y)+p.y)*.5L)) {
                // Rebase cycle integrals at every confirmed peak, avoiding
                // cancellation against a large or long-lived prehistory.
                const Real origin = peak_.area;
                maximum(peak_);
                integral_ -= origin; p.area -= origin; last_peak_.area = 0.;
                direction_ = -1; valley_ = p;
            }
        } else {
            if (p.y < valley_.y) valley_ = p;
            if (Real(p.y)-valley_.y > resolution((Real(p.y)+valley_.y)*.5L)) { direction_ = 1; peak_ = p; }
        }
    }
public:
    explicit Monitor(const Options& options = {}, bool periodic = true) : opt_(options), periodic_(periodic) { opt_.validate(); }
    void reset() {
        have_sample_ = flat_candidate_ = cycle_candidate_ = have_peak_ = accepted_ = false;
        direction_ = flat_confirmed_ = cycle_confirmed_ = 0;
        integral_ = window_area_ = 0.; windows_.clear(); cycles_.clear(); result_ = {};
    }
    bool update(double time_ps, double intensity) {
        if (!std::isfinite(time_ps) || !std::isfinite(intensity) || intensity < 0.) { reset(); return false; }
        if (have_sample_ && time_ps <= previous_.t) reset();
        if (accepted_) return true;
        Point p{time_ps, intensity, 0.};
        if (!have_sample_) { start(p); return false; }
        integral_ += (Real(previous_.y)+intensity)*.5L*(Real(time_ps)-previous_.t); p.area = integral_;
        flat_sample(p);
        if (!have_sample_) return false;
        if (periodic_ && !accepted_) periodic_sample(p);
        previous_ = p;
        if (accepted_) result_.detected_at_ps = time_ps;
        return accepted_;
    }
    // Observation only: reuse the same hysteresis/extrema and cycle integrals,
    // without flat/periodic acceptance or retaining candidate cycle records.
    void set_cycle_observer(std::function<void(const CycleObservation&)> observer) { cycle_observer_ = std::move(observer); }
    void observe_cycle_sample(double time_ps, double intensity) {
        if (!cycle_observer_) throw std::logic_error("Cycle observer is required");
        if (!std::isfinite(time_ps) || !std::isfinite(intensity) || intensity < 0. ||
            (have_sample_ && time_ps <= previous_.t)) { reset(); return; }
        Point p{time_ps, intensity, 0.};
        if (!have_sample_) { start(p); return; }
        integral_ += (Real(previous_.y)+intensity)*.5L*(Real(time_ps)-previous_.t);
        p.area = integral_;
        periodic_sample(p);
        previous_ = p;
    }
    const Result& result() const { return result_; }
    std::size_t stored_records() const { return windows_.size()+cycles_.size(); }
};

inline void write_options(std::initializer_list<std::ostream*> files, const Options& o, bool periodic, const char* pump = "step") {
    for (auto* file : files) {
        const auto precision = file->precision();
        *file << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "# convergence_flat_window_ps = " << o.flat_window_ps
              << "\n# convergence_flat_windows = " << o.flat_windows
              << "\n# convergence_flat_confirmation_windows = " << o.flat_confirmation_windows
              << "\n# convergence_periodic_cycles = " << o.periodic_cycles
              << "\n# convergence_periodic_min_span_ps = " << o.periodic_min_span_ps
              << "\n# convergence_periodic_confirmation_cycles = " << o.periodic_confirmation_cycles
              << "\n# convergence_relative_flatness = " << o.relative_flatness
              << "\n# convergence_absolute_flatness = " << o.absolute_flatness
              << "\n# convergence_relative_period = " << o.relative_period
              << "\n# convergence_relative_amplitude = " << o.relative_amplitude
              << "\n# convergence_relative_cycle_mean = " << o.relative_cycle_mean
              << "\n# convergence_max_periodic_cycles = " << o.max_periodic_cycles
              << "\n# convergence_periodic_enabled = " << periodic
              << "\n# convergence_pump_behavior = " << pump << '\n';
        file->precision(precision);
    }
}
inline void write_result(std::initializer_list<std::ostream*> files, const Result& r) {
    for (auto* file : files) {
        const auto precision = file->precision();
        *file << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "# convergence_classification = " << r.classification
              << "\n# convergence_start_ps = " << r.start_ps
              << "\n# convergence_qualifying_end_ps = " << r.qualifying_end_ps
              << "\n# convergence_confirmation_end_ps = " << r.confirmation_end_ps
              << "\n# convergence_detected_at_ps = " << r.detected_at_ps
              << "\n# convergence_mean_intensity = " << r.mean_intensity
              << "\n# convergence_peak_to_peak = " << r.peak_to_peak
              << "\n# convergence_mean_spread = " << r.mean_spread
              << "\n# convergence_mean_period_ps = " << r.mean_period_ps
              << "\n# convergence_period_spread_ps = " << r.period_spread_ps
              << "\n# convergence_mean_amplitude = " << r.mean_amplitude
              << "\n# convergence_amplitude_spread = " << r.amplitude_spread
              << "\n# convergence_maxima_spread = " << r.maxima_spread
              << "\n# convergence_minima_spread = " << r.minima_spread
              << "\n# convergence_windows_or_cycles = " << r.windows_or_cycles << '\n';
        file->precision(precision);
    }
}
} // namespace cup_convergence
