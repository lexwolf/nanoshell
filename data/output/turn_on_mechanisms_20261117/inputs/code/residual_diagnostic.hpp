#ifndef NANOSHELL_RESIDUAL_DIAGNOSTIC_HPP
#define NANOSHELL_RESIDUAL_DIAGNOSTIC_HPP
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

// Statistics only. Convergence and cycle identification remain in CUP.
struct ResidualDiagnostic {
    struct Window { double start, end, mean, low, high; int cycles; };
    double window_ps=40., range_rtol=.05, range_atol=1e-8, mean_rtol=1e-4;
    int window_count=3;
    std::vector<Window> windows;
    cup_convergence::Monitor cycles;
    bool started=false, periodic=false, done=false, invalid=false;
    double accepted_at=0., target_span=0., detected_period=0., previous_t=0., previous_y=0.;
    double start=0., low=0., high=0.;
    long double area=0.;
    int cycle_count=0, target_cycles=0;
    std::string status="unresolved";
    explicit ResidualDiagnostic(const cup_convergence::Options& options):cycles(options,true) {
        cycles.set_cycle_observer([this](const cup_convergence::CycleObservation& c) {
            if(done) return;
            if(cycle_count==0) { start=c.start_ps; low=c.minimum; high=c.maximum; area=0.; }
            low=std::min(low,c.minimum); high=std::max(high,c.maximum);
            area+=(long double)c.mean*(c.end_ps-c.start_ps); ++cycle_count;
            if(cycle_count>=target_cycles && c.end_ps-start>=target_span) {
                windows.push_back({start,c.end_ps,double(area/(c.end_ps-start)),low,high,cycle_count});
                cycle_count=0;
                if(int(windows.size())>=window_count) done=true;
            }
        });
    }
    bool update(double t, double y, const cup_convergence::Result& accepted) {
        if(!started) {
            started=true; accepted_at=t; detected_period=accepted.mean_period_ps;
            periodic=std::string(accepted.classification)=="sustained_periodic";
            target_span=periodic ? std::max(window_ps,5.*detected_period) : window_ps;
            const double count=periodic ? std::ceil(target_span/detected_period) : 0.;
            if(!std::isfinite(target_span) || count>std::numeric_limits<int>::max() || !std::isfinite(count)) { invalid=true; return true; }
            target_cycles=int(count);
            previous_t=start=t; previous_y=low=high=y; area=0.;
        }
        if(!std::isfinite(t)||!std::isfinite(y)||y<0.) { invalid=true; return true; }
        if(periodic) cycles.observe_cycle_sample(t,y);
        else if(t>previous_t) {
            double left_t=previous_t, left_y=previous_y;
            while(left_t<t && !done) {
                const double right_t=std::min(t,start+window_ps);
                const double right_y=previous_y+(y-previous_y)*(right_t-previous_t)/(t-previous_t);
                area+=(long double)(left_y+right_y)*.5L*(right_t-left_t);
                low=std::min({low,left_y,right_y}); high=std::max({high,left_y,right_y});
                if(right_t>=start+window_ps) {
                    windows.push_back({start,right_t,double(area/(right_t-start)),low,high,0});
                    start=right_t; low=high=right_y; area=0.;
                    if(int(windows.size())>=window_count) done=true;
                }
                left_t=right_t; left_y=right_y;
            }
        }
        previous_t=t; previous_y=y;
        // Missing/extinguished cycles cannot create a periodic zero estimate.
        if(!done && t-accepted_at>3.*window_count*target_span) { invalid=true; return true; }
        return done;
    }
    void write(const std::string& folder) {
        bool valid=done&&!invalid && int(windows.size())==window_count;
        double minr=std::numeric_limits<double>::infinity(),maxr=0.,sumr=0.;
        double minm=minr,maxm=0.,summ=0.;
        bool descending=true;
        for(std::size_t i=0;i<windows.size();++i) {
            auto w=windows[i];
            if(!(w.mean>0.)||!std::isfinite(w.mean)||!std::isfinite(w.high-w.low)) { valid=false; continue; }
            const double r=(w.high-w.low)/w.mean;
            minr=std::min(minr,r);maxr=std::max(maxr,r);sumr+=r;
            minm=std::min(minm,w.mean);maxm=std::max(maxm,w.mean);summ+=w.mean;
            if(i && r>(windows[i-1].high-windows[i-1].low)/windows[i-1].mean+range_atol) descending=false;
        }
        if(valid) {
            const double tolerance=range_atol+range_rtol*sumr/windows.size();
            const bool drift=maxm-minm>mean_rtol*summ/windows.size();
            if(descending && maxr-minr>tolerance) status=drift ? "decaying_variation_and_mean_drift" : "decaying_variation";
            else if(drift) status="drifting_mean";
            else if(maxr-minr<=tolerance) status="stable_interval_variation";
        }
        std::ofstream meta(folder+"/diagnostic_metadata.tsv");
        meta<<std::setprecision(17)<<"diagnostic_status\t"<<status<<"\ndiagnostic_complete\t"<<(done&&!invalid)
            <<"\ndiagnostic_window_ps\t"<<window_ps<<"\ndiagnostic_windows_requested\t"<<window_count
            <<"\ndiagnostic_range_rtol\t"<<range_rtol<<"\ndiagnostic_range_atol\t"<<range_atol
            <<"\ndiagnostic_mean_rtol\t"<<mean_rtol<<"\ndiagnostic_accepted_at_ps\t"<<accepted_at
            <<"\ndiagnostic_target_span_ps\t"<<target_span<<"\ndiagnostic_target_cycles\t"<<target_cycles
            <<"\ndiagnostic_detected_period_ps\t"<<detected_period<<"\ndiagnostic_max_duration_ps\t"<<3.*window_count*target_span<<"\n";
        std::ofstream out(folder+"/diagnostic_windows.tsv");
        out<<std::setprecision(17)<<"start_ps\tend_ps\tmean_intensity\tminimum_intensity\tmaximum_intensity\tpeak_to_peak\trelative_interval_range\tcycles\n";
        for(const auto& w:windows) {
            out<<w.start<<'\t'<<w.end<<'\t'<<w.mean<<'\t'<<w.low<<'\t'<<w.high<<'\t'<<w.high-w.low<<'\t';
            if(w.mean>0.) out<<(w.high-w.low)/w.mean; else out<<"null";
            out<<'\t'<<w.cycles<<'\n';
        }
    }
};
#endif
