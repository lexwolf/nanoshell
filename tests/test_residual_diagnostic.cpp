#include <cup/convergence.hpp>
#include "../src/residual_diagnostic.hpp"
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

int main() {
    cup_convergence::Result result;
    result.classification="sustained_flat";
    {
        ResidualDiagnostic d({}); d.window_ps=1.;
        for(int i=0;i<90 && !d.done;++i) d.update(i*.037,2.,result);
        assert(d.done && d.windows.size()==3);
        for(auto w:d.windows) { assert(std::abs(w.end-w.start-1.)<1e-12); assert(w.mean==2. && w.low==2. && w.high==2.); }
    }
    {
        // Linear ramp with nonaligned boundaries has exact trapezoidal means.
        ResidualDiagnostic d({}); d.window_ps=1.;
        for(int i=0;i<90 && !d.done;++i) d.update(i*.037,2.+i*.037,result);
        assert(d.done);
        for(auto w:d.windows) assert(std::abs(w.mean-(2.+(w.start+w.end)*.5))<1e-12);
    }
    {
        result.classification="sustained_periodic"; result.mean_period_ps=2.;
        ResidualDiagnostic d({}); d.window_ps=4.;
        for(int i=0;i<40000 && !d.done;++i) d.update(i*.001,2.+.2*std::cos(i*.001*std::acos(-1.)),result);
        assert(d.done && d.windows.size()==3);
        for(auto w:d.windows) {
            assert(w.cycles>=5 && w.end-w.start>=10.);
            assert(std::abs(w.mean-2.)<1e-10);
            assert(std::abs((w.high-w.low)-.4)<1e-10);
            assert(std::abs(w.start/2.-std::round(w.start/2.))<1e-10);
            assert(std::abs(w.end/2.-std::round(w.end/2.))<1e-10);
        }
    }
    {
        result.classification="sustained_flat";
        ResidualDiagnostic d({}); d.window_ps=1.;
        assert(d.update(0.,std::numeric_limits<double>::quiet_NaN(),result));
        assert(d.invalid && !d.done);
    }
    {
        // Zero mean is unresolved, rather than a claimed zero persistent amplitude.
        ResidualDiagnostic d({});d.window_ps=1.;
        for(int i=0;i<400 && !d.done;++i) d.update(i*.01,0.,result);
        auto out=std::filesystem::temp_directory_path()/("residual-diagnostic-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(out);d.write(out.string());
        assert(d.status=="unresolved");
        std::filesystem::remove_all(out);
    }
    {
        auto check=[&](auto signal, const std::string& expected) {
            ResidualDiagnostic d({}); d.window_ps=10.;
            for(int i=0;i<=30001 && !d.done;++i) d.update(i*.001,signal(i*.001),result);
            auto out=std::filesystem::temp_directory_path()/("residual-classification-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
            std::filesystem::create_directories(out);d.write(out.string());
            if(!d.done || d.status!=expected) std::cerr<<"expected "<<expected<<", got "<<d.status<<"\n";
            assert(d.done && d.status==expected);
            std::filesystem::remove_all(out);
        };
        check([](double){return 2.;},"stable_interval_variation");
        check([](double t){return 2.+.01*t;},"decaying_variation_and_mean_drift");
        check([](double t){return 2.+.1*std::exp(-t/10.)*std::sin(2.*std::acos(-1.)*t);},"decaying_variation_and_mean_drift");
        check([](double t){return 2.+(.1+.003*t)*std::sin(2.*std::acos(-1.)*t);},"unresolved");
    }
    {
        result.classification="sustained_periodic";result.mean_period_ps=1.;
        ResidualDiagnostic d({});d.window_ps=1.;
        for(int i=0;i<5001 && !d.invalid;++i) d.update(i*.01,2.,result);
        assert(d.invalid && !d.done && d.windows.empty());
    }
    std::cout<<"Diagnostic boundaries, means, shared complete cycles, nonfinite and zero-mean tests passed\n";
}
