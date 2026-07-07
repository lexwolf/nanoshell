/*
 * This file is part of the Nano-Shell Simulation Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <armadillo>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#define CUP_BACKEND_QUASI_STATIC
#include <cup/cup.hpp>

using namespace std;

struct FrohlichRoot {
  double omega;
  double Gth;
  double branch_score;
};

string normalize_branch_mode(const char* mode_arg) {
  string mode = mode_arg ? mode_arg : "bright";
  if (mode == "dark") mode = "secondary";
  if (mode == "primary") mode = "bright";

  if (mode != "bright" && mode != "secondary" &&
      mode != "lowest" && mode != "highest") {
    cerr << "Error: unknown branch mode '" << mode << "'." << endl;
    cerr << "Valid modes are: bright, secondary, dark, lowest, highest." << endl;
    exit(2);
  }

  return mode;
}

double frohlich_residual(nanosphere& ns, double omega, double eps_b, double eps_s, double rho) {
  const double rho3 = rho*rho*rho;
  complex<double> eps1 = ns.metal(omega);
  complex<double> fval = F(eps1, eps_s, rho3);
  return eqn9(omega, omega, fval, eps_b, ns.Dome);
}

FrohlichRoot refine_frohlich_root(nanosphere& ns, double left, double right,
                                  double eps_b, double eps_s, double rho) {
  double fl = frohlich_residual(ns, left, eps_b, eps_s, rho);
  double mid = left;

  for (int ii=0; ii<120 && fabs(right-left)>1.e-10; ii++) {
    mid = 0.5*(left+right);
    double fm = frohlich_residual(ns, mid, eps_b, eps_s, rho);
    if (fl*fm <= 0.) {
      right = mid;
    } else {
      left = mid;
      fl = fm;
    }
  }

  double omega = 0.5*(left+right);
  const double rho3 = rho*rho*rho;
  complex<double> eps1 = ns.metal(omega);
  complex<double> fval = F(eps1, eps_s, rho3);

  FrohlichRoot root;
  root.omega = omega;
  root.Gth = -imag(fval);
  root.branch_score = 0.;
  return root;
}

vector<FrohlichRoot> find_frohlich_roots(nanosphere& ns, double omi, double oma,
                                         double eps_b, double eps_s, double rho) {
  vector<FrohlichRoot> roots;
  const int samples = 20000;
  const double dome = (oma-omi)/samples;
  double x_prev = omi;
  double f_prev = frohlich_residual(ns, x_prev, eps_b, eps_s, rho);

  for (int ii=1; ii<=samples; ii++) {
    double x = omi + ii*dome;
    double f = frohlich_residual(ns, x, eps_b, eps_s, rho);

    if (f_prev == 0. || f_prev*f <= 0.) {
      FrohlichRoot root = refine_frohlich_root(ns, x_prev, x, eps_b, eps_s, rho);
      bool duplicate = false;
      for (const auto& existing : roots) {
        if (fabs(existing.omega-root.omega) < 1.e-5) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate && root.Gth > 0. && root.Gth < 1. && isfinite(root.Gth)) {
        roots.push_back(root);
      }
    }

    x_prev = x;
    f_prev = f;
  }

  return roots;
}

double score_stationary_response(nanosphere ns, const FrohlichRoot& root,
                                 char* mdl, char* mtl, char* hst,
                                 double omi, double oma, char* sol, double rho) {
  ns.ome_g = root.omega;
  ns.G = 1.2*root.Gth;
  vector<pair<double,complex<double>>> response =
    ns.steady_state(mdl, mtl, hst, omi, oma, 10000, sol, rho);

  double peak = 0.;
  for (const auto& point : response) {
    double amplitude = abs(point.second);
    if (amplitude > peak) peak = amplitude;
  }
  return peak;
}

FrohlichRoot select_branch(vector<FrohlichRoot>& roots, const string& branch_mode, double rho) {
  if (roots.empty()) {
    cerr << "Error: no physical Frohlich roots found for rho=" << rho << "." << endl;
    exit(3);
  }

  if (branch_mode == "lowest") {
    return *min_element(
      roots.begin(), roots.end(),
      [](const FrohlichRoot& lhs, const FrohlichRoot& rhs) {
        return lhs.omega < rhs.omega;
      }
    );
  }

  if (branch_mode == "highest") {
    return *max_element(
      roots.begin(), roots.end(),
      [](const FrohlichRoot& lhs, const FrohlichRoot& rhs) {
        return lhs.omega < rhs.omega;
      }
    );
  }

  if (branch_mode == "secondary") {
    return *max_element(
      roots.begin(), roots.end(),
      [](const FrohlichRoot& lhs, const FrohlichRoot& rhs) {
        return lhs.omega < rhs.omega;
      }
    );
  }

  sort(
    roots.begin(), roots.end(),
    [](const FrohlichRoot& lhs, const FrohlichRoot& rhs) {
      return lhs.branch_score > rhs.branch_score;
    }
  );

  return roots[0];
}

/** Compila con: 
Example compilation:

NGM_ROOT="$(realpath ../extern/nano_geo_matrix)" && g++ -I"$NGM_ROOT/include" -I"$NGM_ROOT/modules" -I"$NGM_ROOT/modules/cup" -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
**/

int main(int argc, char ** argv){
    if (argc==1){
        cout<<endl<<" Usage: "<<argv[0]<<" <rho> [bright|secondary|dark|lowest|highest]"<<endl<<endl;
        exit(1);
        }
  double rho = atof(argv[1]);
  string branch_mode = normalize_branch_mode(argc >= 3 ? argv[2] : "bright");
  double   omemi, omema, eps_b, E0, eps3, rho_file;
  double *result;
  FrohlichRoot selected;
  
  char mtl[16], mdl[16], sol[16], hst[16], active[16];
  nanosphere ns;
  ifstream nano("../data/input/nanosphere_eV.dat");
  if (!nano) {
    cerr << "Error: Cannot open input file" << endl;
    return 1;
  }

  nano>>ns.a>>ns.Dome>>ns.ome_g>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0>>rho_file>>hst;
  (void)rho_file;

  ns.init();
  eps3=ns.set_host(sol);
  eps_b=ns.set_host(hst);
  ns.set_metal(mtl,mdl,1);
  ns.set_active(active);

  double search_omi = 1.5;
  double search_oma = 4.5;
  if (string(mtl) == "silver" &&
      (branch_mode == "secondary" || branch_mode == "highest")) {
    search_omi = 3.0;
    search_oma = 6.5;
  }

  vector<FrohlichRoot> roots = find_frohlich_roots(ns, search_omi, search_oma, eps_b, eps3, rho);

  if (roots.empty()) {
    result=ns.frohlich(search_omi, search_oma, eps_b, eps3, rho);
    selected.omega = result[0];
    selected.Gth = result[1];
    selected.branch_score = 0.;
  } else {
    for (auto& root : roots) {
      root.branch_score =
        score_stationary_response(ns, root, mdl, mtl, hst, omemi, omema, sol, rho);
    }

    selected = select_branch(roots, branch_mode, rho);
  }

  cout.precision(7);        //set the precision
  cout.setf(ios::fixed);
  cout<<selected.omega<<" "<<selected.Gth;
  return 0;
  }
