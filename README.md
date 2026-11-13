# Nano-Shell Simulation Project

This project provides tools for simulating the optical and plasmonic properties of nanoshell structures, particularly in the quasi-static regime, both below and above the emission threshold. The core theory behind this project has been published in the open-access paper:

**Reference**: ["Emission and lasing properties of a single nanoshell spaser nanoparticle"](https://www.degruyter.com/document/doi/10.1515/nanoph-2024-0491/html), *Nanophotonics* (2024).

---

## Overview

This repository contains:

- Source codes for simulating nanoshells with gain-enhanced cores.
- Headers defining material models and key functions.
- Tools for calculating steady-state polarizabilities and time-dependent dynamics.
- Example input files and testing scripts.

The code uses the GNU Scientific Library (GSL) and Armadillo for numerical computations.

---

## Compilation

### Requirements

To compile the project, ensure the following dependencies are installed:

- A C++ compiler (e.g., g++)
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- [Armadillo C++ library](http://arma.sourceforge.net/)

### Steps

1. Run the `configure` script to check dependencies:

   ```bash
   ./configure
   ```

2. Compile the project using the Makefile:

   ```bash
   make
   ```

3. To clean compiled binaries:

   ```bash
   make clean
   ```

---

## Input Files

### 1. `data/input/nanosphere_eV.dat`

This file describes the nanoshell structure parameters:

| Parameter                  | Description                                                                                     | Example Value |
|----------------------------|-------------------------------------------------------------------------------------------------|---------------|
| External radius (nm)       | Radius of the nanoshell, \(a\) in the paper                                                     | 10            |
| Emission width (eV)        | \(\hbar \Delta\), amplitude of the emission permittivity of the gain medium                     | 0.15          |
| Emission frequency (eV)    | \(\hbar \omega_g\), emission frequency of the gain medium                                       | 2.8122        |
| Gain level                 | \(G\), gain level (it also works with \(-G\) for retro compatibility)                           | 0.033         |
| Spectrum bounds (eV)       | \(\omega_{min}\), \(\omega_{max}\), bounds of the spectral range of interest                    | 2.0, 4.5      |
| Metal type                 | Material of the shell (e.g., `silver`)                                                          | silver        |
| Metal model                | Model for the metal (`drude`, `spline`)                                                         | drude         |
| Gain medium model          | Model for the active core medium (`lorentz`, `flat`)                                            | lorentz       |
| Host material              | Surrounding medium (e.g., `water`)                                                              | water         |
| Excitation field amplitude | \(E_0\), excitation field (arbitrary units); compare with saturation field amplitude \(E_{sat}\)| 1.e-30        |
| Radius ratio               | \(\rho\), ratio of inner to outer radius                                                        | 0.6           |
| Core material              | Material of the nanoshell core (e.g., `silica`)                                                 | silica        |

### 2. `data/input/time.dat`

This file specifies time-dependent simulation parameters:

| Parameter         | Description                                             | Example Value |
|-------------------|---------------------------------------------------------|---------------|
| Total time (ps)   | Duration of the simulation.                             | 100           |
| Pump-on time (ps) | Time at which the pump is switched on (in ns).          | 10            |

---

## Code Usage

### Including `cup.H`

To include and use the functions defined in `cup.H`, you must include the following headers in your code in the specified order:

```cpp
#include "../src/headers/math33.hpp"
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#include "../src/headers/cup.hpp"
```

### Material Initialization

Before using subroutines in `cup.H`, you must:

1. Define a `nanosphere` instance (e.g., `nanosphere ns`).
2. Initialize it using:

   ```cpp
   ns.init();
   ```

3. Set the metal and active medium using:

   ```cpp
   ns.set_metal("<metal type>", "<model>", 1);
   ns.set_active("<gain medium model>");
   ```

---

## Test Programs

The `tests` directory includes sample programs demonstrating the functionality of the Nano-Shell Simulation Project. These programs can be compiled and executed to verify specific features of the code. Below is an overview of the included test programs:

- **time_behavior.cxx**: Simulates the time-dependent behavior of the complex dipole moment for a given frequency.
- **threshold.cxx**: Calculates the threshold gain \( G_{\text{th}} \) and the threshold frequency \( \omega_{\text{th}} \), as defined in the associated publication.
- **steady_state.cxx**: Computes the steady-state spectrum of the polarizability.

### Compilation and Execution

To compile and run the test programs, navigate to the `tests` directory and use the following commands:

```bash
g++ -Wall -I../src/headers -L../src/lib time_behavior.cxx -o time_behavior -lgsl -lgslcblas -lm -larmadillo
./time_behavior

g++ -Wall -I../src/headers -L../src/lib threshold.cxx -o threshold -lgsl -lgslcblas -lm -larmadillo
./threshold

g++ -Wall -I../src/headers -L../src/lib steady_state.cxx -o steady_state -lgsl -lgslcblas -lm -larmadillo
./steady_state
```

---

## License

This project is released under the GNU GENERAL PUBLIC LICENSE. See `LICENSE` for details.

---

## Contributing

Contributions are welcome! If you find any bugs or have feature requests, feel free to open an issue or submit a pull request.

---

## Contact

For questions or collaboration opportunities, please reach out to Alessandro Veltri at <alessandro.veltri@gmail.com>.

## Validate automatic convergence against Supplementary Fig. 1

`scripts/validate_convergence.bash` builds an optimized
`bin/emission_convergence` driver against the pinned `nano_geo_matrix` submodule.
Run the wrapper from any directory. It requires the saved Fig. 1 input snapshots
and traces in `data/output/[nanoph-2024-0491]_rangeE0/`, plus gnuplot and Python 3.
It preserves those files and all repo inputs.

The workflow runs the four raw seeds at the saved fixed duration, then runs them
with `T=0` at default convergence settings and with all relative tolerances halved.
It uses only the shared CUP detector, with 0.1 ps automatic output spacing; the
analytical solution is not needed for this numerical experiment. Results live in
`data/output/convergence_validation/`. Existing validation output
causes the wrapper to stop rather than overwrite a previous experiment.

The summary script checks fixed-duration traces byte-for-byte against the saved
reproduction and reads the solver's convergence metadata into `summary.tsv`.
Peak-to-peak statistics come from the detector at full integration resolution,
normalized by the saved paper `Isat`. The PNG in `img/output/` shows accepted tails
including confirmation, with default and tighter settings overlaid. That preview
uses sampled output and the solver's existing rounded field columns; very weak
ripples may show serialization noise. It is not used to calculate the statistics.

Flat classifications are tolerance-based: a small periodic ripple or residual
transient may be below flatness resolution. Their reported collective variation
is not necessarily a converged periodic amplitude. Compare tolerance settings
before interpreting weak-seed peak-to-peak values physically.

## Residual modulation scan (no figures)

`scripts/scan_residual_modulation.py` compiles the convergence driver with `-O2`
inside its output directory. Supply `--input INPUT --time TIME --frequencies
FREQUENCIES.json --output NEW_DIRECTORY`; paths are relative to your current
working directory. The frequency JSON must contain exactly the keys `1.05`,
`1.10`, `1.15`, each with an explicit list of positive excitation energies in eV.
An empty list skips a gain. The explicit token `"threshold"` selects the computed
threshold energy for a smoke check. No frequency ranges are generated. Select
actual off-threshold values yourself; `2.811` eV is the existing validation
excitation, not a proposed scan boundary.

Each point runs at both tolerance factors 1 and 0.5, with `E0/Esat=1e-9`.
`--cadence-ps` defaults to 0.1. `--timeout-seconds` defaults to 120 per run;
this is a wall-clock kill limit, **not** a physical-time cap. The solver has no
configurable automatic physical ceiling: its unchanged limit is 100000 ps.
Timeouts and failures leave the later observation unresolved and are retried
on resume. Any original accepted convergence evidence already saved is retained. Completed points (including maximum-duration results) are
reused only when the manifest matches. Use a new directory to change settings.
Do not run concurrent scans into the same directory.

Normalization is now verified against the gain-switching manuscript and the
working units in [Recalde et al. (2023)](https://doi.org/10.3390/ma16051911).
The corrected scan uses `Esat_num=sqrt(Dome/(10*abs(G)*Ome_p))` and
`raw_E0=1e-9*Esat_num`. The earlier scan omitted `Ome_p`. Historical reproduction
utilities and raw-input calls retain their behavior; old scan results require
regeneration. [The derivation and diagnostic definitions](docs/residual_modulation.md)
explain the density-dependent numerical unit and record all SI conversion
factors. The manuscript assumption is mu=10 D.

The convergence-plus-confirmation range is now explicitly an **interval
variation**, stored as `accepted_interval_peak_to_peak` and
`accepted_relative_interval_range`. It is not automatically a persistent
modulation amplitude. Means and ranges use the shared detector's full-resolution
`norm(dip)` statistics. Flat classification remains tolerance dependent.

The scan continues the same physical trajectory after acceptance by default,
observing three later windows of at least 40 ps each. Configure
`--diagnostic-window-ps`, `--diagnostic-windows`, `--range-rtol`, `--range-atol`,
and `--mean-rtol`; use `--diagnostic-window-ps 0` to disable observation. Durations above the
solver's 100000 ps ceiling are rejected.
Periodic observations use complete cycles identified by CUP's existing extrema
tracker and span at least five detected periods per window. Flat observations
use the configured duration and full-resolution trapezoidal means. The original
convergence decision/time remains separate from the later terminal time.

`summary.json` and `summary.tsv` contain normalization, convergence and diagnostic
settings/results. `late_windows.tsv` records each late window's endpoints,
mean, min/max, range, relative range and cycle count. The diagnostic reports
stable interval variation, decaying variation, mean drift, or unresolved evidence
using the documented comparison tolerances. Results from each convergence
tolerance remain separate. Zero-mean or incomplete diagnostics stay unresolved;
timeouts never imply zero modulation. No persistent-amplitude extrapolation is
made.

The detector may accept a finite plateau before later evolution or a beat whose
variation falls below its resolution across the 40 ps default flat evidence and
confirmation interval. Confirmation reduces but cannot exclude that ambiguity;
flat acceptance is not proof of absence of a slow beat. Tightening tolerances
is a sensitivity check, not such a proof. No new convergence rule is added.

Snapshots include input/time files, frequency lists, driver/scan sources,
repository and submodule commits and patches, source hashes, build log and each
run's command/log. Shared inputs, paper outputs, and earlier validation results
are untouched. The driver retains its old five-argument interface (2.811 eV);
optional trailing arguments are excitation energy (or `threshold`), gain ratio,
normalized field ratio, output cadence, diagnostic window duration (0 disables),
window count, range relative tolerance, range absolute tolerance, and mean
relative tolerance. `prepare` mode writes normalization
metadata without integration. Existing fixed/automatic callers remain valid.

For an auditable scan, `--manuscript PATH` also snapshots the manuscript source.
