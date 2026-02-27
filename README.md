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
#include "../src/headers/nanoshell.hpp"
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

For questions or collaboration opportunities, please reach out to Alessandro Veltri at alessandro.veltri@gmail.com.
# Tiny sync fix
