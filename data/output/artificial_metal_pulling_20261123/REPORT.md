# Artificial-metal stationary pulling test

At both radius ratios, the artificial metal's imaginary permittivity rises
across its emitting reference-energy band. The pulled fin becomes more
symmetric and retains its high-energy tail in the reduced stationary model.

| Radius ratio | Artificial emitting reference band (eV) | Silver asymmetry | Artificial asymmetry |
| --- | --- | --- | --- |
| 0.6 | 2.12925–2.91800 | +1.000 | +0.1013 |
| 0.8 | 1.71350–2.37650 | +1.000 | +0.1420 |

Here the asymmetry is (right half-width − left half-width) divided by their
sum, measured on the upper envelope of the pulled intensity. Halving the
reference grid spacing changed the asymmetry by less than 0.0015. The
rho=0.6 silver benchmark matched the archived normalized intensity to
4.16e-8 and the emitted energy to 3.43e-9 eV at the checked interior points.

The artificial response adds a Lorentz oscillator at 3.8 eV, with width
1.0 eV and strength 1.5 eV², to silver Drude with damping reduced to 1.5%
of its original value. A real offset cancels the oscillator's real
contribution at 2.8 eV. Each run aligns its gain center to its bright-mode
threshold and sets G/Gth to 1.20. The artificial material also changes
absolute loss and threshold gain.

The Lorentz contribution is sampled at the reference energy. This experiment
uses the reduced envelope equation and introduces no dynamical Lorentz
polarization state. The result describes this material and this stationary
construction. See the
[extended report](https://chatgpt.com/library?file_id=libfile_f6d0f39b7b788191ba4d4e3434b29cc7&account_id=personal)
for the figures, equations and implementation excerpts.

## Reproduce

Run from the nanoshell repository root with its pinned submodule initialized.
The build copies headers into an experiment-local overlay. It leaves the
tracked library and solver sources untouched.

```bash
cd data/output/artificial_metal_pulling_20261123
bash build.sh
for rho in 0.6 0.8; do
    for strength in 0 1.5; do
        for intervals in 4700 9400; do
            stem="rho${rho}_S${strength}_${intervals}"
            ./driver "$rho" "$strength" "$intervals" > "$stem.csv" 2> "$stem.txt" || exit 1
        done
    done
done
python3 report.py
```

The original run used nanoshell revision `6f38bb7`, shared library revision
`9689dad9c276cf107a8e26d3b93baeeade7e8df7` and GCC 15.2.0.
`report.py` reads the archived silver benchmark from the adjacent corrected
gain-switching output directory.
