"""Render the linear-scale variant from existing gain_dependence_plot.dat metrics.

Run gain_dependence.py first only if the underlying simulation data have changed.
This entry point leaves the original figure, metrics, and LaTeX folder unchanged.
"""
from pathlib import Path
import subprocess

if __name__ == '__main__':
    base = Path(__file__).resolve().parents[1]
    if not (base / 'Pub_Ready/gain_dependence_plot.dat').is_file():
        raise SystemExit('Missing metrics: run gain_dependence.py first.')
    subprocess.run(['gnuplot', str(base / 'gnuplot/gain_dependence_linear.gp')], check=True)
