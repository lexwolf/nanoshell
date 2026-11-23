"""Plot and measure the dominant stationary pulled envelope."""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = Path(__file__).resolve().parent
COLORS = ("#264d78", "#ad482a")
NAMES = ("Silver Drude", "Artificial Drude–Lorentz")


def load(rho, strength, count):
    stem = HERE / f"rho{rho}_S{strength}_{count}"
    data = np.loadtxt(stem.with_suffix(stem.suffix + ".csv"), delimiter=",", skiprows=1)
    params = np.loadtxt(stem.with_suffix(stem.suffix + ".txt"))
    positive = np.flatnonzero(data[:, 2] > 0)
    assert len(positive) > 10 and np.all(np.diff(positive) == 1)
    assert positive[0] > 0 and positive[-1] < len(data)-1
    return data, data[positive], params


def envelope(points):
    # Interpolate each monotonic piece independently. Take the upper curve
    # where distinct reference frequencies map to the same emitted energy.
    emitted, intensity = points[:, 1:3].T
    turns = 1 + np.flatnonzero(np.diff(np.sign(np.diff(emitted))))
    grid = np.unique(np.r_[np.linspace(emitted.min(), emitted.max(), 8001), emitted])
    upper = np.zeros_like(grid)
    for start, stop in zip(np.r_[0, turns], np.r_[turns, len(emitted)-1]):
        indices = np.argsort(emitted[start:stop+1])
        x = emitted[start:stop+1][indices]
        y = intensity[start:stop+1][indices]
        upper = np.maximum(upper, np.interp(grid, x, y, left=0, right=0))
    return grid, upper


def asymmetry(x, y):
    peak = np.argmax(y)
    half = y[peak]/2
    left = np.flatnonzero(y[:peak] < half)
    right = np.flatnonzero(y[peak+1:] < half)
    # A band edge acts as a vertical closure if the peak extends to it.
    l = x[0]
    r = x[-1]
    if left.size:
        j = left[-1]
        l = np.interp(half, y[j:j+2], x[j:j+2])
    if right.size:
        j = peak + right[0]
        r = np.interp(half, y[j:j+2][::-1], x[j:j+2][::-1])
    return (r+l-2*x[peak])/(r-l), l, r


def report():
    baseline = np.loadtxt(HERE.parent/"gain_switching_corrected_20261114/tables/pulling_corrected.tsv")
    silver, _, _ = load("0.6", "0", 9400)
    interior = baseline[:, 2] > .01
    reference = baseline[interior]
    errors = np.abs(np.interp(reference[:, 0], silver[:, 0], silver[:, 2])-reference[:, 2])
    emitted_errors = np.abs(np.interp(reference[:, 0], silver[:, 0], silver[:, 1])-reference[:, 1])
    print("Silver benchmark max absolute errors (J, emitted eV):", errors.max(), emitted_errors.max())
    assert errors.max() < 1e-4 and emitted_errors.max() < 1e-4

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    rows = []
    with PdfPages(HERE/"Artificial_metal_pulling.pdf") as pdf:
        for rho in ("0.6", "0.8"):
            fig = plt.figure(figsize=(8.27, 11.69))
            top, bottom = fig.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1.3]})
            fig.subplots_adjust(left=.13, right=.93, top=.87, bottom=.24, hspace=.53)
            fig.text(.13, .943, "Does the emission fin reverse?", fontsize=24, weight="bold")
            fig.text(.13, .906, rf"$\rho={rho}$   |   $G/G_{{\mathrm{{th}}}}=1.20$   |   Gain center aligned to each threshold",
                     fontsize=12, color="#445166")
            top.set_title("Metal permittivity", loc="left", weight="bold")
            bottom.set_title("Frequency pulling", loc="left", weight="bold")
            zoom = bottom.inset_axes([.56, .55, .42, .41])
            real = top.inset_axes([.37, .47, .33, .36])
            limits = []
            summary = []
            for name, strength, color in zip(NAMES, ("0", "1.5"), COLORS):
                data, active, (eth, gth, isat) = load(rho, strength, 9400)
                grid, upper = envelope(active)
                coarse_a = asymmetry(*envelope(load(rho, strength, 4700)[1]))[0]
                a, left, right = asymmetry(grid, upper)
                assert abs(a-coarse_a) < .015
                slope = np.gradient(data[:, 4], data[:, 0])
                band = (data[:, 0] >= active[0, 0]) & (data[:, 0] <= active[-1, 0])
                if strength != "0":
                    assert slope[band].min() > 0
                    top.axvspan(active[0, 0], active[-1, 0], color=color, alpha=.09)
                top.plot(data[:, 0], data[:, 4], lw=2, color=color, label=name)
                top.axvline(eth, color=color, alpha=.5, lw=.9, ls=":")
                real.plot(data[:, 0], data[:, 3], color=color, lw=1.3)
                visible = (data[:, 0] >= active[0, 0]) & (data[:, 0] <= active[-1, 0])
                bottom.plot(data[visible, 0], data[visible, 2], color=color, ls="--", lw=1.4)
                for axis in (bottom, zoom):
                    axis.fill_between(grid, upper, color=color, alpha=.16)
                    axis.plot(grid, upper, color=color, lw=1.8)
                zoom.hlines(upper.max()/2, left, right, colors=color, linestyles=":", lw=1.2)
                limits.extend([grid[0], grid[-1]])
                summary.append((name, eth, gth, a))
                rows.append([rho, name, eth, gth, isat, active[0, 0], active[-1, 0],
                             grid[0], grid[-1], slope[band].min(), slope[band].max(), a, coarse_a])
            lo, hi = ((1.95, 3.35) if rho == "0.6" else (1.65, 2.73))
            for axis in (top, bottom, real):
                axis.set_xlim(lo, hi)
            top.set_ylim(0, max(np.loadtxt(HERE/f"rho{rho}_S0_9400.csv", delimiter=",", skiprows=1)[
                (data[:, 0] >= lo) & (data[:, 0] <= hi), 4])*1.15)
            top.set(xlabel="Energy (eV)", ylabel=r"$\mathrm{Im}\,\varepsilon_m$")
            real.set_title(r"$\mathrm{Re}\,\varepsilon_m$", fontsize=9)
            real.tick_params(labelsize=8)
            top.legend(loc="upper right", frameon=False)
            top.text(.025, .06, "Shaded band: emitting references for the artificial metal",
                     transform=top.transAxes, fontsize=8, color="#445166")
            bottom.set(xlabel="Reference or emitted energy (eV)", ylabel=r"$I_{\rm em}/I_{\rm sat}$", ylim=(0, .5))
            zoom.set(xlim=(min(limits)-.006, max(limits)+.006), ylim=(0, .47))
            zoom.set_title("Pulled fins at emitted energy", fontsize=9)
            zoom.tick_params(labelsize=8)
            bottom.text(.025, .96, "Dashed: reference energy\nFilled: emitted energy",
                        va="top", transform=bottom.transAxes, fontsize=9)
            for axis in (top, bottom, zoom):
                axis.grid(alpha=.15)
                axis.set_axisbelow(True)
            verdict = "reversal" if summary[0][3]*summary[1][3] < 0 else "no reversal"
            fig.text(.13, .183, f"Half-width asymmetry: {verdict}", weight="bold", fontsize=13)
            for index, (name, eth, gth, a) in enumerate(summary):
                fig.text(.13, .153-index*.023,
                         f"{name}    threshold {eth:.5f} eV    $G_{{th}}$ = {gth:.5f}    $A$ = {a:+.3f}",
                         color=COLORS[index], fontsize=10)
            fig.text(.13, .080,
                     "A = (right − left)/(right + left), using half-widths of the upper pulled curve.\n"
                     "Both curves use the reduced stationary pulling model. The Lorentz response is\n"
                     "evaluated at the reference energy. The comparison tests this model’s fin orientation.",
                     fontsize=9, linespacing=1.5, color="#445166")
            fig.text(.13, .043,
                     "Artificial metal: oscillator at 3.8 eV, width 1.0 eV, strength 1.5 eV².\n"
                     "Drude damping is 1.5% of silver’s. The real permittivity is matched at 2.8 eV.",
                     fontsize=8.5, linespacing=1.5, color="#445166")
            pdf.savefig(fig)
            plt.close(fig)
    with (HERE/"summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["rho", "metal", "threshold_eV", "Gth", "Isat", "reference_min_eV", "reference_max_eV",
                         "emitted_min_eV", "emitted_max_eV", "slope_min_per_eV", "slope_max_per_eV", "A", "A_coarse"])
        writer.writerows(rows)


if __name__ == "__main__":
    report()
