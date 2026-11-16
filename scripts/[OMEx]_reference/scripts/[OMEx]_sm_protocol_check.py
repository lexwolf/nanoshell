"""Audit existing records used in the Supplement 1 numerical checks.

Run with python3 Gain_Switching/scripts/sm_protocol_check.py from any folder.
No simulation is run and no source data or figure is changed. The only output
is Pub_Ready/sm_protocol_check.json, with source hashes and extraction rules.
"""
from pathlib import Path
import hashlib
import json
import re

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "Gain_Switching/Pub_Ready/sm_protocol_check.json"
NOTEBOOK = ROOT / "NEW_stuff/TimeBehavior.ipynb"
REPRESENTATIVE = ROOT / "NEW_stuff/Data/Representative"
PROBE = ROOT / "Gain_Switching/Data/ProbeSweep/middle"
ESAT = 0.317931
END_PS = 200.0
LATE_START_PS = 150.0


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def header(path):
    values = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            match = re.match(
                r"#\s*(G|ome|ome_em|Esat)\s*=\s*([+\-\d.eE]+)", line
            )
            if match:
                key, value = match.groups()
                if key in values:
                    raise ValueError(f"Repeated header field {key}: {path}")
                values[key] = float(value)
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError(f"Nonfinite header: {path}")
    return values


def first_crossing(t, intensity, level):
    indices = np.flatnonzero(
        (intensity[:-1] < level) & (intensity[1:] >= level)
    )
    if not len(indices):
        raise ValueError("No ascending intensity crossing")
    i = int(indices[0])
    if t[i + 1] <= t[i]:
        raise ValueError("Rounded timestamps do not resolve this crossing")
    crossing = t[i] + (level - intensity[i]) * (t[i + 1] - t[i]) / (
        intensity[i + 1] - intensity[i]
    )
    return float(crossing), [float(t[i]), float(t[i + 1])]


def representative_checks():
    notebook = json.loads(NOTEBOOK.read_text())
    records, arrays = [], []
    for cell_index, name, step in [(0, "1.1Gth.dat", 1.0), (1, "1.1Gthadj.dat", 0.5)]:
        cell = "".join(notebook["cells"][cell_index]["source"])
        label = re.search(r"^# dt\s*=\s*([\d.]+)", cell, re.MULTILINE)
        commented_path = f'# file_path = "Data/Representative/{name}"'
        if label is None or float(label.group(1)) != step or commented_path not in cell:
            raise ValueError(f"Notebook step-label provenance changed for {name}")
        path = REPRESENTATIVE / name
        values = header(path)
        if not {"G", "ome", "ome_em"} <= values.keys():
            raise ValueError(f"Missing representative header fields: {path}")
        data = np.loadtxt(path)
        if data.ndim != 2 or data.shape[1] != 9 or not np.isfinite(data).all():
            raise ValueError(f"Unexpected or nonfinite representative data: {path}")
        t = data[:, 0]
        if np.any(np.diff(t) < 0) or t[0] < 0 or t[0] >= 0.001 or t[-1] != END_PS:
            raise ValueError(f"Unexpected time coverage/order: {path}")
        records.append(dict(
            source=relative(path), sha256=digest(path), header=values,
            notebook_cell_index=cell_index, notebook_step_label=step,
            saved_rows=int(len(data)), first_saved_time_ps=float(t[0]),
            last_saved_time_ps=float(t[-1]),
        ))
        arrays.append(data)
    if records[0]["header"] != records[1]["header"]:
        raise ValueError("Representative records have different parameter headers")
    if records[0]["sha256"] == records[1]["sha256"]:
        raise ValueError("Representative records are duplicate files")
    if records[0]["header"]["G"] != 0.148397 or records[0]["header"]["ome"] != 2.8122:
        raise ValueError("Representative parameters no longer match the reported case")

    for record, data in zip(records, arrays):
        selected = data[(data[:, 0] >= 0) & (data[:, 0] <= END_PS)]
        t = selected[:, 0]
        intensity = (selected[:, 3] ** 2 + selected[:, 4] ** 2) / ESAT**2
        peak = float(intensity.max())
        onset, bracket = first_crossing(t, intensity, 0.81 * peak)
        late = t >= LATE_START_PS
        if late.sum() < 2 or t[late][0] > LATE_START_PS + 0.001:
            raise ValueError("Insufficient late-window coverage")
        # Dip/alpha is only a consistency diagnostic for E0, not run metadata.
        alpha = complex(selected[0, 1], selected[0, 2])
        inferred_seed = complex(selected[0, 3], selected[0, 4]) / alpha
        record.update(
            Esat=ESAT, t_on_ps=onset, crossing_bracket_ps=bracket,
            peak_I_over_Isat=peak, threshold_I_over_Isat=0.81 * peak,
            late_mean_I_over_Isat=float(intensity[late].mean()),
            late_saved_samples=int(late.sum()),
            late_first_saved_time_ps=float(t[late][0]),
            last_100000_sample_mean_I_over_Isat=float(intensity[-100000:].mean()),
            last_100000_sample_start_ps=float(t[-100000]),
            final_N=float(selected[-1, 6]),
            first_saved_dip_over_alpha_real=float(inferred_seed.real),
            first_saved_dip_over_alpha_imag=float(inferred_seed.imag),
        )
    difference = abs(records[1]["t_on_ps"] - records[0]["t_on_ps"])
    return dict(
        notebook_source=relative(NOTEBOOK), notebook_sha256=digest(NOTEBOOK),
        extraction=dict(
            selected_time_interval_ps=[0.0, END_PS],
            common_available_time_interval_ps=[max(a[0, 0] for a in arrays), END_PS],
            intensity="(column4^2 + column5^2)/Esat^2; columns counted from 1",
            onset="First ascending crossing of 0.81 times each record's own maximum intensity, linearly interpolated",
            late_mean="Arithmetic mean of saved intensity samples in 150 <= t <= 200 ps; not a time-weighted mean",
            normalization_source="Esat=0.317931 in TimeBehavior.ipynb and the gain-matched probe-sweep header",
        ),
        records=records,
        absolute_onset_difference_ps=difference,
        onset_difference_relative_to_dt1=difference / records[0]["t_on_ps"],
        absolute_late_mean_difference=abs(records[1]["late_mean_I_over_Isat"] - records[0]["late_mean_I_over_Isat"]),
        limitations=[
            "Step labels come from commented file paths in notebook cells 0 and 1; those cells actively read different TimeBehavior files. The present metrics are freshly computed from the Representative records.",
            "The headers do not record integration steps or the full generating configuration. The exact modified simulation source is unavailable.",
            "Neither record includes the t=0 state. Saved timestamps and field components are rounded; printed sample cadence is not the integration step.",
            "This single fixed-pump comparison does not establish convergence order, a numerical error bound, or the accuracy of gain-switched pulses.",
            "The 100000-sample tail spans different physical intervals in the two records; the primary comparison uses the common 150-200 ps interval.",
        ],
    )


def probe_checks():
    records = []
    for gain, name in [(1.05, "1.05Gth.dat"), (1.10, "1.1Gth.dat"), (1.15, "1.15Gth.dat")]:
        path = PROBE / name
        values = header(path)
        if not {"G", "ome", "Esat"} <= values.keys() or values["Esat"] <= 0:
            raise ValueError(f"Missing/invalid probe header: {path}")
        data = np.loadtxt(path)
        if data.shape != (500, 4) or not np.isfinite(data).all():
            raise ValueError(f"Unexpected or nonfinite probe data: {path}")
        if np.any(data[:, 2] <= 0) or np.any(np.diff(data[:, 2]) <= 0) or np.any(data[:, 3] <= 0):
            raise ValueError(f"Invalid probe amplitudes or maximum response: {path}")
        if not np.allclose(data[:, 1] / values["Esat"], data[:, 2], rtol=3e-5, atol=0):
            raise ValueError(f"Probe normalization disagrees with header: {path}")
        intensity = (data[:, 3] / values["Esat"]) ** 2
        def point(index):
            return dict(
                saved_row_index=int(index), E0_over_Esat=float(data[index, 2]),
                maximum_dipolar_amplitude=float(data[index, 3]),
                maximum_I_over_Isat=float(intensity[index]),
                supplied_t_on_ps=float(data[index, 0]),
            )
        records.append(dict(
            source=relative(path), sha256=digest(path), header=values,
            gain_ratio=gain, rows=len(data), weak_probe_endpoint=point(0),
            largest_probe_endpoint=point(len(data) - 1),
            largest_recorded_response=point(int(np.argmax(intensity))),
        ))
    if len({record["header"]["ome"] for record in records}) != 1:
        raise ValueError("Probe records do not share study frequency")
    thresholds = [record["header"]["G"] / record["gain_ratio"] for record in records]
    if not np.allclose(thresholds, np.median(thresholds), rtol=1e-5, atol=0):
        raise ValueError("Probe gain headers disagree with the intended gain ratios")
    if records[1]["header"]["Esat"] != ESAT:
        raise ValueError("Representative and probe saturation fields disagree")
    return dict(
        extraction="E0/Esat is column3; the normalized maximum dipolar response is (column4/header Esat)^2; columns counted from 1",
        physical_scope="Maximum dipolar response under a continuous coherent probe; not a calibrated emitted power, autonomous emission, or necessarily a stationary intensity",
        records=records,
    )


def main():
    result = dict(
        script=relative(Path(__file__).resolve()), schema_version=1,
        operation="Post-processing and provenance checks of existing data only; no simulation or figure generation",
        representative_fixed_pump=representative_checks(),
        probe_peak_sweeps=probe_checks(),
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    comparison = result["representative_fixed_pump"]
    print(f"Wrote {relative(OUTPUT)}")
    print(f"Representative onset difference: {comparison['absolute_onset_difference_ps']:.9g} ps")
    print("Checked two representative records and three probe-peak sweeps.")


if __name__ == "__main__":
    main()
