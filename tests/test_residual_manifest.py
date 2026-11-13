#!/usr/bin/env python3
"""Check provenance snapshots and resume rejection without compiling or integrating."""
import hashlib
import io
from contextlib import redirect_stderr
import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/scan_residual_modulation.py'


def main():
    spec = importlib.util.spec_from_file_location('residual_scan', SCRIPT)
    scan = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scan)
    real_run = subprocess.run
    def compile_only_mock(command, *args, **kwargs):
        if command[0] == 'g++':
            return subprocess.CompletedProcess(command, 0)
        return real_run(command, *args, **kwargs)

    with tempfile.TemporaryDirectory(prefix='residual-manifest-') as temp:
        base = Path(temp)
        frequencies = base / 'frequencies.json'
        frequencies.write_text(json.dumps({gain: [] for gain in scan.GAINS}))
        source = ROOT / 'data/output/convergence_validation/inputs'
        output = base / 'output'
        argv = [str(SCRIPT), '--frequencies', str(frequencies), '--input', str(source / 'input08.dat'),
                '--time', str(source / 'time.dat'), '--output', str(output)]
        # Empty frequency lists and mocked compilation guarantee no solver execution.
        with patch.object(sys, 'argv', argv), patch.object(scan.subprocess, 'run', side_effect=compile_only_mock) as build:
            scan.main()
            assert sum(call.args[0][0] == 'g++' for call in build.call_args_list) == 1
        manifest = json.loads((output / 'manifest.json').read_text())
        header = 'extern/nano_geo_matrix/include/nano_geo_matrix/core/saturation_scale.hpp'
        expected_hash = hashlib.sha256((ROOT / header).read_bytes()).hexdigest()
        assert manifest['source_sha256'][header] == expected_hash
        assert (output / 'saturation_scale.hpp').read_bytes() == (ROOT / header).read_bytes()
        assert manifest['saturation_scale']['tau1_policy'] == 'require_match'
        assert manifest['saturation_scale']['tau1_over_tau2'] == 5.
        assert manifest['saturation_scale']['convention'] == 'inversion_one_plasma_scaled_saturation_v1'
        # Simulate an old manifest without the new provenance contract.
        del manifest['saturation_scale']
        (output / 'manifest.json').write_text(json.dumps(manifest))
        errors = io.StringIO()
        with patch.object(sys, 'argv', argv), patch.object(scan.subprocess, 'run', side_effect=compile_only_mock) as build, redirect_stderr(errors):
            try:
                scan.main()
            except SystemExit as error:
                assert error.code == 2
            else:
                raise AssertionError('incompatible resume was accepted')
            assert all(call.args[0][0] != 'g++' for call in build.call_args_list)
        assert 'resume configuration/provenance mismatch' in errors.getvalue()
    print('Shared utility snapshot/hash, authority provenance and incompatible resume rejection: PASS')


if __name__ == '__main__':
    main()
