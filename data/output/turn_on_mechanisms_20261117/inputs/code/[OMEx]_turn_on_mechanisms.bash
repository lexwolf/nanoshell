#!/usr/bin/env bash
set -euo pipefail
omex_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-help}" in
 build|linear|edges|refine|analyze|report|validate|package)
  exec python3 "$omex_dir/[OMEx]_turn_on_$1.py"
  ;;
 help) printf '%s\n' 'Turn-on diagnostics: build, linear, edges, refine, analyze, report, validate, package.' 'Only refine launches new nonlinear trajectories; existing matching cases are reused.' ;;
 *) exit 2 ;;
esac
