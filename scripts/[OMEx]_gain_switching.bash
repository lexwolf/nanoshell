#!/usr/bin/env bash
set -euo pipefail
omex_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-help}" in
  build|pilot|batch|residual_batch|process_plot|residual_summary|validate|prepare_maps|pack_native|index_and_bundle)
    omex_task="$1"
    shift
    exec python3 "$omex_script_dir/[OMEx]_${omex_task}.py" "$@"
    ;;
  help)
    printf '%s\n' 'OMEx workflow: build, prepare_maps, pilot, batch, residual_batch,' \
      'process_plot, residual_summary, validate, pack_native, index_and_bundle.' \
      'Only batch/pilot/residual_batch launch time-dependent calculations.'
    ;;
  *) printf 'Unknown OMEx action: %s\n' "$1" >&2; exit 2 ;;
esac
