"""Shared paths for the OMEx gain-switching analysis; no scientific settings here."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/output/gain_switching_corrected_20261114'
PREFIX = '[OMEx]_'
IMAGES = ROOT / 'img/output/[OMEx]_gain_switching'
DIAGNOSTICS = IMAGES / 'diagnostics'
DRIVER_SOURCE = ROOT / 'src/[OMEx]_gain_switching.cxx'
BINARY = ROOT / 'bin/[OMEx]_gain_switching'
