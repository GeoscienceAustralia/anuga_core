"""Shared fixtures and helpers for behaviour_only regression tests."""

import os
import sys
import subprocess

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Directory fixture — each test runs in its own subdirectory so that
# simulation scripts find their input files and write .sww files locally.
# ---------------------------------------------------------------------------

@pytest.fixture()
def simdir(request):
    """Change cwd to the directory containing the test file, restore after."""
    test_dir = os.path.dirname(os.path.abspath(request.fspath))
    orig = os.getcwd()
    os.chdir(test_dir)
    yield test_dir
    os.chdir(orig)


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def run_sim(script, alg='DE1', np_procs=1):
    """Run a simulation script as a subprocess.  Returns the exit code."""
    cmd = [sys.executable, script, '-alg', alg, '-np', str(np_procs)]
    result = subprocess.run(cmd)
    return result.returncode


def ensure_sww(sww_file, script, alg='DE1', np_procs=1):
    """Clean up stale files and run the simulation to produce *sww_file*.

    Always re-runs the simulation so that the .sww used for regression
    comparison is always freshly generated from the current code.  Old .msh /
    .sww files are removed first to avoid mesh-generation interference.
    """
    for f in os.listdir('.'):
        if f.endswith(('.sww', '.msh', '.stdout', '.png')):
            os.remove(f)

    rc = run_sim(script, alg=alg, np_procs=np_procs)
    assert rc == 0, f'{script} failed with exit code {rc}'


# ---------------------------------------------------------------------------
# SWW metric extraction
# ---------------------------------------------------------------------------

def sww_metrics(sww_file):
    """Extract scalar regression metrics from a .sww file.

    Returns a dict of 1-D numpy arrays suitable for num_regression.check().
    """
    from anuga.utilities import plot_utils as util

    p = util.get_output(sww_file, 0.001)
    pc = util.get_centroids(p, velocity_extrapolation=True)

    stage = pc.stage          # shape (T, N)
    xmom  = pc.xmom           # shape (T, N)

    return {
        'peak_max_stage':  np.array([float(stage.max())]),
        'peak_min_stage':  np.array([float(stage.min())]),
        'final_max_stage': np.array([float(stage[-1].max())]),
        'final_min_stage': np.array([float(stage[-1].min())]),
        'peak_max_xmom':   np.array([float(xmom.max())]),
        'peak_min_xmom':   np.array([float(xmom.min())]),
    }


# ---------------------------------------------------------------------------
# HEC-RAS comparison helper
# ---------------------------------------------------------------------------

def hecras_stage_at_gauges(pc, gauges_csv, x_channel=15.0):
    """Return ANUGA stage timeseries at the y-positions of HEC-RAS gauges.

    Parameters
    ----------
    pc : centroids object from plot_utils.get_centroids()
    gauges_csv : path to the HEC-RAS gauges CSV file (3-row header)
    x_channel : x-coordinate of the ANUGA channel centreline (default 15 m)

    Returns
    -------
    ras_stations : list of station names (strings)
    ras_stage    : (T_ras, n_stations) array of HEC-RAS stage
    anuga_stage  : (T_anuga, n_stations) array of ANUGA stage
    """
    # Read HEC-RAS gauge file
    with open(gauges_csv) as f:
        headers   = f.readline().strip().split(',')
        _units    = f.readline()   # units row (skip)
        _blank    = f.readline()   # blank row (skip)
    ras_data = np.genfromtxt(gauges_csv, skip_header=3, delimiter=',')
    ras_time = np.linspace(0.0, 60.0 * (ras_data.shape[0] - 1), ras_data.shape[0])

    # Identify columns that have stage data (non-empty header, skip first 2 cols)
    stations = []
    col_indices = []
    for i, h in enumerate(headers):
        h = h.strip()
        if h and i >= 2:
            stations.append(h)
            col_indices.append(i)

    # For each station extract the y-position from "RIVER1 REACH1 <y>"
    ras_stage   = []
    anuga_stage = []

    for col_idx, station in zip(col_indices, stations):
        parts = station.split()
        try:
            y_station = float(parts[-1].rstrip('.*'))
        except ValueError:
            continue

        # Nearest ANUGA centroid at x=x_channel, y=y_station
        dist = (pc.x - x_channel)**2 + (pc.y - (1000.0 - y_station))**2
        anuga_idx = int(dist.argmin())

        ras_stage.append(ras_data[:, col_idx])
        anuga_stage.append(pc.stage[:, anuga_idx])

    return (stations,
            np.array(ras_stage).T,    # (T_ras,   n_stations)
            np.array(anuga_stage).T)  # (T_anuga, n_stations)


def hecras_correlation_metrics(pc, gauges_csv, x_channel=15.0):
    """Compute mean Pearson correlation and RMS error vs HEC-RAS gauges.

    Interpolates the shorter timeseries onto the longer one before comparing.
    Returns a dict of 1-D numpy arrays for num_regression.
    """
    _, ras_stage, anuga_stage = hecras_stage_at_gauges(pc, gauges_csv, x_channel)

    n_stations = ras_stage.shape[1]
    t_ras   = np.linspace(0.0, 60.0 * (ras_stage.shape[0] - 1), ras_stage.shape[0])
    t_anuga = np.linspace(0.0, t_ras[-1], anuga_stage.shape[0])

    correlations = []
    rms_errors   = []

    for i in range(n_stations):
        # Interpolate ANUGA onto HEC-RAS time points
        anuga_interp = np.interp(t_ras, t_anuga, anuga_stage[:, i])
        ref = ras_stage[:, i]

        # Pearson correlation
        if ref.std() > 1e-10 and anuga_interp.std() > 1e-10:
            r = float(np.corrcoef(ref, anuga_interp)[0, 1])
        else:
            r = 1.0  # both flat — trivially correlated
        correlations.append(r)

        # RMS error
        rms = float(np.sqrt(np.mean((anuga_interp - ref)**2)))
        rms_errors.append(rms)

    return {
        'mean_hecras_correlation': np.array([float(np.mean(correlations))]),
        'mean_hecras_rms_error':   np.array([float(np.mean(rms_errors))]),
        'min_hecras_correlation':  np.array([float(np.min(correlations))]),
    }
