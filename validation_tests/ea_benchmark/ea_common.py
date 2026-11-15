"""
Shared helpers for the Environment Agency 2D benchmark cases (SC120002,
Neelz and Pender 2013): data access, input readers, gauge sampling,
baselines and plots.

Each case directory (test3_obstruction, ...) imports this module through
a sys.path entry for its parent directory.
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def data_path(test, filename):
    """Path of a file in data/<test>/, downloading the data set if needed."""
    path = os.path.join(DATA_DIR, test, filename)
    if not os.path.exists(path):
        sys.path.insert(0, HERE)
        import data_download
        data_download.download()
        data_download.unpack()
    return path


def read_table(path):
    """Read a two-column CSV with a header row; return (header, array)."""
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    header = rows[0]
    values = np.array([[float(v) for v in row[:2]] for row in rows[1:] if row],
                      dtype=float)
    return header, values


def read_inflow(path):
    """Inflow hydrograph as (time [s], Q [m^3/s]); minutes are converted."""
    header, values = read_table(path)
    t, q = values[:, 0], values[:, 1]
    if 'min' in header[0].lower():
        t = 60.0 * t
    return t, q


def inflow_function(t_table, q_table):
    """Q(t) by linear interpolation, zero outside the table (as specified)."""
    def Q(t):
        return float(np.interp(t, t_table, q_table, left=0.0, right=0.0))
    return Q


def inflow_volume(t_table, q_table, t_end):
    """Volume let in by the hydrograph up to t_end [m^3]."""
    t = np.linspace(0.0, t_end, 20001)
    return np.trapezoid(np.interp(t, t_table, q_table, left=0.0, right=0.0), t)


def read_points(path):
    """Output points as an ordered dict {name: (x, y)}."""
    points = {}
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    for row in rows[1:]:
        if row:
            points[row[0].strip()] = (float(row[1]), float(row[2]))
    return points


def read_polygon(path):
    """The first polygon (or polyline) of an ESRI shapefile, as a list of [x, y].

    A polygon is returned unclosed; a polyline keeps all its points.

    Reads the .shp directly (format: ESRI Shapefile Technical Description,
    1998), so no shapefile package is needed.
    """
    import struct
    with open(path, 'rb') as f:
        data = f.read()
    off = 100 + 8                                  # file header, record header
    shape_type, = struct.unpack('<i', data[off:off + 4])
    assert shape_type in (3, 13, 23, 5, 15, 25), 'not a polygon or polyline shapefile: %s' % path
    n_parts, n_points = struct.unpack('<2i', data[off + 36:off + 44])
    parts = struct.unpack('<%di' % n_parts, data[off + 44:off + 44 + 4 * n_parts])
    start = off + 44 + 4 * n_parts
    xy = np.frombuffer(data, '<f8', 2 * n_points, start).reshape(-1, 2)
    end = parts[1] if n_parts > 1 else n_points
    pts = [list(map(float, p)) for p in xy[parts[0]:end]]
    if shape_type in (5, 15, 25) and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def read_asc(path):
    """ESRI ASCII grid as (x_centres, y_centres ascending, z[ny, nx])."""
    header = {}
    with open(path) as f:
        for _ in range(6):
            key, value = f.readline().split()
            header[key.lower()] = float(value)
    z = np.loadtxt(path, skiprows=6)[::-1, :]          # row 0 = south
    nx, ny = int(header['ncols']), int(header['nrows'])
    cs = header['cellsize']
    x = header['xllcorner'] + cs * (np.arange(nx) + 0.5)
    y = header['yllcorner'] + cs * (np.arange(ny) + 0.5)
    return x, y, z


def grid_function(x_grid, y_grid, z):
    """Bilinear interpolation of a grid, for domain.set_quantity."""
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((y_grid, x_grid), z, bounds_error=False,
                                     fill_value=None)

    def f(x, y):
        return interp(np.column_stack((y, x)))
    return f


def absolute(domain, f):
    """Wrap f(x, y) in absolute coordinates for domain.set_quantity.

    ANUGA passes quantity functions coordinates relative to the mesh's
    geo-reference origin, which is the bounding polygon's lower-left
    corner, not (0, 0).
    """
    xll = domain.geo_reference.get_xllcorner()
    yll = domain.geo_reference.get_yllcorner()

    def g(x, y):
        return f(np.asarray(x) + xll, np.asarray(y) + yll)
    return g


# ---------------------------------------------------------------------------
# Gauges, sampled during evolve (serial runs)
# ---------------------------------------------------------------------------

class Gauges:
    """Record stage, depth and speed at the triangles holding each point."""

    def __init__(self, domain, points):
        self.domain = domain
        self.names = list(points)
        # Points are absolute, as get_triangle_containing_point expects.
        self.ids = [domain.get_triangle_containing_point(points[n])
                    for n in self.names]
        self.time = []
        self.stage = []
        self.depth = []
        self.speed = []

    def sample(self):
        d = self.domain
        stage = d.quantities['stage'].centroid_values[self.ids]
        elev = d.quantities['elevation'].centroid_values[self.ids]
        xmom = d.quantities['xmomentum'].centroid_values[self.ids]
        ymom = d.quantities['ymomentum'].centroid_values[self.ids]
        h = np.maximum(stage - elev, 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            speed = np.where(h > 1.0e-3, np.hypot(xmom, ymom) / h, 0.0)
        self.time.append(d.get_time())
        self.stage.append(stage.copy())
        self.depth.append(h)
        self.speed.append(speed)

    def save(self, path):
        """One CSV: time, then stage/depth/speed per gauge."""
        cols = ['time']
        for n in self.names:
            cols += ['stage_%s' % n, 'depth_%s' % n, 'speed_%s' % n]
        data = [np.array(self.time)]
        for k in range(len(self.names)):
            data += [np.array(self.stage)[:, k], np.array(self.depth)[:, k],
                     np.array(self.speed)[:, k]]
        np.savetxt(path, np.column_stack(data), delimiter=',',
                   header=','.join(cols), comments='', fmt='%.6g')


def load_gauges(path):
    """Read a gauge CSV written by Gauges.save into a dict of arrays."""
    with open(path) as f:
        cols = f.readline().strip().split(',')
    data = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    return {c: data[:, i] for i, c in enumerate(cols)}


def domain_volume(domain):
    """Water volume in the domain [m^3]."""
    stage = domain.quantities['stage'].centroid_values
    elev = domain.quantities['elevation'].centroid_values
    return float(np.sum(np.maximum(stage - elev, 0.0) * domain.areas))


# ---------------------------------------------------------------------------
# Baselines and plots
# ---------------------------------------------------------------------------

def compare_baseline(result, baseline, columns, atol):
    """Max abs difference per column between two gauge dicts, on result times.

    Returns {column: max |result - baseline|} and the columns over atol.
    """
    diffs = {}
    for c in columns:
        b = np.interp(result['time'], baseline['time'], baseline[c])
        diffs[c] = float(np.max(np.abs(result[c] - b)))
    failed = [c for c in columns if diffs[c] > atol[c.split('_')[0]]]
    return diffs, failed


def check_baseline(gauge_csv, stem, alg, columns, atol, indent=''):
    """Compare a run's gauges with baseline_<stem>_<alg>.csv.

    EA_REGENERATE_BASELINE=1 writes (or rewrites) the baseline instead.
    With no baseline for this algorithm, nothing is compared. Returns the
    columns that moved by more than atol (keyed by quantity: 'stage', ...).
    """
    import shutil
    baseline = 'baseline_%s_%s.csv' % (stem, alg)
    if os.environ.get('EA_REGENERATE_BASELINE') == '1':
        shutil.copy(gauge_csv, baseline)
        print(indent + 'Baseline written: %s' % baseline)
        return []
    if not os.path.exists(baseline):
        print(indent + 'No baseline for %s; not compared' % alg)
        return []
    diffs, failed = compare_baseline(load_gauges(gauge_csv), load_gauges(baseline),
                                     columns, atol)
    for c, d in diffs.items():
        print(indent + 'Baseline max |diff| %s: %.4f' % (c, d))
    return failed


def plot_gauges(result, names, quantity, filename, title, ylabel,
                time_scale=1.0, time_label='Time (s)', baseline=None, ncol=2,
                tmax=None, label='ANUGA'):
    """One panel per gauge: ANUGA (and the stored baseline, if given)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    n = len(names)
    ncol = min(ncol, n)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 3.0 * nrow),
                             squeeze=False, sharex=True)
    t = result['time'] / time_scale
    for k, name in enumerate(names):
        ax = axes[k // ncol][k % ncol]
        col = '%s_%s' % (quantity, name)
        if baseline is not None:
            ax.plot(baseline['time'] / time_scale, baseline[col], color='0.6',
                    lw=2.5, label='stored baseline')
        ax.plot(t, result[col], color='C0', lw=1.2, label=label)
        ax.set_title('Point %s' % name if name.isdigit() else name, fontsize=10)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel(time_label)
        if tmax is not None:
            ax.set_xlim(0, tmax / time_scale)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)
    axes[0][0].legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(filename, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Checks: one list, shared by validate_*.py and the report
# ---------------------------------------------------------------------------

class Checks:
    """Named pass/fail checks with the value found and the criterion.

    validate_<test>.py asserts on them; plot_results.py draws them as a
    table for the report and exits non-zero if any failed, so the release
    report marks the case failed.
    """

    def __init__(self):
        self.rows = []          # (name, value, criterion, ok)

    def add(self, name, value, criterion, ok):
        self.rows.append((name, value, criterion, bool(ok)))

    def failed(self):
        return [r for r in self.rows if not r[3]]

    def report(self, indent=''):
        for name, value, criterion, ok in self.rows:
            print(indent + '%-4s %s: %s (%s)' % ('ok' if ok else 'FAIL', name, value, criterion))

    def figure(self, filename, title):
        """The checks as a table image (one per algorithm in the report)."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        n = len(self.rows)
        fig, ax = plt.subplots(figsize=(9.0, 0.3 * n + 0.5))
        ax.axis('off')
        cells = [[name, value, criterion, 'pass' if ok else 'FAIL']
                 for name, value, criterion, ok in self.rows]
        table = ax.table(cellText=cells, colLabels=['Check', 'ANUGA', 'Criterion', ''],
                         colWidths=[0.36, 0.22, 0.34, 0.08], loc='center', cellLoc='left')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.auto_set_column_width(range(4))     # fit the text, which can be long
        table.scale(1.0, 1.35)
        for (r, c), cell in table.get_celld().items():
            cell.set_edgecolor('0.8')
            if r == 0:
                cell.set_text_props(weight='bold')
            elif c == 3:
                cell.set_text_props(color='#1a7f37' if self.rows[r - 1][3] else '#c62828',
                                    weight='bold')
        ax.set_title(title, fontsize=11)
        fig.savefig(filename, dpi=120, bbox_inches='tight')
        plt.close(fig)


def baseline_check(checks, gauge_csv, stem, alg, columns, atol):
    """Add the stored-baseline comparison to checks (see check_baseline)."""
    baseline = 'baseline_%s_%s.csv' % (stem, alg)
    tol = list(atol.values())[0]
    if os.environ.get('EA_REGENERATE_BASELINE') == '1':
        check_baseline(gauge_csv, stem, alg, columns, atol)
        checks.add('Levels against stored baseline', 'rewritten',
                   '|diff| < %g m' % tol, True)
    elif not os.path.exists(baseline):
        checks.add('Levels against stored baseline', 'no baseline for %s' % alg,
                   '|diff| < %g m' % tol, True)
    else:
        diffs, failed = compare_baseline(load_gauges(gauge_csv), load_gauges(baseline),
                                         columns, atol)
        checks.add('Levels against stored baseline', 'max |diff| %.4f m' % max(diffs.values()),
                   '|diff| < %g m' % tol, not failed)


def read_volume(path):
    """The key/value file each run writes next to its gauges."""
    with open(path) as f:
        return {k: float(v) for k, v in (line.split() for line in f)}
