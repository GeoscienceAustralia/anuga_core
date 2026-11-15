"""Draw the figures on the "ARR design rainfall data" appendix.

Writes each figure as .svg (HTML) and .pdf (LaTeX). The rendered files are
committed, so building the docs does not need the data; only rerunning this does::

    python make_arr_rainfall_figures.py Point_ECsouth.zip

The two pattern figures use the ARR Data Hub's temporal patterns, downloaded as the
appendix describes. The two IFD figures use invented depths and an invented coast.
They are written to a zip file in the Bureau's format and read back through
Arr_ifd_rain, so what they show is what the class does. They are not real rainfall,
and no Bureau grid is used here.
"""
import argparse
import tempfile
import warnings
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from anuga.rain import ARR_point_rainfall_patterns, Arr_ifd_rain

DURATION = 1440               # the 24-hour storm in the appendix's examples
AEP = 1.0                     # %
DURATIONS = [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 45, 60, 90, 120, 180, 270, 360,
             540, 720, 1080, 1440, 1800, 2160, 2880, 4320, 5760, 7200, 8640, 10080]

# The invented grid. The depth at the reference point is REF_DEPTH for 24 h, 1% AEP
XLL, YLL, CELL, NCOLS, NROWS = 150.0, -30.0, 0.025, 20, 20
REF_LON, REF_LAT = 150.2, -29.8
REF_DEPTH = 300.0

AEPS = [50, 10, 2, 1, 0.1]    # %, with a validated light-to-dark ordinal ramp
RAMP = ['#86b6ef', '#3987e5', '#256abf', '#184f95', '#0d366b']
SEQUENTIAL = ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
ACCENT = '#2a78d6'
CONTEXT = '#a8a79f'
INK = '#0b0b0b'
INK_2 = '#52514e'
GRID = '#e4e3df'
OCEAN = '#e4e3df'

DURATION_TICKS = {1: '1 min', 10: '10 min', 60: '1 h', 360: '6 h',
                  1440: '24 h', 10080: '7 days'}


def style():
    plt.rcParams.update({
        'font.size': 9, 'axes.edgecolor': GRID, 'axes.labelcolor': INK_2,
        'axes.linewidth': 0.8, 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': GRID, 'grid.linewidth': 0.6, 'xtick.color': INK_2,
        'ytick.color': INK_2, 'xtick.major.size': 0, 'ytick.major.size': 0,
        'xtick.minor.size': 0, 'ytick.minor.size': 0, 'text.color': INK,
        'svg.hashsalt': 'anuga', 'svg.fonttype': 'path',
    })


def save(fig, stem, out_dir):
    fig.savefig(out_dir / f'{stem}.svg', metadata={'Date': None})
    fig.savefig(out_dir / f'{stem}.pdf', metadata={'CreationDate': None})
    plt.close(fig)


def read_patterns(zip_path):
    """Return {EventID: percent of the depth in each step} for the rare 24 h set."""
    prp = ARR_point_rainfall_patterns(str(zip_path), 'ECsouth')
    out = {}
    for i, row in enumerate(prp.linesInc):
        if i > 0 and int(row[1]) == DURATION and row[4].strip() == 'rare':
            if int(row[2]) != 60:
                raise ValueError(f'pattern {row[0]}: expected hourly steps')
            out[int(row[0])] = np.array([float(v) for v in row[5:] if v.strip()])
    if len(out) != 10:
        raise ValueError(f'expected ten patterns, found {len(out)}')
    return out


def coast_lon(lat):
    """Longitude of the invented coast."""
    return 150.38 + 0.05 * np.sin(14.0 * (lat + 30.0))


def synthetic_depth(dur, aep, lon, lat):
    """Invented depth (mm), smooth in duration, AEP and position."""
    x = np.log(dur / DURATION)
    gumbel = -np.log(-np.log(1.0 - aep / 100.0))
    gumbel_ref = -np.log(-np.log(1.0 - AEP / 100.0))
    hill = np.exp(-((lon - 150.28) ** 2 + (lat + 29.62) ** 2) / 0.02)
    hill_ref = np.exp(-((REF_LON - 150.28) ** 2 + (REF_LAT + 29.62) ** 2) / 0.02)
    return (REF_DEPTH * np.exp(0.38 * x - 0.015 * x ** 2)
            * (2.0 + gumbel) / (2.0 + gumbel_ref)
            * (1.0 + 0.6 * hill) / (1.0 + 0.6 * hill_ref))


def write_synthetic_zip(directory):
    """Write the invented grids in the format of the Bureau's download."""
    lons = XLL + (np.arange(NCOLS) + 0.5) * CELL
    lats = YLL + (np.arange(NROWS) + 0.5) * CELL
    with zipfile.ZipFile(Path(directory) / 'synthetic.zip', 'w') as archive:
        for dur in DURATIONS:
            for aep in AEPS:
                lines = [f'ncols {NCOLS}', f'nrows {NROWS}', f'xllcorner {XLL}',
                         f'yllcorner {YLL}', f'cellsize {CELL}', 'NODATA_value -9999']
                for lat in lats[::-1]:     # the first row is the northernmost
                    depth = synthetic_depth(dur, aep, lons, lat)
                    depth[lons > coast_lon(lat)] = 0.0
                    lines.append(' '.join(f'{d:.1f}' for d in depth))
                name = f'catchment_depth_{dur}min_{round(aep * 1000):05d}aep.txt.asc'
                archive.writestr(name, '\n'.join(lines) + '\n')


def ifd_curves(ifd, out_dir):
    durations = sorted(int(d) for d in set(ifd.dur_list))
    fig, ax = plt.subplots(figsize=(7, 4.2), layout='constrained')
    ends = []
    for aep, colour in zip(AEPS, RAMP):
        depth = [ifd.open_grd(Dur=d, aep_percent=aep).get_rain_at_point(REF_LON, REF_LAT)
                 for d in durations]
        label = f'{aep:g}% AEP'
        ax.plot(durations, depth, color=colour, lw=1.8, label=label)
        ends.append((depth[-1], label))
    # the line ends can be closer than the text is tall, so spread the labels
    previous = 0
    for end, label in ends:
        y = max(end, previous * 1.2)
        ax.annotate(label, (durations[-1], y), xytext=(6, 0),
                    textcoords='offset points', va='center', color=INK_2)
        previous = y
    depth = ifd.open_grd(Dur=DURATION, aep_percent=AEP).get_rain_at_point(REF_LON, REF_LAT)
    ax.plot([DURATION], [depth], 'o', ms=6, color=RAMP[AEPS.index(AEP)],
            mec='white', mew=1.2, zorder=3)
    ax.annotate(f'24 h, {AEP:g}% AEP: {depth:.0f} mm', (DURATION, depth),
                xytext=(DURATION, 45), ha='center', va='top',
                arrowprops={'arrowstyle': '-', 'color': INK_2, 'lw': 0.8,
                            'shrinkB': 4})
    ax.set(xscale='log', yscale='log', xlabel='Storm duration',
           ylabel='Rainfall depth (mm)')
    ax.set_xticks(list(DURATION_TICKS), list(DURATION_TICKS.values()))
    ax.set_yticks([5, 10, 20, 50, 100, 200, 500, 1000],
                  ['5', '10', '20', '50', '100', '200', '500', '1000'])
    ax.minorticks_off()
    ax.grid(True, which='major')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], frameon=False, loc='upper left',
              labelcolor=INK_2)
    save(fig, 'arr_ifd_curves', out_dir)


def ifd_grid(ifd, out_dir):
    grd = ifd.open_grd(Dur=DURATION, aep_percent=AEP)
    lons = XLL + (np.arange(NCOLS) + 0.5) * CELL
    lats = YLL + (np.arange(NROWS) + 0.5) * CELL
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')    # the ocean cells warn, as they should
        depth = grd.get_rain_at_points(lon_grid.ravel(), lat_grid.ravel())
    depth = depth.reshape(lon_grid.shape)

    fig, ax = plt.subplots(figsize=(6.2, 4.8), layout='constrained')
    ax.set_facecolor(OCEAN)
    ax.grid(False)
    edges_x = XLL + np.arange(NCOLS + 1) * CELL
    edges_y = YLL + np.arange(NROWS + 1) * CELL
    mesh = ax.pcolormesh(edges_x, edges_y, np.ma.masked_equal(depth, 0.0),
                         cmap=LinearSegmentedColormap.from_list('depth', SEQUENTIAL),
                         edgecolors='white', linewidth=0.3)
    coast_lats = np.linspace(YLL, YLL + NROWS * CELL, 200)
    ax.plot(coast_lon(coast_lats), coast_lats, color=INK, lw=1.2)

    # A point on land, in the cell whose centre is in the ocean with the most land
    first_ocean = [int(np.argmax(depth[i] == 0.0)) for i in range(NROWS)]
    land_width = [coast_lon(lats[i]) - edges_x[first_ocean[i]] for i in range(NROWS)]
    row = int(np.argmax(land_width))
    col = first_ocean[row]
    shore_lon = edges_x[col] + 0.5 * land_width[row]
    points = [(REF_LON, REF_LAT, (-8, -16), 'right'),
              (shore_lon, lats[row], (10, 12), 'left')]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for lon, lat, offset, align in points:
            value = grd.get_rain_at_point(lon, lat)
            ax.plot(lon, lat, 'o', ms=6, mfc='white', mec=INK, mew=1.2, zorder=3)
            ax.annotate(f'{value:.0f} mm', (lon, lat), xytext=offset, ha=align,
                        textcoords='offset points',
                        bbox={'fc': 'white', 'ec': 'none', 'pad': 1.5})
    assert depth[row, col] == 0.0 and edges_x[col] < shore_lon < coast_lon(lats[row])

    ax.set(xlabel='Longitude', ylabel='Latitude',
           xlim=(edges_x[0], edges_x[-1]), ylim=(edges_y[0], edges_y[-1]),
           aspect=1.0 / np.cos(np.radians(-29.75)))
    ax.tick_params(axis='x', pad=6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(mesh, ax=ax, shrink=0.8, label='Depth (mm)')
    bar.outline.set_visible(False)
    save(fig, 'arr_ifd_grid', out_dir)


def pattern_multiples(patterns, out_dir):
    fig, axes = plt.subplots(2, 5, figsize=(7.4, 3.4), sharex=True, sharey=True,
                             layout='constrained')
    hours = np.arange(24) + 0.5
    for ax, (event, percent) in zip(axes.flat, patterns.items()):
        ax.bar(hours, percent, width=0.8, color=ACCENT)
        ax.set_title(f'{event}', fontsize=9, color=INK_2, loc='left', pad=3)
        ax.grid(axis='x', visible=False)
        ax.set_xlim(0, 24)
        ax.set_xticks([0, 12, 24])
    fig.supxlabel('Hours from the start of the storm', fontsize=9, color=INK_2)
    fig.supylabel('Rain in the hour (% of the storm depth)', fontsize=9, color=INK_2)
    save(fig, 'arr_pattern_ensemble', out_dir)


def pattern_cumulative(patterns, out_dir):
    fig, ax = plt.subplots(figsize=(7, 3.6), layout='constrained')
    hours = np.arange(25)
    events = list(patterns)
    for event in events[1:]:
        ax.plot(hours, np.r_[0, np.cumsum(patterns[event])], lw=1, color=CONTEXT,
                label='The other nine rare patterns' if event == events[1] else '_')
    ax.plot(hours, np.r_[0, np.cumsum(patterns[events[0]])], lw=2, color=ACCENT,
            label=f'Pattern {events[0]}, rows[0] in the example above')
    ax.set(xlim=(0, 24), ylim=(0, 100), xlabel='Hours from the start of the storm',
           ylabel='Rain so far (% of the storm depth)')
    ax.set_xticks(range(0, 25, 6))
    ax.legend(frameon=False, loc='upper left', labelcolor=INK_2)
    save(fig, 'arr_pattern_cumulative', out_dir)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('patterns_zip')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    style()
    with tempfile.TemporaryDirectory() as directory:
        write_synthetic_zip(directory)
        ifd = Arr_ifd_rain(directory, Lat=REF_LAT, Lon=REF_LON)
        ifd_curves(ifd, args.out)
        ifd_grid(ifd, args.out)
    patterns = read_patterns(args.patterns_zip)
    pattern_multiples(patterns, args.out)
    pattern_cumulative(patterns, args.out)
    # numbers quoted in the appendix text
    peak = [p.max() for p in patterns.values()]
    first_half = [p[:12].sum() for p in patterns.values()]
    print(f'wettest hour: {min(peak):.1f} to {max(peak):.1f} % of the depth')
    print(f'first 12 h: {min(first_half):.1f} to {max(first_half):.1f} %')


if __name__ == '__main__':
    main()
