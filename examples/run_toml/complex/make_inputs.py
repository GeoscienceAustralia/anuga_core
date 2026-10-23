#!/usr/bin/env python3
"""Generate the input data files for the complex anuga_toml_run example.

The committed CSVs in this directory were produced by this script; re-run it
(``python make_inputs.py``) to regenerate them or to experiment with the
geometry.  No ANUGA import is needed — it only writes small CSV files.

Domain: a 400 m x 200 m plain that slopes from 3 m elevation in the west down
to 0 m at the eastern outlet.  Water enters from an inlet hydrograph in the
west, rain falls over the whole domain, and everything drains out through a
prescribed-stage boundary on the east.
"""

import numpy as np

LENX, LENY = 400.0, 200.0


def write_polygon(path, points):
    """Plain x,y polygon/line CSV (no header) — the format read_polygon expects."""
    with open(path, 'w') as f:
        for x, y in points:
            f.write(f'{x:g},{y:g}\n')


def write_timeseries(path, header, rows):
    """Timeseries CSV WITH a header row (genfromtxt skip_header=1)."""
    with open(path, 'w') as f:
        f.write(header + '\n')
        for t, v in rows:
            f.write(f'{t:g},{v:g}\n')


# --- Bounding polygon: rectangle, vertices counter-clockwise from the origin --
# Edge order (0-based, between consecutive vertices, plus the closing edge):
#   0: (0,0)->(400,0)     south
#   1: (400,0)->(400,200) east   (outlet)
#   2: (400,200)->(0,200) north
#   3: (0,200)->(0,0)     west
write_polygon('bounding_polygon.csv',
              [(0, 0), (LENX, 0), (LENX, LENY), (0, LENY)])

# --- Sloped bed as x,y,z points (read via nearest-neighbour interpolation) ----
# z = 3 m in the west, dropping linearly to 0 m at the eastern outlet.
xs = np.linspace(0, LENX, 21)
ys = np.linspace(0, LENY, 11)
with open('elevation.csv', 'w') as f:
    f.write('x,y,z\n')
    for x in xs:
        z = 3.0 * (1.0 - x / LENX)
        for y in ys:
            f.write(f'{x:g},{y:g},{z:.4f}\n')

# --- Interior refinement region: a central band resolved more finely ----------
write_polygon('refine_channel.csv',
              [(40, 70), (360, 70), (360, 130), (40, 130)])

# --- Inlet cross-section line in the west (water is injected across this line) -
write_polygon('inlet_line.csv', [(20, 60), (20, 140)])

# --- Inlet hydrograph: ramp up to 8 m^3/s, hold, then recede ------------------
write_timeseries('inlet_hydrograph.csv', 'time,discharge',
                 [(0, 0), (120, 8), (600, 8), (900, 2), (1800, 0)])

# --- Rainfall: a 40 mm/hr burst over the first 10 minutes ---------------------
write_timeseries('rain.csv', 'time,rain_mm_hr',
                 [(0, 0), (60, 40), (600, 40), (660, 0), (1800, 0)])

# --- Downstream (east) boundary stage: hold the sea/outlet at 0.5 m -----------
write_timeseries('downstream_stage.csv', 'time,stage',
                 [(0, 0.5), (1800, 0.5)])

print('Wrote: bounding_polygon.csv, elevation.csv, refine_channel.csv, '
      'inlet_line.csv, inlet_hydrograph.csv, rain.csv, downstream_stage.csv')
