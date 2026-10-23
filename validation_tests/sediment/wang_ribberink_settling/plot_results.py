import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import setup_diagrams          # noqa: E402  the shared flume schematics

setup_diagrams.settling_flume('setup.png')

"""Depth-averaged concentration along the test section: measured against
every adaptation closure that has been run."""
import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from validate_settling import compare, measured

LABELS = {'none': 'instantaneous exchange [D-1], the old default',
          'armanini': 'Galappatti lag [D-3]',
          'carried': 'carried near-bed ratio [D-4] (same as [D-1] here)',
          'two_layer': 'two-layer + velocity split [D-5], the default'}
for run in (1, 2):
    files = sorted(glob.glob('settling_run%d_*.npz' % run))
    if not files:
        continue
    fig, ax = plt.subplots(figsize=(9, 4.5))
    xm, cm = measured(run)
    ax.plot(xm, cm, 'kx', ms=8, label='measured (Wang & Ribberink 1986, Table 1)')
    for f in files:
        adapt = f[len('settling_run%d_' % run):-4]
        r = np.load(f)
        x, c = r['x'], r['c'] * 2650.0e3
        bins = np.arange(-0.05, 20.05, 0.1); idx = np.digitize(x, bins) - 1
        xc = np.array([x[idx == k].mean() for k in range(len(bins) - 1) if (idx == k).any()])
        cc = np.array([c[idx == k].mean() for k in range(len(bins) - 1) if (idx == k).any()])
        ax.plot(xc, cc, lw=1.5, label=LABELS.get(adapt, adapt))
    ax.axvline(16.0, color='0.7', ls=':')
    ax.set_xlim(0, 17); ax.set_ylim(0, 160); ax.set_xlabel('x from the start of the perforated bed (m)')
    ax.set_ylabel('depth-averaged concentration (ppm)'); ax.grid(alpha=0.3); ax.legend(fontsize='small')
    ax.set_title('Wang & Ribberink settling flume, run %d' % run)
    fig.tight_layout(); fig.savefig('concentration_run%d.png' % run, dpi=130); plt.close(fig)
