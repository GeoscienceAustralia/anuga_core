"""
EA 2D benchmark Test 6A: dam break past a building, laboratory scale.

The IMPACT isolated-building flume of Soares-Frazao and Zech (2002): a
3.6 m wide channel with a dam that has a 1 m opening, 0.4 m of water
behind it and 0.02 m in front, and a 0.8 m x 0.4 m building set at an
angle downstream. All boundaries are vertical walls. Manning n = 0.01,
0.1 m resolution, run to 2 min. The set-up is shared with Test 6B in
../ea_dambreak.py, which describes the geometry and the EA_6A_MAX_AREA and
EA_6A_BUILDING settings.

Writes test6a_gauges.csv (level, depth and speed at G1-G6 every 0.1 s),
test6a_volume.txt (volumes and the largest speed in the domain),
test6a_geometry.json and test6a.sww every second.
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_dambreak  # noqa: E402

args = anuga.get_args()
ea_dambreak.run('6A', args.alg, args.verbose)
