"""
EA 2D benchmark Test 6B: dam break past a building, at field scale.

The IMPACT laboratory case (Soares-Frazao and Zech 2002) scaled up by 20:
a channel with a dam that has a 20 m opening, 8 m of water behind it and
0.4 m in front, and a 16 m x 8 m building set at an angle downstream. All
boundaries are vertical walls. Manning n = 0.05, 2 m resolution, run to
30 min. The set-up is shared with Test 6A in ../ea_dambreak.py, which
describes the geometry and the EA_6B_MAX_AREA and EA_6B_BUILDING settings.

Writes test6b_gauges.csv (level, depth and speed at G1-G6 every 1 s),
test6b_volume.txt (volumes and the largest speed in the domain),
test6b_geometry.json and test6b.sww every minute.
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_dambreak  # noqa: E402

args = anuga.get_args()
ea_dambreak.run('6B', args.alg, args.verbose)
