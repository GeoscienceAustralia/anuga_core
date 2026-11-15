"""
EA 2D benchmark Test 5: valley flooding after a dam failure.

A flood from a failed embankment dam runs 17 km down a real valley, which
falls at about 1:100 in its upper part and 1:1000 in its lower part, and
settles in a pond at its foot. The inflow (a trapezoid peaking at
3000 m^3/s from 10 to 20 min, back to zero at 100 min) enters along a
~260 m line across the upstream end. The modelled area is the EA's outline
polygon (19.1 km^2), all of its boundaries closed, on a dry bed. Manning
n = 0.04, 50 m resolution, run to 30 h.

EA_5_MAX_AREA sets the triangle area (default 2500 m^2, the 50 m grid);
the outputs are then named test5_a<area>_*. Finer meshes are best run in
'unified' compute mode (ANUGA_DEFAULT_COMPUTE_MODE=unified) on a GPU build.

Writes test5_gauges.csv (level, depth and speed at points 1-7 every
1 min), test5_volume.txt and test5.sww every hour.
"""

import os
import sys

import numpy as np
import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_common as ea  # noqa: E402

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

DEFAULT_AREA = 2500.0      # the 50 m grid
MAX_AREA = float(os.environ.get('EA_5_MAX_AREA', DEFAULT_AREA))
NAME = 'test5' if MAX_AREA == DEFAULT_AREA else 'test5_a%g' % MAX_AREA
FINAL_TIME = 30 * 3600.0
YIELD_STEP = 60.0
OUTPUT_STEP = 3600.0
# The inflow line (Test5BC_backup.txt)
INFLOW_LINE = [[232595.0, 830480.0], [232785.0, 830290.0]]

x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test5', 'Test5DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test5', 'Test5BC.csv'))
points = ea.read_points(ea.data_path('Test5', 'Test5output.csv'))
outline = ea.read_polygon(ea.data_path('Test5', 'Test5ActiveArea_region.shp'))

domain = anuga.create_domain_from_regions(
    outline,
    boundary_tags={'wall': list(range(len(outline)))},
    maximum_triangle_area=MAX_AREA,
    breaklines=[INFLOW_LINE],
    use_cache=False, verbose=False)

domain.set_name(NAME)
domain.set_flow_algorithm(alg)
domain.set_quantity('elevation',
                    ea.absolute(domain, ea.grid_function(x_dem, y_dem, z_dem)),
                    location='centroids')
domain.set_quantity('friction', 0.04, location='centroids')
domain.set_quantity('stage', expression='elevation', location='centroids')

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'wall': Br})

from anuga.validation_utilities import save_parameters_tex  # noqa: E402
save_parameters_tex(domain)

anuga.Inlet_operator(domain, INFLOW_LINE, Q=ea.inflow_function(t_bc, q_bc))

gauges = ea.Gauges(domain, points)
for t in domain.evolve(yieldstep=YIELD_STEP, outputstep=OUTPUT_STEP,
                       finaltime=FINAL_TIME):
    gauges.sample()
    if verbose and abs(t % 3600.0) < 1.0e-6:
        print(domain.timestepping_statistics())

gauges.save(NAME + '_gauges.csv')
with open(NAME + '_volume.txt', 'w') as f:
    f.write('inflow %.6f\n' % ea.inflow_volume(t_bc, q_bc, FINAL_TIME))
    f.write('final %.6f\n' % ea.domain_volume(domain))
    f.write('triangles %d\n' % domain.number_of_triangles)
print('compute mode %s' % domain.compute_mode)
