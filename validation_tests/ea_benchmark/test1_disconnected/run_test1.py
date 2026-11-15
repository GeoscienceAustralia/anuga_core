"""
EA 2D benchmark Test 1: flooding a disconnected water body.

A 700 m x 100 m strip slopes up from x = 0 to a sill (10.25 m) at x = 300 m,
beyond which lies a depression (floor 9.75 m) rising to 10.39 m at the far
end. The water level along the left-hand edge rises from 9.7 m to 10.35 m
(5-60 min), holds until 11 h and falls back to 9.7 m by 12 h, so the
depression fills to 10.35 m and then drains to the sill. The other edges
are walls; the initial level is 9.7 m (only the first ~25 m are wet).
Manning n = 0.03, 10 m resolution, run to 20 h. The level is imposed with
a transmissive normal-momentum boundary, so water can leave as freely as
it enters.

Writes test1_gauges.csv (level, depth and speed at points 1 and 2 every
1 min), test1_volume.txt and test1.sww every hour.
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

MAX_AREA = 100.0           # the 10 m grid (ANUGA used 714 triangles in 2009)
FINAL_TIME = 20 * 3600.0
YIELD_STEP = 60.0
OUTPUT_STEP = 3600.0
INITIAL_LEVEL = 9.7

x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test1', 'test1DEM.asc'))
t_bc, h_bc = ea.read_inflow(ea.data_path('Test1', 'Test1BC.csv'))   # level, not flow
points = ea.read_points(ea.data_path('Test1', 'Test1output.csv'))

domain = anuga.create_domain_from_regions(
    [[0.0, 0.0], [700.0, 0.0], [700.0, 100.0], [0.0, 100.0]],
    boundary_tags={'bottom': [0], 'right': [1], 'top': [2], 'left': [3]},
    maximum_triangle_area=MAX_AREA,
    use_cache=False, verbose=False)

domain.set_name('test1')
domain.set_flow_algorithm(alg)
domain.set_quantity('elevation',
                    ea.absolute(domain, ea.grid_function(x_dem, y_dem, z_dem)),
                    location='centroids')
domain.set_quantity('friction', 0.03, location='centroids')
elev = domain.quantities['elevation'].centroid_values
domain.set_quantity('stage', np.maximum(elev, INITIAL_LEVEL), location='centroids')

Br = anuga.Reflective_boundary(domain)
Bt = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
    domain, function=lambda t: float(np.interp(t, t_bc, h_bc)))
domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

from anuga.validation_utilities import save_parameters_tex  # noqa: E402
save_parameters_tex(domain)

gauges = ea.Gauges(domain, points)
for t in domain.evolve(yieldstep=YIELD_STEP, outputstep=OUTPUT_STEP,
                       finaltime=FINAL_TIME):
    gauges.sample()
    if verbose and abs(t % 3600.0) < 1.0e-6:
        print(domain.timestepping_statistics())

gauges.save('test1_gauges.csv')
with open('test1_volume.txt', 'w') as f:
    f.write('final %.6f\n' % ea.domain_volume(domain))
    f.write('triangles %d\n' % domain.number_of_triangles)
print('compute mode %s' % domain.compute_mode)
