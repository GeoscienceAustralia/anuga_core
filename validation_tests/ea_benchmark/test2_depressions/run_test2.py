"""
EA 2D benchmark Test 2: filling of floodplain depressions.

A 2000 m x 2000 m square holds a 4 x 4 grid of identical ~0.5 m deep
depressions (a product of sinusoids) on an average slope of 1:1500 to the
south and 1:3000 to the east. A hydrograph rising to 20 m^3/s (5-10 min),
holding until 86 min and back to zero at 91 min (97,200 m^3 in all) enters
along the 100 m of the western edge below the north-western corner. All
boundaries are closed, the bed starts dry, Manning n = 0.03, 20 m
resolution, run to 48 h.

EA_2_MAX_AREA sets the triangle area (default 400 m^2, the 20 m grid); the
outputs are then named test2_a<area>_*.

Writes test2_gauges.csv (level, depth and speed at the centres of the 16
depressions every 5 min), test2_volume.txt and test2.sww every 2 h.
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_common as ea  # noqa: E402

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

DEFAULT_AREA = 400.0       # the 20 m grid
MAX_AREA = float(os.environ.get('EA_2_MAX_AREA', DEFAULT_AREA))
NAME = 'test2' if MAX_AREA == DEFAULT_AREA else 'test2_a%g' % MAX_AREA
FINAL_TIME = 48 * 3600.0
YIELD_STEP = 300.0
OUTPUT_STEP = 7200.0
# The inflow: 100 m of the western edge (Test2BC_polyline), as a line just
# inside the domain so that it picks the triangles along the edge.
INFLOW_LINE = [[1.0, 2000.0], [1.0, 1900.0]]

x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test2', 'test2DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test2', 'Test2_BC.csv'))
points = ea.read_points(ea.data_path('Test2', 'Test2output.csv'))

domain = anuga.create_domain_from_regions(
    [[0.0, 0.0], [2000.0, 0.0], [2000.0, 2000.0], [0.0, 2000.0]],
    boundary_tags={'wall': [0, 1, 2, 3]},
    maximum_triangle_area=MAX_AREA,
    use_cache=False, verbose=False)

domain.set_name(NAME)
domain.set_flow_algorithm(alg)
domain.set_quantity('elevation',
                    ea.absolute(domain, ea.grid_function(x_dem, y_dem, z_dem)),
                    location='centroids')
domain.set_quantity('friction', 0.03, location='centroids')
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
    if verbose and abs(t % 7200.0) < 1.0e-6:
        print(domain.timestepping_statistics())

gauges.save(NAME + '_gauges.csv')
with open(NAME + '_volume.txt', 'w') as f:
    f.write('inflow %.6f\n' % ea.inflow_volume(t_bc, q_bc, FINAL_TIME))
    f.write('final %.6f\n' % ea.domain_volume(domain))
    f.write('triangles %d\n' % domain.number_of_triangles)
print('compute mode %s' % domain.compute_mode)
