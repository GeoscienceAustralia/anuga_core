"""
EA 2D benchmark Test 3: momentum conservation over a small obstruction.

A 300 m x 100 m channel slopes down at 1:200 into two depressions, at
X = 150 m and X = 250 m, separated by a low crest at X = 200 m. A short
inflow (1,310 m^3) at the left end just fills the first depression; a
model that conserves momentum carries some water over the crest into the
second. Manning n = 0.01, 5 m resolution (~1,200 elements), 900 s.

Writes test3_gauges.csv (level and speed at points 1 and 2 every 2 s) and
test3_volume.txt (inflow and final volumes).
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_common as ea  # noqa: E402

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

L, W = 300.0, 100.0
MAX_AREA = 25.0            # 5 m grid: ~1,200 triangles over 30,000 m^2
FINAL_TIME = 900.0
OUTPUT_STEP = 2.0
POINTS = {'1': (150.0, 50.0), '2': (250.0, 50.0)}

x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test3', 'test3DEM.asc'))
t_bc, q_bc = ea.read_inflow(ea.data_path('Test3', 'Test3BC.csv'))

domain = anuga.create_domain_from_regions(
    [[0.0, 0.0], [L, 0.0], [L, W], [0.0, W]],
    boundary_tags={'bottom': [0], 'right': [1], 'top': [2], 'left': [3]},
    maximum_triangle_area=MAX_AREA,
    use_cache=False, verbose=False)

domain.set_name('test3')
domain.set_flow_algorithm(alg)
domain.set_quantity('elevation',
                    ea.absolute(domain, ea.grid_function(x_dem, y_dem, z_dem)),
                    location='centroids')
domain.set_quantity('friction', 0.01, location='centroids')
domain.set_quantity('stage', expression='elevation', location='centroids')

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

# Inflow across the full 100 m width, one cell in from the left wall.
inflow_line = [[2.5, 0.0], [2.5, W]]
anuga.Inlet_operator(domain, inflow_line, Q=ea.inflow_function(t_bc, q_bc))

gauges = ea.Gauges(domain, POINTS)
for t in domain.evolve(yieldstep=OUTPUT_STEP, finaltime=FINAL_TIME):
    gauges.sample()
    if verbose and abs(t % 60.0) < 1.0e-6:
        print(domain.timestepping_statistics())

gauges.save('test3_gauges.csv')
with open('test3_volume.txt', 'w') as f:
    f.write('inflow %.6f\n' % ea.inflow_volume(t_bc, q_bc, FINAL_TIME))
    f.write('final %.6f\n' % ea.domain_volume(domain))
    f.write('triangles %d\n' % domain.number_of_triangles)
