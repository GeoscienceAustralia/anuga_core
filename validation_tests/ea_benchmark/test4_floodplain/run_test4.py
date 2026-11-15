"""
EA 2D benchmark Test 4: speed of flood propagation over an extended floodplain.

A flat, horizontal 1000 m x 2000 m floodplain (bed at 0 m) is flooded
through a 20 m breach in the middle of its western side, with a peak
inflow of 20 m^3/s. Manning n = 0.05, 5 m resolution, run to 5 h.

Writes test4_gauges.csv (level, depth and speed at the six output points
every 20 s), test4_volume.txt, and test4.sww every 30 min.
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_common as ea  # noqa: E402

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

L, W = 1000.0, 2000.0
MAX_AREA = 25.0            # 5 m grid
FINAL_TIME = 5 * 3600.0
YIELD_STEP = 20.0
OUTPUT_STEP = 1800.0

t_bc, q_bc = ea.read_inflow(ea.data_path('Test4', 'Test4BC.csv'))
points = ea.read_points(ea.data_path('Test4', 'Test4output.csv'))

domain = anuga.create_domain_from_regions(
    [[0.0, 0.0], [L, 0.0], [L, W], [0.0, W]],
    boundary_tags={'bottom': [0], 'right': [1], 'top': [2], 'left': [3]},
    maximum_triangle_area=MAX_AREA,
    use_cache=False, verbose=False)

domain.set_name('test4')
domain.set_flow_algorithm(alg)
domain.set_quantity('elevation', 0.0, location='centroids')
domain.set_quantity('friction', 0.05, location='centroids')
domain.set_quantity('stage', 0.0, location='centroids')

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

# The 20 m breach, in the middle of the western side, one cell in.
inflow_line = [[2.5, W / 2 - 10.0], [2.5, W / 2 + 10.0]]
anuga.Inlet_operator(domain, inflow_line, Q=ea.inflow_function(t_bc, q_bc))

gauges = ea.Gauges(domain, points)
for t in domain.evolve(yieldstep=YIELD_STEP, outputstep=OUTPUT_STEP,
                       finaltime=FINAL_TIME):
    gauges.sample()
    if verbose and abs(t % 600.0) < 1.0e-6:
        print(domain.timestepping_statistics())

gauges.save('test4_gauges.csv')
with open('test4_volume.txt', 'w') as f:
    f.write('inflow %.6f\n' % ea.inflow_volume(t_bc, q_bc, FINAL_TIME))
    f.write('final %.6f\n' % ea.domain_volume(domain))
    f.write('triangles %d\n' % domain.number_of_triangles)
