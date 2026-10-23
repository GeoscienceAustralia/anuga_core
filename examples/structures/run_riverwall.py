"""
Riverwall example
=================

Demonstrates how to:

  1. Build a mesh with breaklines so that triangle edges align with the walls.
  2. Create two named riverwalls with different crest heights and hydraulic
     parameters.
  3. Inspect and modify wall state *during* the evolve loop using the
     RiverWall runtime interface:
       - ``get_elevation`` / ``set_elevation``
       - ``set_elevation_offset``
       - ``get_hydraulic_parameter`` / ``set_hydraulic_parameter``

Domain layout (plan view, 60 m × 20 m)
---------------------------------------

  x=0      x=20        x=40        x=60
  |         |           |           |
  |  pool   | leveeA    |  channel  | leveeB  | outlet |
  |  h≈1.5  | z=1.0     |           | z=0.8   |        |
  |_________|___________|___________|_________|________|

Water starts at stage 1.5 m on the left, 0 m elsewhere.
leveeA (x=20) holds the water back initially; after 30 s we raise
leveeB (x=40) to 2 m to demonstrate mid-simulation elevation changes.

Run
---
    python run_riverwall.py

Output
------
    riverwall_example.sww   — viewable in ANUGA viewer or ParaView
"""

import numpy as np
import anuga

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

bounding_polygon = [
    [0.0,  0.0],
    [60.0, 0.0],
    [60.0, 20.0],
    [0.0,  20.0],
]

boundary_tags = {
    'bottom': [0],
    'right':  [1],
    'top':    [2],
    'left':   [3],
}

# Walls are defined as 3-D polylines: [[x, y, z], ...]
# x,y are plan-view coordinates; z is the crest elevation in metres.
river_walls = {
    'leveeA': [[20.0,  0.0, 1.0],
               [20.0, 20.0, 1.0]],
    'leveeB': [[40.0,  0.0, 0.8],
               [40.0, 20.0, 0.8]],
}

# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
# Pass the wall polylines as *breaklines* so that triangle edges are forced
# to align exactly with the wall — this is required for riverwalls to work.

anuga.create_pmesh_from_regions(
    bounding_polygon,
    boundary_tags=boundary_tags,
    maximum_triangle_area=4.0,
    breaklines=list(river_walls.values()),
    use_cache=False,
    verbose=False,
    filename='riverwall_example.msh',
)

domain = anuga.create_domain_from_file('riverwall_example.msh')
domain.set_flow_algorithm('DE0')   # riverwalls require a DE algorithm
domain.set_name('riverwall_example')

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

def initial_stage(x, y):
    """Pool on the left, dry elsewhere."""
    return np.where(x < 20.0, 1.5, 0.0)

domain.set_quantity('elevation', 0.0, location='centroids')
domain.set_quantity('friction',  0.03, location='centroids')
domain.set_quantity('stage', initial_stage, location='centroids')

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------

Br = anuga.Reflective_boundary(domain)
Bout = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])   # open outflow on the right

domain.set_boundary({'left': Br, 'top': Br, 'bottom': Br, 'right': Bout})

# ---------------------------------------------------------------------------
# Riverwalls
# ---------------------------------------------------------------------------
# ``riverwallPar`` sets hydraulic parameters per wall; any parameter not
# specified here uses the default values (Qfactor=1, s1=0.9, s2=0.95, …).

riverwall_par = {
    'leveeA': {'Qfactor': 1.0, 's1': 0.8, 's2': 0.95,
               'h1': 1.0, 'h2': 1.5, 'Cd_through': 0.0},
    'leveeB': {'Qfactor': 0.9, 'Cd_through': 0.0},
}

domain.riverwallData.create_riverwalls(
    river_walls,
    riverwallPar=riverwall_par,
    verbose=True,
)

rw = domain.riverwallData   # shorthand for the runtime interface

print('\nInitial crest elevations:')
for name in rw.get_wall_names():
    elev = rw.get_elevation(name)
    print(f'  {name}: min={elev.min():.3f}  max={elev.max():.3f}  '
          f'n_edges={len(elev)}')

print('\nHydraulic parameters (leveeA):')
for param in rw.hydraulic_variable_names:
    print(f'  {param} = {rw.get_hydraulic_parameter("leveeA", param):.4f}')

# ---------------------------------------------------------------------------
# Evolve
# ---------------------------------------------------------------------------

print('\nEvolving ...\n')

for t in domain.evolve(yieldstep=5.0, finaltime=60.0):
    domain.print_timestepping_statistics()

    # --- runtime modifications inside the evolve loop ---

    if abs(t - 30.0) < 1e-6:
        # Raise leveeB crest by 1.2 m at t=30 s to slow drainage
        print('\n  [t=30] Raising leveeB crest by +1.2 m')
        rw.set_elevation_offset('leveeB', 1.2)
        print(f'         New leveeB elevation: '
              f'{rw.get_elevation("leveeB").mean():.3f} m')

    if abs(t - 45.0) < 1e-6:
        # Lower leveeA to 0.5 m to release the pool
        print('\n  [t=45] Lowering leveeA crest to 0.5 m')
        rw.set_elevation('leveeA', 0.5)
        print(f'         New leveeA elevation: '
              f'{rw.get_elevation("leveeA").mean():.3f} m')

    if abs(t - 50.0) < 1e-6:
        # Increase leveeA Qfactor to amplify flow over the lowered crest
        print('\n  [t=50] Setting leveeA Qfactor to 1.5')
        rw.set_hydraulic_parameter('leveeA', 'Qfactor', 1.5)

print('\nDone.  Output written to riverwall_example.sww')
