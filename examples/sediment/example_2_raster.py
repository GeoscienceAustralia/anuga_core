"""Example 2 -- sediment transport over a real DEM.

Ported from anugaSed's `run_raster_sed_transport.py` (Perignon 2016, MIT), using
their `topo.asc` and `outline.csv` unchanged. Same mesh resolution, boundary
tags, inflow stage and duration; our API and the non-cohesive erosion route --
see README.md.

Demonstrates the module on an unstructured mesh built from a raster, rather
than a rectangular test domain.
"""
import os
import numpy as np
import anuga

FILENAME_ROOT = 'topo'
INFLOW_CONC = 0.01
D50 = 65.0e-6
FINALTIME = 30.0
INFLOW_STAGE = 1528.0        # their Dirichlet_boundary([1528, 0, 0])

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

# --- their mesh, from the raster -------------------------------------------
if not os.path.exists(FILENAME_ROOT + '.pts'):
    anuga.asc2dem(FILENAME_ROOT + '.asc', use_cache=False, verbose=False)
    anuga.dem2pts(FILENAME_ROOT + '.dem', use_cache=False, verbose=False)

bounding_polygon = anuga.read_polygon('outline.csv')
boundary_tags = {'bottom': [0], 'side1': [1], 'side2': [2],
                 'top': [3], 'side3': [4], 'side4': [5]}

from anuga.pmesh.mesh_interface import create_mesh_from_regions
create_mesh_from_regions(bounding_polygon=bounding_polygon,
                         boundary_tags=boundary_tags,
                         maximum_triangle_area=200,
                         filename=FILENAME_ROOT + '.msh')

domain = anuga.Domain(FILENAME_ROOT + '.msh')
domain.set_name('sed_raster')
domain.set_flow_algorithm('DE0')
print('Triangles: %d      extent: %s' % (len(domain), domain.get_extent()))

domain.set_quantity('elevation', filename=FILENAME_ROOT + '.pts',
                    use_cache=False, verbose=False, alpha=0.1)
domain.set_quantity('stage', expression='elevation')     # dry start
domain.set_quantity('friction', 0.03)

min_elev = domain.quantities['elevation'].vertex_values.min()
Bd = anuga.Dirichlet_boundary([INFLOW_STAGE, 0.0, 0.0])
Bi = anuga.Dirichlet_boundary([min_elev - 1.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)
# Their original also mapped an 'exterior' tag; this mesh has no such tag, so
# it is dropped rather than carried over.
domain.set_boundary({'bottom': Bi, 'side1': Br, 'side2': Br,
                     'top': Bd, 'side3': Br, 'side4': Br})

# --- our API ----------------------------------------------------------------
domain.sediment_porosity = 0.3
domain.add_sediment_fraction(name='sand', diameter=D50, rho_s=2650.0, initial_concentration=0.0)
domain.set_tracer_boundary('sand', 'top', INFLOW_CONC)    # 'top' carries INFLOW_STAGE

areas = domain.areas
z0 = domain.quantities['elevation'].centroid_values.copy()
lam = domain.sediment_porosity

print('v_s = %.3e m/s   D50 = %.0f um   inflow c = %.3f   inflow stage = %.0f m'
      % (domain.sediment_settling_velocity[0], D50 * 1e6, INFLOW_CONC,
         INFLOW_STAGE))
print()
print('   time    wet cells   water vol      suspended      bed change')
for t in domain.evolve(yieldstep=10.0, finaltime=FINALTIME):
    z = domain.quantities['elevation'].centroid_values
    h = np.maximum(domain.quantities['stage'].centroid_values - z, 0.0)
    m = float((domain.tracer_conserved_values[0] * areas).sum())
    bed = float(((1.0 - lam) * (z - z0) * areas).sum())
    print('  %5.1f s   %7d   %10.1f m3  %.6e   %+.4e'
          % (t, int((h > 0.01).sum()), float((h * areas).sum()), m, bed))

dz = domain.quantities['elevation'].centroid_values - z0
print()
print('Bed change: min %+.4e m   max %+.4e m   over %d triangles'
      % (dz.min(), dz.max(), dz.size))
print('Output written to sed_raster.sww')
