"""Run riverwall simulation (sequentially or in parallel) to support test_parallel_riverwall.py
"""

# ------------------------
# Import necessary modules
# ------------------------
from math import exp
from anuga import create_mesh_from_regions, create_domain_from_file
from anuga import Reflective_boundary, Transmissive_momentum_set_stage_boundary
from anuga import myid, distribute, barrier, numprocs, finalize

#----------------------------------
# set up MPI to abort on error
#----------------------------------
from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

verbose = False

alg = 'DE0'
scale_me = 1.0
mesh_filename = 'riverwall.msh'

# -----------
# Set up mesh
# -----------
boundaryPolygon = [[0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
midResPolygon = [[30., 30.], [30., 70.], [70., 70.], [70., 30.]]
higherResPolygon = [[40., 40.], [40., 60.], [60., 60.], [60., 40.]]

# Riverwall = list of lists, each with a set of x,y,z
# (and optional QFactor) values
riverWall = {'centralWall':
             [[50., 0.0, -0.0],
              [50., 45., -0.0],
              [50., 46., -0.2],
              [50., 100.0, -0.0]]
             }

riverWall_Par = {'centralWall': {'Qfactor': 1.0}}

# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
regionPtAreas = [[99., 99., 10.0*10.0*0.5],
                 [1., 1., 10.0*10.0*0.5],
                 [45., 45., 1.0*1.0*0.5],
                 [55., 55., 1.0*1.0*0.5],
                 [65., 65., 3.0*3.0*0.5],
                 [35., 35., 3.0*3.0*0.5]]

# --------------------------------------------------------------------------
# Setup computational domain and quantities
# --------------------------------------------------------------------------
if myid == 0:
    create_mesh_from_regions(boundaryPolygon,
                             boundary_tags={'left': [0],
                                            'top': [1],
                                            'right': [2],
                                            'bottom': [3]},
                             maximum_triangle_area=1.0e+20,
                             minimum_triangle_angle=28.0,
                             filename=mesh_filename,
                             interior_regions=[[higherResPolygon, 1.*1.*0.5],
                                               [midResPolygon, 3.0*3.0*0.5]],
                             breaklines=list(riverWall.values()),
                             use_cache=False,
                             verbose=verbose,
                             regionPtArea=regionPtAreas)

    base_domain = create_domain_from_file(mesh_filename)
    base_domain.set_flow_algorithm(alg)
    base_domain.set_datadir('.')
    base_domain.set_store_vertices_uniquely()

    # ----------------------------------------
    # Define topography and initial conditions
    # ----------------------------------------
    def topography(x, y):
        return -x/150. * scale_me

    def stagefun(x, y):
        stg = -0.5 * scale_me
        return stg

    base_domain.set_quantity('elevation', topography)
    base_domain.set_quantity('friction', 0.03)
    base_domain.set_quantity('stage', stagefun)
else:
    base_domain = None


# ----------------------------------------------
# Decide if this is a sequential or parallel run
# ----------------------------------------------
if numprocs == 1:
    # This is a sequential run
    domain = base_domain
    domain.set_name('s_riverwall')
else:
    # This is a parallel run
    domain = distribute(base_domain, verbose=verbose)
    domain.set_name('p_riverwall')
domain.set_store_vertices_uniquely()
domain.riverwallData.create_riverwalls(riverWall,
                                       riverWall_Par,
                                       verbose=verbose)

# -------------------------
# Setup boundary conditions
# -------------------------
Br = Reflective_boundary(domain)  # Solid reflective wall


def boundaryFun(t):
    output = -0.4 * exp(-t/100.) - 0.1
    output = min(output, -0.11)
    return output


Bin_tmss = Transmissive_momentum_set_stage_boundary(domain=domain,
                                                    function=boundaryFun)

# ---------------------------------------------
# Associate boundary tags with boundary objects
# ---------------------------------------------
domain.set_boundary({'left': Br, 'right': Bin_tmss, 'top': Br, 'bottom': Br})

# ------------------------------
# Evolve the system through time
# ------------------------------
for t in domain.evolve(yieldstep=10.0, finaltime=150.0):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())

# Wrap up parallel matters if required.
domain.sww_merge(delete_old=True)
finalize()
