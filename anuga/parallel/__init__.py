"""
This is the public API to ANUGA_PARALLEL.

Ideally, all tools needed to run parallel simulations should be
imported from this module
"""



from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester


from .parallel_api import distribute
from .parallel_api import myid, numprocs, get_processor_name
from .parallel_api import send, receive
from .parallel_api import pypar_available, barrier, finalize

if pypar_available:
    from .parallel_meshes import parallel_rectangle
    from .parallel_shallow_water import Parallel_domain as Parallel_shallow_water_domain
    from .parallel_advection     import Parallel_domain as Parallel_advection_domain
    from .parallel_operator_factory import Inlet_operator, Boyd_box_operator, Boyd_pipe_operator
    from .parallel_operator_factory import Weir_orifice_trapezoid_operator
else:
    from anuga import rectangular_cross as parallel_rectangle
    from anuga import Domain as Parallel_shallow_water_domain
    from anuga.advection import Advection_Domain as Parallel_advection_domain



