"""This is the public API to ANUGA_PARALLEL.

Ideally, all tools needed to run paralllel simulations should be 
imported from this module or anuga
"""


#FIXME (SR) Should put this into __init__.py

from anuga_parallel.parallel_api import distribute
from anuga_parallel.parallel_api import myid, numprocs, get_processor_name
from anuga_parallel.parallel_api import send, receive
from anuga_parallel.parallel_api import pypar_available, barrier, finalize

from anuga_parallel.parallel_meshes import parallel_rectangle

from parallel_shallow_water import Parallel_domain as Parallel_shallow_water_domain
from parallel_advection     import Parallel_domain as Parallel_advection_domain
