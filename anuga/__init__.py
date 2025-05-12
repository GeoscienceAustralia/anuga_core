""" ANUGA models the effect of tsunamis and flooding upon a terrain mesh.
    In typical usage, a Domain class is created for a particular piece of
    terrain. Boundary conditions are specified for the domain, such as inflow
    and outflow, and then the simulation is run.

    This is the public API to ANUGA. It provides a toolkit of often-used
    modules, which can be used directly by including the following line in
    the user's code:

    >>> import anuga

    This usage pattern abstracts away the internal heirarchy of the ANUGA
    system, allowing the user to concentrate on writing simulations without
    searching through the ANUGA source tree for the functions that they need.

    Also, it isolates the user from "under-the-hood" refactorings.
"""

# -----------------------------------------------------
# Make selected classes available directly
# -----------------------------------------------------



from .revision import  __git_sha__
from .revision import __git_committed_datetime__
from .revision import __version__

# ----------------------------------
# NetCDF changes stdout to terminal
# Causes trouble when using jupyter
# ---------------------------------
import sys
_stdout = sys.stdout

# ---------------------------------
# Setup the tester from numpy
# ---------------------------------
from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester

#from anuga.__config__ import show as show_config

# --------------------------------
# Important basic classes
# --------------------------------
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.abstract_2d_finite_volumes.region import Region
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.operators.base_operator import Operator
from anuga.structures.structure_operator import Structure_operator

from anuga.utilities.animate import SWW_plotter
from anuga.utilities.animate import Domain_plotter


from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

# ------------------------------------------------------------------------------
# Miscellaneous
# ------------------------------------------------------------------------------
from anuga.abstract_2d_finite_volumes.util import file_function, \
                                        sww2timeseries, sww2csv_gauges, \
                                        csv2timeseries_graphs

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                                    rectangular

from anuga.file.csv_file import load_csv_as_building_polygons,  \
                                load_csv_as_polygons

from anuga.file.sts import create_sts_boundary

from anuga.file.ungenerate import load_ungenerate

from anuga.geometry.polygon import read_polygon
from anuga.geometry.polygon import plot_polygons
from anuga.geometry.polygon import inside_polygon
from anuga.geometry.polygon import polygon_area
from anuga.geometry.polygon_function import Polygon_function

from anuga.coordinate_transforms.lat_long_UTM_conversion import LLtoUTM, UTMtoLL

from anuga.abstract_2d_finite_volumes.pmesh2domain import \
                                            pmesh_to_domain_instance

from anuga.fit_interpolate.fit import fit_to_mesh_file
from anuga.fit_interpolate.fit import fit_to_mesh

from anuga.utilities.system_tools import file_length
from anuga.utilities.sww_merge import sww_merge_parallel as sww_merge
from anuga.utilities.file_utils import copy_code_files
from anuga.utilities.numerical_tools import safe_acos as acos
import anuga.utilities.plot_utils as plot_utils


from anuga.caching import cache
from os.path import join
from anuga.config import indent

from anuga.utilities.parse_time import parse_time

# ----------------------------
# Parallel api
# ----------------------------
from anuga.parallel.parallel_api import distribute
from anuga.parallel.parallel_api import myid, numprocs, get_processor_name
from anuga.parallel.parallel_api import send, receive, reduce
from anuga.parallel.parallel_api import pypar_available, barrier, finalize
from anuga.parallel.parallel_api import collect_value
from anuga.parallel.parallel_api import mpicmd
from anuga.parallel.parallel_api import mpi_extra_options


from anuga.parallel.parallel_api import sequential_distribute_dump
from anuga.parallel.parallel_api import sequential_distribute_load
        

# -----------------------------
# Checkpointing
# -----------------------------
from anuga.shallow_water.checkpoint import load_checkpoint_file

# -----------------------------
# SwW Standard Boundaries
# -----------------------------
from anuga.shallow_water.boundaries import File_boundary
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.shallow_water.boundaries import Characteristic_stage_boundary
from anuga.shallow_water.boundaries import Field_boundary
from anuga.shallow_water.boundaries import \
                    Time_stage_zero_momentum_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_stage_zero_momentum_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_momentum_set_stage_boundary
from anuga.shallow_water.boundaries import \
                    Transmissive_n_momentum_zero_t_momentum_set_stage_boundary
from anuga.shallow_water.boundaries import \
                    Flather_external_stage_zero_velocity_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
                    Compute_fluxes_boundary

# -----------------------------
# General Boundaries
# -----------------------------
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Dirichlet_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Time_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Time_space_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
                            import Transmissive_boundary

# -----------------------------
# Shallow Water Tsunamis
# -----------------------------
from anuga.tsunami_source.smf import slide_tsunami, slump_tsunami

# -----------------------------
# Forcing
# These are old, should use operators
# -----------------------------
from anuga.shallow_water.forcing import Inflow, Rainfall, Wind_stress

# -----------------------------
# File conversion utilities
# -----------------------------
from anuga.file_conversion.file_conversion import sww2obj
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.file_conversion.file_conversion import tsh2sww
from anuga.file_conversion.urs2nc import urs2nc
from anuga.file_conversion.urs2sww import urs2sww
from anuga.file_conversion.urs2sts import urs2sts
from anuga.file_conversion.dem2pts import dem2pts
from anuga.file_conversion.esri2sww import esri2sww
from anuga.file_conversion.sww2dem import sww2dem, sww2dem_batch
from anuga.file_conversion.asc2dem import asc2dem
from anuga.file_conversion.xya2pts import xya2pts
from anuga.file_conversion.ferret2sww import ferret2sww
from anuga.file_conversion.dem2dem import dem2dem
from anuga.file_conversion.sww2array import sww2array
from anuga.file_conversion.llasc2pts import llasc2pts

# -----------------------------
# Parsing arguments
# -----------------------------
from anuga.utilities.argparsing import create_standard_parser
from anuga.utilities.argparsing import parse_standard_args


def get_args():
    """ Explicitly parse the argument list using standard anuga arguments

    Don't use this if you want to setup your own parser
    """
    parser = create_standard_parser()
    return parser.parse_args()

# -----------------------------
# Running Script
# -----------------------------
from anuga.utilities.run_anuga_script import run_script as run_anuga_script

# ---------------------------
# Simulation and Excel mesh_interface
# ---------------------------
from anuga.simulation.simulation import Simulation

# -----------------------------
# Mesh API
# -----------------------------
from anuga.pmesh.mesh_interface import create_mesh_from_regions

# -----------------------------
# SWW file access
# -----------------------------
from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section

# ---------------------------
# Operators
# ---------------------------
from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator

from anuga.operators.rate_operators import Rate_operator
from anuga.operators.set_friction_operators import Set_depth_friction_operator

from anuga.operators.set_elevation_operator import Set_elevation_operator
from anuga.operators.set_quantity_operator import Set_quantity_operator
from anuga.operators.set_stage_operator import Set_stage_operator


from anuga.operators.set_elevation import Set_elevation
from anuga.operators.set_quantity import Set_quantity
from anuga.operators.set_stage import Set_stage

from anuga.operators.sanddune_erosion_operator import Sanddune_erosion_operator
from anuga.operators.erosion_operators import Bed_shear_erosion_operator
from anuga.operators.erosion_operators import Flat_slice_erosion_operator
from anuga.operators.erosion_operators import Flat_fill_slice_erosion_operator

# ---------------------------
# Structure Operators
# ---------------------------



if pypar_available:
    from anuga.parallel.parallel_operator_factory import Inlet_operator
    from anuga.parallel.parallel_operator_factory import Boyd_box_operator
    from anuga.parallel.parallel_operator_factory import Boyd_pipe_operator
    from anuga.parallel.parallel_operator_factory import Weir_orifice_trapezoid_operator
    from anuga.parallel.parallel_operator_factory import Internal_boundary_operator
else:
    from anuga.structures.inlet_operator import Inlet_operator
    from anuga.structures.boyd_box_operator import Boyd_box_operator
    from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
    from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator
    from anuga.structures.internal_boundary_operator import Internal_boundary_operator

from anuga.structures.internal_boundary_functions import pumping_station_function



# ----------------------------
# Parallel distribute
# ----------------------------

# ----------------------------
#
# Added by Petar Milevski 10/09/2013

from anuga.utilities.model_tools import get_polygon_from_single_file
from anuga.utilities.model_tools import get_polygons_from_Mid_Mif
from anuga.utilities.model_tools import get_polygon_list_from_files
from anuga.utilities.model_tools import get_polygon_dictionary
from anuga.utilities.model_tools import get_polygon_value_list
from anuga.utilities.model_tools import read_polygon_dir
from anuga.utilities.model_tools import read_hole_dir_multi_files_with_single_poly
from anuga.utilities.model_tools import read_multi_poly_file
from anuga.utilities.model_tools import read_hole_dir_single_file_with_multi_poly
from anuga.utilities.model_tools import read_multi_poly_file_value
from anuga.utilities.model_tools import Create_culvert_bridge_Operator
from anuga.utilities.model_tools import get_WCC_2002_Blockage_factor
from anuga.utilities.model_tools import get_WCC_2016_Blockage_factor

# ---------------------------
# User Access Functions
# ---------------------------

from anuga.utilities.system_tools import get_user_name
from anuga.utilities.system_tools import get_host_name
from anuga.utilities.system_tools import get_version

from anuga.utilities.system_tools import get_revision_number
from anuga.utilities.system_tools import get_revision_date
from anuga.utilities.mem_time_equation import estimate_time_mem

# -------------------------
# create domain functions
# -------------------------
from anuga.extras import create_domain_from_regions
from anuga.extras import create_domain_from_file
from anuga.extras import rectangular_cross_domain

from anuga.utilities import log as log

from anuga.config import g
from anuga.config import velocity_protection

# --------------------------------------
# NetCDF changes stdout to the terminal
# This resets it
# --------------------------------------
try:
    from importlib import reload
except:
    pass
reload(sys)
sys.stdout = _stdout
