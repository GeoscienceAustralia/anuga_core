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



from ._version import __git_sha__
from ._version import __git_committed_datetime__
from ._version import __version__

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
from anuga.coordinate_transforms.redfearn import epsg_to_ll, ll_to_epsg

from anuga.abstract_2d_finite_volumes.pmesh2domain import \
                                            pmesh_to_domain_instance, \
                                            pmesh_to_mesh, \
                                            pmesh_to_basic_mesh

from anuga.abstract_2d_finite_volumes.basic_mesh import (
    Basic_mesh, rectangular_basic_mesh, rectangular_cross_basic_mesh,
    basic_mesh_from_mesh_file)

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
from anuga.utilities.parse_time import seconds_to_hhmmss

# ----------------------------
# Parallel api
# ----------------------------
from anuga.parallel.parallel_api import distribute
from anuga.parallel.parallel_api import distribute_collaborative
from anuga.parallel.parallel_api import distribute_basic_mesh
from anuga.parallel.parallel_api import distribute_basic_mesh_collaborative
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
from anuga.pmesh.mesh_interface import create_pmesh_from_regions
from anuga.pmesh.mesh_interface import create_mesh_from_regions  # deprecated alias

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
from anuga.utilities.system_tools import memory_stats
from anuga.utilities.system_tools import print_memory_stats
from anuga.utilities.mem_time_equation import estimate_time_mem

# -------------------------
# create domain functions
# -------------------------
from anuga.extras import create_domain_from_regions
from anuga.extras import create_domain_from_file
from anuga.extras import rectangular_cross_domain
from anuga.extras import create_basic_mesh_from_regions

from anuga.utilities import log as log
from anuga.utilities.log import set_logfile, TeeStream, file_only
from anuga.utilities.log import verbose as log_verbose

from anuga.config import g
from anuga.config import velocity_protection
from anuga.config import MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU
from anuga.config import LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2

# --------------------------------------
# Public API — names exported by `import anuga`
# --------------------------------------
__all__ = [
    # Core classes
    'Basic_mesh',
    'basic_mesh_from_mesh_file',
    'Domain',
    'Domain_plotter',
    'Generic_Domain',
    'Geo_reference',
    'Geospatial_data',
    'Mesh',
    'Operator',
    'Quantity',
    'Region',
    'Simulation',
    'Structure_operator',
    'SWW_plotter',
    # Boundaries
    'Characteristic_stage_boundary',
    'Compute_fluxes_boundary',
    'Dirichlet_boundary',
    'Field_boundary',
    'File_boundary',
    'Flather_external_stage_zero_velocity_boundary',
    'Reflective_boundary',
    'Time_boundary',
    'Time_space_boundary',
    'Time_stage_zero_momentum_boundary',
    'Transmissive_boundary',
    'Transmissive_momentum_set_stage_boundary',
    'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary',
    'Transmissive_stage_zero_momentum_boundary',
    # Operators
    'Bed_shear_erosion_operator',
    'Flat_fill_slice_erosion_operator',
    'Flat_slice_erosion_operator',
    'Kinematic_viscosity_operator',
    'Rate_operator',
    'Sanddune_erosion_operator',
    'Set_depth_friction_operator',
    'Set_elevation',
    'Set_elevation_operator',
    'Set_quantity',
    'Set_quantity_operator',
    'Set_stage',
    'Set_stage_operator',
    # Structure operators (parallel or serial depending on pypar)
    'Boyd_box_operator',
    'Boyd_pipe_operator',
    'Inlet_operator',
    'Internal_boundary_operator',
    'Weir_orifice_trapezoid_operator',
    # Forcing (legacy)
    'Inflow',
    'Rainfall',
    'Wind_stress',
    # Geometry and polygon utilities
    'inside_polygon',
    'plot_polygons',
    'Polygon_function',
    'polygon_area',
    'read_polygon',
    # Mesh factory helpers
    'rectangular',
    'rectangular_basic_mesh',
    'rectangular_cross',
    'rectangular_cross_basic_mesh',
    'rectangular_cross_domain',
    # pmesh helpers
    'pmesh_to_basic_mesh',
    'pmesh_to_domain_instance',
    'pmesh_to_mesh',
    # File and format conversion
    'asc2dem',
    'create_sts_boundary',
    'dem2dem',
    'dem2pts',
    'esri2sww',
    'ferret2sww',
    'llasc2pts',
    'load_checkpoint_file',
    'load_csv_as_building_polygons',
    'load_csv_as_polygons',
    'load_ungenerate',
    'sww2array',
    'sww2csv_gauges',
    'sww2dem',
    'sww2dem_batch',
    'sww2obj',
    'sww2timeseries',
    'sww_merge',
    'timefile2netcdf',
    'tsh2sww',
    'urs2nc',
    'urs2sts',
    'urs2sww',
    'xya2pts',
    # Fitting and interpolation
    'fit_to_mesh',
    'fit_to_mesh_file',
    'file_function',
    'csv2timeseries_graphs',
    # Coordinate transforms
    'LLtoUTM',
    'UTMtoLL',
    'epsg_to_ll',
    'll_to_epsg',
    # Parallel API
    'barrier',
    'collect_value',
    'distribute',
    'distribute_basic_mesh',
    'distribute_basic_mesh_collaborative',
    'distribute_collaborative',
    'finalize',
    'mpi_extra_options',
    'mpicmd',
    'myid',
    'numprocs',
    'get_processor_name',
    'pypar_available',
    'receive',
    'reduce',
    'send',
    'sequential_distribute_dump',
    'sequential_distribute_load',
    # Model tools / polygon utilities
    'Create_culvert_bridge_Operator',
    'get_polygon_dictionary',
    'get_polygon_from_single_file',
    'get_polygon_list_from_files',
    'get_polygon_value_list',
    'get_polygons_from_Mid_Mif',
    'get_WCC_2002_Blockage_factor',
    'get_WCC_2016_Blockage_factor',
    'pumping_station_function',
    'read_hole_dir_multi_files_with_single_poly',
    'read_hole_dir_single_file_with_multi_poly',
    'read_multi_poly_file',
    'read_multi_poly_file_value',
    'read_polygon_dir',
    # Domain creation
    'create_basic_mesh_from_regions',
    'create_domain_from_file',
    'create_domain_from_regions',
    'create_mesh_from_regions',   # deprecated alias for create_pmesh_from_regions
    'create_pmesh_from_regions',
    # System / utility
    'acos',
    'cache',
    'copy_code_files',
    'estimate_time_mem',
    'file_length',
    'get_args',
    'get_host_name',
    'get_revision_date',
    'get_revision_number',
    'get_user_name',
    'get_version',
    'memory_stats',
    'print_memory_stats',
    'log',
    'set_logfile',
    'TeeStream',
    'parse_standard_args',
    'parse_time',
    'plot_utils',
    'run_anuga_script',
    'seconds_to_hhmmss',
    # Tsunami source
    'slide_tsunami',
    'slump_tsunami',
    # Interrogation
    'get_flow_through_cross_section',
    # Config constants
    'g',
    'indent',
    'LOW_FROUDE_1',
    'LOW_FROUDE_2',
    'LOW_FROUDE_OFF',
    'MULTIPROCESSOR_GPU',
    'MULTIPROCESSOR_OPENMP',
    'velocity_protection',
    # Parsing
    'create_standard_parser',
    # Testing
    'test',
]

# --------------------------------------
# NetCDF changes stdout to the terminal
# This resets it
# --------------------------------------
try:
    from importlib import reload
except ImportError:
    pass
reload(sys)
sys.stdout = _stdout
