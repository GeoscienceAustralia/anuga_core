"""Module where global ANUGA model parameters and default values are set
"""

import os
import sys


################################################################################
# Numerical constants
################################################################################

epsilon = 1.0e-12                   # Smallest number - used for safe division
max_float = 1.0e36                  # Largest number - used to initialise
                                    # (max, min) ranges
default_smoothing_parameter = 0.001 # Default alpha for penalised
                                    # least squares fitting
single_precision = 1.0e-6           # Smallest single precision number
velocity_protection = 1.0e-6        # Used to compute velocity from momentum
                                    # See section 7.4 on Flux limiting
                                    # in the user manual


################################################################################
# Standard filenames, directories and system parameters used by ANUGA
################################################################################

pmesh_filename = '.\\pmesh'
revision_filename = 'revision.py'
default_datadir = '.'
time_format = '%d/%m/%y %H:%M:%S'    # Used with timefile2netcdf
umask = 0o02  # Controls file and directory permission created by anuga (UNIX)
default_boundary_tag = 'exterior'
institution = 'Geosciences Australia'

# Major revision number for use with create_distribution
# and update_anuga_user_guide
#major_revision = '1.3.1'

################################################################################
# Physical constants
################################################################################

manning = 0.03  # Manning's friction coefficient
#g = 9.80665    # Gravity - FIXME reinstate this and fix unit tests.
g = 9.8
#g(phi) = 9780313 * (1 + 0.0053024 sin(phi)**2 - 0.000 0059 sin(2*phi)**2)
# micro m/s**2, where phi is the latitude
#The 'official' average is 9.80665

eta_w = 3.0e-3 # Wind stress coefficient
rho_a = 1.2e-3 # Atmospheric density
rho_w = 1023   # Fluid density [kg/m^3] (rho_w = 1023 for salt water)

################################################################################
# Limiters - used with linear reconstruction of vertex
# values from centroid values
################################################################################
# Note the individual beta values are set in domain.set_flow_method which also sets
# the timestepping method

beta_w = 1.0

# Alpha_balance controls how limiters are balanced between deep and shallow.
# A large value will favour the deep water limiters, allowing the a closer hug
# to the coastline.  This will minimise 'creep' but at the same time cause
# smaller time steps
# Range:
alpha_balance = 2.0

# Flag use of new limiters.
# tight_slope_limiters = 0 means use old limiters (e.g. for some tests)
# tight_slope_limiters = 1 means use new limiters that hug the bathymetry closer
tight_slope_limiters = True

use_edge_limiter = False    # The edge limiter is better, but most runs have been
                            # using vertex limiting. Validations passed with this
                            # one True 9th May 2008, but many unit tests need
                            # backward compatibility flag set FIXME(Ole).

# Use centroid velocities to reconstruct momentum at vertices in
# very shallow water
# This option has a first order flavour to it, but we still have second order
# reconstruction of stage and this option only applies in
# balance_deep_and_shallow when
# alpha < 1 so in deeper water the full second order scheme is used.
#
# This option is good with tight_slope_limiters, especially for large domains.
use_centroid_velocities = True

# FIXME (Ole) Maybe get rid of order altogether and use beta_w
default_order = 2

# Option to use velocity extrapolation instead of momentum extrapolation in the
# routine domain.extrapolate_second_order_sw
extrapolate_velocity_second_order=True

# Option to setup compute_fluxes_method
# Currently "original' and 'wb_1' to 'wb_3' and 'tsunami'
compute_fluxes_method = 'wb_2'

# Option to setup distribute_to_vertices_and_edges_method
# Currently "original' and 'tsunami'
distribute_to_vertices_and_edges_method = 'original'

# Option to turn on low damping for low Froude flows
low_froude = 0 

################################################################################
# Friction Method
################################################################################

sloped_mannings_function = False

################################################################################
# Timestepping
################################################################################

CFL = 1.0  # CFL condition assigned to domain.CFL - controls timestep size

# Choose type of timestepping and spatial reconstruction method

timestepping_method = 1

# For shallow water we have a method that sets both timestepping and spatial reconstruction and
# beta values. In this case the settings for timestepping_method will be overriden

#flow_algorithm = '1_0'    # 1st order euler and conservative piecewise constant spatial reconstruction
#flow_algorithm = '1_5'  # 1st order euler and conservative piecewise linear spatial reconstruction
#flow_algorithm = '1_75' # 1st order euler and more aggressive piecewise linear spatial reconstruction
#flow_algorithm = '2_0'    # 2nd order TVD scheme and more aggressive piecewise linear spatial reconstruction
#flow_algorithm = '2.5'  # 3rd order TVD scheme and more aggressive piecewise linear spatial reconstruction
#flow_algorithm = 'tsunami' # 2nd order space and time, well balanced inc at wet-dry fronts, porosity-type alg
flow_algorithm = 'DE0' # 1st order time 2nd order space, discontinuous elevation, well balanced + better shallow flows than 'tsunami'
#flow_algorithm = 'DE1' # 2nd order space and time, discontinuous elevation, well balanced + better shallow flows than 'tsunami'



# rk2 is a little more stable than euler, so rk2 timestepping
# can deal with a larger beta when slope limiting the reconstructed
# solution. The large beta is needed if solving problems sensitive
# to numerical diffusion, like a small forced wave in an ocean
beta_euler = 1.0
beta_rk2   = 1.6

# Option to search for signatures where isolated triangles are
# responsible for a small global timestep.
# Treating these by limiting their momenta may help speed up the
# overall computation.
# This facility is experimental.
# protect_against_isolated_degenerate_timesteps = False
protect_against_isolated_degenerate_timesteps = False

min_timestep = 1.0e-6 # Minimal timestep accepted in ANUGA
max_timestep = 1.0e+3
max_smallsteps = 50   # Max number of degenerate steps allowed b4
                      # trying first order

# Perhaps minimal timestep could be based on the geometry as follows:
# Define maximal possible speed in open water v_max, e.g. 500m/s (soundspeed?)
# Then work out minimal internal distance in mesh r_min and set
# min_timestep = r_min/v_max
#
# Max speeds are calculated in the flux function as
#
# lambda = v +/- sqrt(gh)
#
# so with 500 m/s, h ~ 500^2/g = 2500 m well out of the domain of the
# shallow water wave equation
#
# The actual soundspeed can be as high as 1530m/s
# (see http://staff.washington.edu/aganse/public.projects/
#            clustering/clustering.html),
# but that would only happen with h>225000m in this equation. Why ?
# The maximal speed we specify is really related to the max speed
# of surface pertubation
#
# v_max = 100 #For use in domain_ext.c
# sound_speed = 500

################################################################################
# Ranges specific to the shallow water wave equation
# These control maximal and minimal values of quantities
################################################################################

# Water depth below which it is considered to be 0 in the model
minimum_allowed_height = 1.0e-05

# Water depth below which it is *stored* as 0
minimum_storable_height = 1.0e-03

# FIXME (Ole): Redefine this parameter to control maximal speeds in general
# and associate it with protect_against_isolated_degenerate_timesteps = True
maximum_allowed_speed = 0.0 # Maximal particle speed of water
#maximum_allowed_speed = 1.0 # Maximal particle speed of water
                             # Too large (100) creates 'flopping' water
                             # Too small (0) creates 'creep'

maximum_froude_number = 100.0 # To be used in limiters.

################################################################################
# Performance parameters used to invoke various optimisations
################################################################################

use_psyco = False      # Use psyco optimisations

optimise_dry_cells = True # Exclude dry and still cells from flux computation
optimised_gradient_limiter = True # Use hardwired gradient limiter

points_file_block_line_size = 1e6 # Number of lines read in from a points file
                                  # when blocking

################################################################################
# NetCDF-specific type constants.  Used when defining NetCDF file variables.
################################################################################

netcdf_char = 'c'
netcdf_byte = 'b'
netcdf_int = 'i'
netcdf_float = 'd'
netcdf_float64 = 'd'
netcdf_float32 = 'f'

################################################################################
# Dynamically-defined constants.
################################################################################

# Determine if we can read/write large NetCDF files
netcdf_mode_w = 'w'
netcdf_mode_a = 'a'
netcdf_mode_r = 'r'


indent = '    '

# Code to set the write mode depending on
# whether Scientific.IO supports large NetCDF files
s = """
import os, tempfile
from anuga.file.netcdf import NetCDFFile

filename = tempfile.mktemp('.nc')

fid = NetCDFFile(filename, 'wl')
fid.close()
os.remove(filename)
"""

"""
# Need to run in a separate process due an
# error with older versions of Scientific.IO
if sys.platform == 'win32':
    null = 'NUL'
else:
    null = '/dev/null'
cmd = 'python -c "%s" 2> %s' % (s, null)
err = os.system(cmd)

if err != 0:
    # The Python script s failed e.g. with a segfault
    # which means that large file support is
    # definitely not supported
    pass
else:
    # Try the import within this process
    try:
        exec(s)
    except IOError:
        # NetCDFFile does not segfault but it does not
        # support large file support
        pass
    else:
        # Set the default mode to large file support
        #netcdf_mode_w = 'w4' # Future use of HDF5
        netcdf_mode_w = 'wl' # Large NetCDF (new default 30/6/2009)
        #netcdf_mode_w = 'w'   # Old style NetCDF used by OSG viewer

"""
