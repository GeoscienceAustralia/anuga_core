"""Module where global ANUGA model parameters and default values are set
"""

#--------------------
# Numerical constants
#--------------------
epsilon = 1.0e-12 # Smallest number - used for safe division
max_float = 1.0e36 # Largest number - used to initialise (max, min) ranges


#-------------------------------------------
# Standard filenames, directories and system 
# parameters used by ANUGA
#-------------------------------------------
pmesh_filename = '.\\pmesh'
version_filename = 'stored_version_info.py'
default_datadir = '.'
time_format = '%d/%m/%y %H:%M:%S'
umask = 002  # Controls file and directory permission created by anuga
default_boundary_tag = 'exterior' 

# Major revision number for use with create_distribution
# and update_anuga_user_guide
major_revision = '1.0beta'


#-------------------
# Physical constants
#-------------------
manning = 0.03  # Manning's friction coefficient
#g = 9.80665    # Gravity - FIXME reinstate this and fix unit tests.
g = 9.8
#g(phi) = 9780313 * (1 + 0.0053024 sin(phi)**2 - 0.000 0059 sin(2*phi)**2) micro m/s**2, where phi is the latitude
#The 'official' average is 9.80665

eta_w = 3.0e-3 #Wind stress coefficient
rho_a = 1.2e-3 #Atmospheric density
rho_w = 1023   #Fluid density [kg/m^3] (rho_w = 1023 for salt water)


#-----------------------------------------------------
# Limiters - used with linear reconstruction of vertex 
# values from centroid values
#-----------------------------------------------------

# Betas [0;1] control the allowed steepness of gradient for second order
# extrapolations. Values of 1 allow the steepes gradients while
# lower values are more conservative. Values of 0 correspond to
# 1'st order extrapolations.
#
# Large values of beta_h may cause simulations to require more timesteps
# as surface will 'hug' closer to the bed.
# Small values of beta_h will make code faster, but one may experience
# artificial momenta caused by discontinuities in water depths in
# the presence of steep slopes. One example of this would be
# stationary water 'lapping' upwards to a higher point on the coast.
#
# NOTE (Ole): I believe this was addressed with the introduction of
# tight_slope_limiters. I wish to retire the beta_? parameters.
# Can you please let me know if you disagree?

# There are separate betas for the w, uh, vh and h limiters
# I think these are better SR but they conflict with the unit tests!
beta_w      = 1.0
beta_w_dry  = 0.2
beta_uh     = 1.0
beta_uh_dry = 0.2
beta_vh     = 1.0
beta_vh_dry = 0.2
beta_h      = 0.2

# beta_h can be safely put to zero esp if we are using
# tight_slope_limiters = 1. This will
# also speed things up in general
beta_h = 0.0


# Alpha_balance controls how limiters are balanced between deep and shallow.
# A large value will favour the deep water limiters, allowing the a closer hug to the coastline.
# This will minimise 'creep' but at the same time cause smaller time steps
# Range:

alpha_balance = 2.0 

# Flag use of new limiters.
# tight_slope_limiters = 0 means use old limiters (e.g. for some tests)
# tight_slope_limiters = 1 means use new limiters that hug the bathymetry closer
tight_slope_limiters = 0



#-------------
# Timestepping
#-------------

CFL = 1.0  # CFL condition assigned to domain.CFL - controls timestep size
      
# Choose type of timestepping,
timestepping_method = 'euler' # 1st order euler
#timestepping_method = 'rk2'   # 2nd Order TVD scheme

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
max_smallsteps = 50  # Max number of degenerate steps allowed b4 trying first order

#Perhaps minimal timestep could be based on the geometry as follows:
#Define maximal possible speed in open water v_max, e.g. 500m/s (soundspeed?)
#Then work out minimal internal distance in mesh r_min and set
#min_timestep = r_min/v_max
#
#Max speeds are calculated in the flux function as
#
#lambda = v +/- sqrt(gh)
#
# so with 500 m/s, h ~ 500^2/g = 2500 m well out of the domain of the
# shallow water wave equation
#
#The actual soundspeed can be as high as 1530m/s
#(see http://staff.washington.edu/aganse/public.projects/clustering/clustering.html),
#but that would only happen with h>225000m in this equation. Why ?
#The maximal speed we specify is really related to the max speed
#of surface pertubation
#
#v_max = 100 #For use in domain_ext.c
#sound_speed = 500


#---------------------------------------------------
# Ranges specific to the shallow water wave equation
# These control maximal and minimal values of 
# quantities
#---------------------------------------------------

# Water depth below which it is considered to be 0 in the model
minimum_allowed_height = 1.0e-3 

# Water depth below which it is *stored* as 0
minimum_storable_height = 1.0e-5

# FIXME (Ole): Redefine this parameter to control maximal speeds in general
# and associate it with protect_against_isolated_degenerate_timesteps = True
maximum_allowed_speed = 0.0 # Maximal particle speed of water
#maximum_allowed_speed = 1.0 # Maximal particle speed of water
                            # Too large (100) creates 'flopping' water
                            # Too small (0) creates 'creep'
			    
maximum_froude_number = 100.0 # To be used in limiters.


#------------------------------------------------------------
# Performance parameters used to invoke various optimisations
#------------------------------------------------------------

use_extensions = True # Use C-extensions
use_psyco = True # Use psyco optimisations

optimise_dry_cells = True # Exclude dry and still cells from flux computation
optimised_gradient_limiter = True # Use hardwired gradient limiter

points_file_block_line_size = 500 # Number of lines read in from a points file
                                  # when blocking






    
