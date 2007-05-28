"""Module where global pyvolution model parameters are set
"""


#FIXME (Ole): Temporary access to global config file
from anuga_config import epsilon, default_boundary_tag



#FIXME (Ole): More of these may need to be moved to anuga_config.py
time_format = '%d/%m/%y %H:%M:%S'

min_timestep = 1.0e-6 #Should be computed based on geometry
max_timestep = 1.0e+3
#This is how:
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
sound_speed = 500


max_smallsteps = 50  #Max number of degenerate steps allowed b4 trying first order

manning = 0.03  #Manning's friction coefficient
#g = 9.80665       #Gravity
g = 9.8
#g(phi) = 9780313 * (1 + 0.0053024 sin(phi)**2 - 0.000 0059 sin(2*phi)**2) micro m/s**2, where phi is the latitude
#The 'official' average is 9.80665




eta_w = 3.0e-3 #Wind stress coefficient
rho_a = 1.2e-3 #Atmospheric density
rho_w = 1023   #Fluid density [kg/m^3] (rho_w = 1023 for salt water)


#Betas [0;1] control the allowed steepness of gradient for second order
#extrapolations. Values of 1 allow the steepes gradients while
#lower values are more conservative. Values of 0 correspond to
#1'st order extrapolations.
#
# Large values of beta_h may cause simulations to require more timesteps
# as surface will 'hug' closer to the bed.
# Small values of beta_h will make code faster, but one may experience
# artificial momenta caused by discontinuities in water depths in
# the presence of steep slopes. One example of this would be
# stationary water 'lapping' upwards to a higher point on the coast.
#
#
#
#There are separate betas for the w, uh, vh and h limiters
#
#Good values are:


# I think these are better SR but they conflict with the unit tests!
beta_w      = 1.0
beta_w_dry  = 0.2
beta_uh     = 1.0
beta_uh_dry = 0.2
beta_vh     = 1.0
beta_vh_dry = 0.2
beta_h      = 0.2

# beta_h can be safely put to zero esp if we are using limit2007 = 1. This will
# also speed things up.
beta_h = 0.0


# Alpha_balance controls how limiters are balanced between deep and shallow.
# A large value will favour the deep water limiters, allowing the a closer hug to the coastline.
# This will minimise 'creep' but at the same time cause smaller time steps
# Range:

alpha_balance = 2.0 

# Flag use of new limiters.
# limit2007 = 0 means use old limiters (e.g. for some tests)
# limit2007 = 1 means use new limiters that hug the bathymetry closer
limit2007 = 0



CFL = 1.0  #FIXME (ole): Is this in use yet??
           #(Steve) yes, change domain.CFL to
           #make changes


pmesh_filename = '.\\pmesh'
version_filename = 'stored_version_info.py'


import os, sys

if sys.platform == 'win32':
    default_datadir = '.'
else:
    default_datadir = '.'


use_extensions = True    #Try to use C-extensions
#use_extensions = False   #Do not use C-extensions

use_psyco = True  #Use psyco optimisations
#use_psyco = False  #Do not use psyco optimisations


optimised_gradient_limiter = True #Use hardwired gradient limiter

#Specific to shallow water W.E.
minimum_allowed_height = 1.0e-3 #Water depth below which it is considered to be 0 in the model

maximum_allowed_speed = 1.0 # Maximal particle speed of water
                            # Too large (100) creates 'flopping' water
                            # Too small (0) creates 'creep'


minimum_storable_height = 1.0e-5 #Water depth below which it is *stored* as 0

points_file_block_line_size = 500 # Number of lines read in from a points file
                                  # when blocking

umask = 002  # used to set file and directory permission created by anuga
