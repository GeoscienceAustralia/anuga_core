"""Module where global model parameters are set for anuga_1d
"""

epsilon = 1.0e-12
h0 = 1.0e-12

default_boundary_tag = 'exterior'


time_format = '%d/%m/%y %H:%M:%S'

min_timestep = 1.0e-6 #Should be computed based on geometry
max_timestep = 1.0e3
max_smallsteps = 50  # Max number of degenerate steps allowed b4 trying first order
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


v_max = 100 #For use in domain_ext.c
sound_speed = 500


max_smallsteps = 50  #Max number of degenerate steps allowed b4 trying first order

manning = 0.0  #Manning's friction coefficient
g = 9.8       #Gravity
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
#There are separate betas for the w-limiter and the h-limiter
#
#
#
#
#Good values are:
#beta_w = 0.9
#beta_h = 0.2



beta_w = 1.5
beta_h = 0.2
timestepping_method = 'euler'

CFL = 1.0  #FIXME (ole): Is this in use yet??
           #(Steve) yes, change domain.CFL to
           #make changes


pmesh_filename = '.\\pmesh'


import os, sys

if sys.platform == 'win32':
    #default_datadir = 'C:\grohm_output'
    default_datadir = '.'
else:
    #default_datadir = os.path.expanduser('~'+os.sep+'grohm_output')
    default_datadir = '.'


use_extensions = True    #Try to use C-extensions
#use_extensions = False   #Do not use C-extensions

use_psyco = True  #Use psyco optimisations
#use_psyco = False  #Do not use psyco optimisations


optimised_gradient_limiter = True #Use hardwired gradient limiter

#Specific to shallow water W.E.
h_min = minimum_allowed_height = 1.0e-6 #1.0e-6 #Water depth below which it is considered to be 0
maximum_allowed_speed = 1000.0
