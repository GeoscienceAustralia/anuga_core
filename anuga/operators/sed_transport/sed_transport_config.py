"""Module where shared sediment transport parameters and default values are set
"""

import os
import sys
from numpy import sqrt

################################################################################
# Values set by user
################################################################################

momentum_sink_set = set()
# distances = []

# default values of operator function arguments
erosion = True
deposition = True
vegfile = None
sed_dispersion = False
turbulence = False
momentum_sinks = False
verbose = False

max_rate_e = 0.001               # maximum rate of sediment entrainment, m/s
max_rate_d = 0.001               # maximum rate of deposition, m/s

min_depth = 0.1                # minimum depth for calculating transport
max_conc = 0.5                  # maximum allowed conc before failure
default_boundary_flux = 0		# Default boundary flux, in Ch

graindiameter = 0.0001			# grain diameter of moving sediment, meters
porosity = 0.3					# porosity of the surface material

frict = 0.05					# friction factor
criticalshear_star = 0.15		# dimensionless critical shear stress
erosioncoeff_star = 2			# coefficient of dimensionless erosion equation
rousecoeff = 1					# coefficient, sed C profile in water
beta = 0.1						# effective speed of sediment in the flow

c1 = 18							# constants from settling velocity equation
c2 = 0.4

Cb = 0.001						# bed drag coefficient
Cd = 1.2                        # drag coefficient of cylinders

# operator_in_use = 0             # flag to know if an operator in sed_transport_operator has already been called
# operators_added = 0             # counter for number of operators added that might need to be removed

################################################################################
# Numerical constants
################################################################################

rho_w = 1000						# density of water, kg / m^3
rho_s = 2650						# density of sediment, kg / m^3
g = 9.81							# gravitational constant, m / s^2
mu = 1.e-6							# kinematic viscosity of water, m^2/s




################################################################################
# Calculated values
################################################################################

s = rho_s/rho_w - 1
R = graindiameter * sqrt(s * g * graindiameter) / mu
rhoo = (rho_w * porosity + rho_s * (1-porosity))

# dimensionless -> SI erosion coefficient
ero_coeff = erosioncoeff_star * (graindiameter * sqrt(((rho_s - rho_w) / rho_w) * g * graindiameter))

# coefficient for dimensionless shear stress
ss_coeff = rho_w / ((rho_s - rho_w) * g * graindiameter)

# settling velocity, m/s
settlingvelocity = (s * g * graindiameter**2) / ((c1 * mu) + (0.75 * c2 * (s * g * graindiameter**3)**0.5))

# coefficient for momentum sink/source from concentration gradients
cg_coeff = - (rho_s - rho_w) * g / 2