"""
Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga
import numpy
from math import sin, pi, exp, sqrt
from anuga import Domain

#---------
#Setup computational domain
#---------
domain = anuga.rectangular_cross_domain(15,15, len1=1., len2=1.)


#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)

domain.set_name('dimensional_lid_driven')   # Output to file runup.sww
domain.set_datadir('.')         # Use current folder
#domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
domain.set_minimum_allowed_height(0.01)


#------------------
# Define topography
#------------------
domain.set_quantity('elevation',0.0)      # Use function for elevation
domain.set_quantity('friction',0.0)       # Constant friction
domain.set_quantity('stage',1.)          # Constant negative initial stage


#--------------------------
# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
Bd=anuga.Dirichlet_boundary([1., 1., 0.])       # Constant boundary values -- not used in this example

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Br, 'top': Bd, 'bottom':Br})


#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
parameter_file=open('parameters.tex', 'w')
parameter_file.write('\\begin{verbatim}\n')
from pprint import pprint
pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
parameter_file.write('\\end{verbatim}\n')
parameter_file.close()

#------------------------------
#Evolve the system through time
#------------------------------
for t in domain.evolve(yieldstep=10.,finaltime=1000.0):
    print domain.timestepping_statistics()
print 'Finished'
