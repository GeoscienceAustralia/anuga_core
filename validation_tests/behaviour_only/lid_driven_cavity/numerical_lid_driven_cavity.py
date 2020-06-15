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
from anuga import myid, finalize, distribute

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    #---------
    #Setup computational domain
    #---------
    domain = anuga.rectangular_cross_domain(15,15, len1=1., len2=1.)
    domain.set_flow_algorithm(alg)
 
    domain.set_name('dimensional_lid_driven')   # Output to file runup.sww
    domain.set_datadir('.')         # Use current folder
   
    
    #------------------
    # Define topography
    #------------------
    domain.set_quantity('elevation',0.0)      # Use function for elevation
    domain.set_quantity('friction',0.0)       # Constant friction
    domain.set_quantity('stage',1.)          # Constant negative initial stage
else:
    domain = None
    
#--------------------------
# Create Parallel Domain
#--------------------------
domain = distribute(domain)

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
if myid == 0:
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()

#------------------------------
#Evolve the system through time
#------------------------------
for t in domain.evolve(yieldstep=10.,finaltime=200.0):
    if myid == 0 and verbose: print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)

finalize()
