"""
Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga
import numpy
from math import sin, pi, exp
from anuga import Domain
from anuga import myid, finalize, distribute


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    #---------
    #Setup computational domain
    #---------
    points, vertices, boundary = anuga.rectangular_cross(100,3, len1=1.0,len2=0.03)
    domain=Domain(points,vertices,boundary)    # Create Domain
    domain.set_name('runup')                   # Output to file runup.sww
    domain.set_datadir('.')                    # Use current folder
    domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
    domain.set_flow_algorithm(alg)
        
    #------------------
    # Define topography
    #------------------
    def topography(x,y):
            return -x/2 #Linear bed slope
    
    def stagefun(x,y):
        return -0.45 #Stage
    
    domain.set_quantity('elevation',topography)     # Use function for elevation
    domain.get_quantity('elevation').smooth_vertex_values() # Steve's fix -- without this, substantial artificial velcities are generated everywhere in the domain. With this fix, there are artificial velocities near the coast, but not elsewhere.
    domain.set_quantity('friction',0.0)             # Constant friction
    domain.set_quantity('stage', stagefun)          # Constant negative initial stage
else:
    domain = None
    
#--------------------------
# create Parallel Domain    
#--------------------------
domain = distribute(domain)

# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
Bt=anuga.Transmissive_boundary(domain)          # Continue all values of boundary -- not used in this example
Bd=anuga.Dirichlet_boundary([-0.1, 0., 0.])     # Constant boundary values

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

#------------------------------
#Evolve the system through time
#------------------------------
#xwrite=open("xvel.out","wb")
#ywrite=open("yvel.out","wb")
## Set print options to be compatible with file writing via the 'print' statement 
#numpy.set_printoptions(threshold=numpy.nan, linewidth=numpy.nan)


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

for t in domain.evolve(yieldstep=0.2,finaltime=30.0):
    if myid == 0 and verbose: print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)

finalize()

