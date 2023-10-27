import anuga
import os
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

import numpy as num
import warnings
import time
import math

from anuga.shallow_water.sw_domain_cuda import nvtxRangePush, nvtxRangePop


nx = 10
ny = 10

def create_domain(name='domain'):

    domain = anuga.rectangular_cross_domain(nx, ny, len1=1., len2=1.)

    domain.set_flow_algorithm('DE0')
    domain.set_low_froude(0)

    domain.set_name(name)  
    domain.set_datadir('.')

    #------------------
    # Define topography
    #------------------
    scale_me=1.0

    #def topography(x,y):
    #    return (-x/2.0 +0.05*num.sin((x+y)*50.0))*scale_me

    def topography(x,y):
        return 0.0

    #def stagefun(x,y):
    #    stage=-0.2*scale_me #+0.01*(x>0.9)
    #    return stage

    def stagefun(x,y):
        stage=1.0-0.5*x
        return stage

    domain.set_quantity('elevation',topography)     # Use function for elevation
    domain.set_quantity('friction',0.03)            # Constant friction
    domain.set_quantity('stage', stagefun)          # Constant negative initial stage

    #--------------------------
    # Setup boundary conditions
    #--------------------------
    Br=anuga.Reflective_boundary(domain)                 # Solid reflective wall
    Bd=anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])   # Constant boundary values -- not used in this example

    #----------------------------------------------
    # Associate boundary tags with boundary objects
    #----------------------------------------------
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})
    #print(domain.__dict__)
    #print(dir(domain))
    return domain


print('')
print(70*'=')
print('Test Runup for update_conserved_quantities')
print(70*'=')
import time
start = time.time()

nvtxRangePush('create domain1')
domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(1)
nvtxRangePop()

nvtxRangePush('create domain1')
domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(1) # will change to 2 once burn in
nvtxRangePop()

end_time =time.time()
print('total creation time', (end_time - start))

#------------------------------
#Evolve the system through time
#------------------------------
yieldstep = 0.0002
finaltime = 0.0002
nvtxRangePush('evolve domain1')
print('Evolve domain1')
print('domain1 number of triangles ',domain1.number_of_elements)
for t in domain1.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain1.print_timestepping_statistics()
nvtxRangePop()

nvtxRangePush('update conserved quantities : domain1')
domain1.update_conserved_quantities()
nvtxRangePop()

nvtxRangePush('evolve domain2')
print('Evolve domain2')
print('domain2 number of triangles ',domain2.number_of_elements)
for t in domain2.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain2.print_timestepping_statistics()
nvtxRangePop()


nvtxRangePush('initialise gpu_interface : domain2')
domain2.set_multiprocessor_mode(4)
nvtxRangePop()

# import pdb; pdb.set_trace()
nvtxRangePush('update conserved quantities kernal : domain2')

nvtxRangePop()

