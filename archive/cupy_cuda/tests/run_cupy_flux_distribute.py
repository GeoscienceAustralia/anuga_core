import anuga
import os
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

import numpy as num
import warnings
import time
import math

from pprint import pprint

from anuga.shallow_water.sw_domain_cuda import nvtxRangePush, nvtxRangePop


nx = 50
ny = 50

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
print('Test Runup for extrapolate')
print(70*'=')

nvtxRangePush('create domain1')
domain1 = create_domain('domain_original')
domain1.set_multiprocessor_mode(1)

quantities1 = domain1.quantities
stage1 = quantities1["stage"]
xmom1 = quantities1["xmomentum"]
ymom1 = quantities1["ymomentum"]
nvtxRangePop()

nvtxRangePush('create domain2')
domain2 = create_domain('domain_cuda')
domain2.set_multiprocessor_mode(2)


quantities2 = domain2.quantities
stage2 = quantities2["stage"]
xmom2 = quantities2["xmomentum"]
ymom2 = quantities2["ymomentum"]
nvtxRangePop()

import time
start = time.time()

#------------------------------
#Evolve the system through time
#------------------------------
yieldstep = 0.002
finaltime = 0.02
nvtxRangePush('evolve domain1')
print('Evolve domain1')
print('domain1 number of triangles ',domain1.number_of_elements)
for t in domain1.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain1.print_timestepping_statistics()
nvtxRangePop()

end = time.time()
print('DOMAIN 1 time ' + str(end - start))



#-----------------------------------------
# Test the kernel version of compute fluxes
#----------------------------------------

start = time.time()


nvtxRangePush('evolve domain2')
print('Evolve domain2')
print('domain2 number of triangles ',domain2.number_of_elements)
for t in domain2.evolve(yieldstep=yieldstep,finaltime=finaltime):
    domain2.print_timestepping_statistics()
nvtxRangePop()

end = time.time()
print('DOMAIN 2 time ' + str(end - start))