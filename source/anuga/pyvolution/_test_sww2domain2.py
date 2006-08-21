"""Simple example of shallow water wave equation using Pyvolution

Water driven by linear slope and Dirichlet boundary

"""

######################
# Module imports 
#
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height, Time_boundary, Transmissive_boundary
from Numeric import array

#Create basic mesh
points, vertices, boundary = rectangular(10, 10)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.visualise = False
domain.store = True
domain.filename = 'bedslope'
domain.default_order=2
domain.quantities_to_be_stored=['stage']

#######################
#Bed-slope and friction
domain.set_quantity('elevation', lambda x,y: -x/3)
domain.set_quantity('friction', 0.1)

######################
# Boundary conditions
from math import sin, pi
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
Bd = Dirichlet_boundary([0.2,0.,0.])
Bw = Time_boundary(domain=domain,
                   f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

######################
#Initial condition
h = 0.05
elevation = domain.quantities['elevation'].vertex_values
domain.set_quantity('stage', elevation + h)
#elevation = domain.get_quantity('elevation',location='unique vertices')
#domain.set_quantity('stage', elevation + h,location='unique vertices')

domain.check_integrity()
######################
#Evolution
for t in domain.evolve(yieldstep = 1, finaltime = 2.0):
    pass
    #domain.write_time()


##NOW TEST IT!!!

from data_manager import sww2domain
from Numeric import allclose

filename = domain.datadir+'\\'+domain.filename+'.sww'

try:
    domain2 = sww2domain(filename,verbose=False)
    assert True == False
except:
    filler = 0
    domain2 = sww2domain(filename,fail_if_NaN=False,NaN_filler = filler,verbose=False)

bits = ['xllcorner','yllcorner','vertex_coordinates','time','starttime']

for quantity in ['elevation']+domain.quantities_to_be_stored:
    bits.append('get_quantity("%s")'%quantity)

for bit in bits:
#    print 'testing that domain.'+bit+' has been restored'
    assert allclose(eval('domain.'+bit),eval('domain2.'+bit))

#print max(max(domain2.get_quantity('xmomentum')))
#print min(min(domain2.get_quantity('xmomentum')))
#print max(max(domain2.get_quantity('ymomentum')))
#print min(min(domain2.get_quantity('ymomentum')))

assert max(max(domain2.get_quantity('xmomentum')))==filler
assert min(min(domain2.get_quantity('xmomentum')))==filler
assert max(max(domain2.get_quantity('ymomentum')))==filler
assert min(min(domain2.get_quantity('ymomentum')))==filler

#print 'passed'

#cleanup
#import os
#os.remove(domain.datadir+'/'+domain.filename+'.sww')