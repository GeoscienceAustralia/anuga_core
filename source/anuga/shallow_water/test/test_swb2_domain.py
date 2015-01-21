"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os

import anuga

from anuga import Reflective_boundary
from anuga import rectangular_cross_domain

from anuga import Domain

#from anuga_tsunami import Domain

import numpy as num
import warnings
import time



class Test_swb2_domain(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['runup_sinusoid_v2.sww']:
            try:
                os.remove(file)
            except:
                pass 


    def test_runup_sinusoid(self):
        """ Run a version of the validation test runup_sinusoid
        to ensure limiting solution has small velocity
        """

        points, vertices, boundary = anuga.rectangular_cross(20,20, len1=1., len2=1.)


        domain=Domain(points,vertices,boundary)    # Create Domain
        domain.set_flow_algorithm('DE0')
        
        domain.set_name('runup_sinusoid_v2')                         # Output to file runup.sww
        domain.set_datadir('.')                          # Use current folder
        domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
        #domain.set_store_vertices_uniquely(True)
        #------------------
        # Define topography
        #------------------
        scale_me=1.0

        def topography(x,y):
            return (-x/2.0 +0.05*num.sin((x+y)*50.0))*scale_me

        def stagefun(x,y):
            stge=-0.2*scale_me #+0.01*(x>0.9)
            return stge

        domain.set_quantity('elevation',topography)     # Use function for elevation
        domain.get_quantity('elevation').smooth_vertex_values()
        domain.set_quantity('friction',0.03)            # Constant friction


        domain.set_quantity('stage', stagefun)             # Constant negative initial stage
        domain.get_quantity('stage').smooth_vertex_values()


        #--------------------------
        # Setup boundary conditions
        #--------------------------
        Br=anuga.Reflective_boundary(domain)                 # Solid reflective wall
        Bd=anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])   # Constant boundary values -- not used in this example

        #----------------------------------------------
        # Associate boundary tags with boundary objects
        #----------------------------------------------
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})

        #------------------------------
        #Evolve the system through time
        #------------------------------

        for t in domain.evolve(yieldstep=7.0,finaltime=7.0):
            #print domain.timestepping_statistics()
            xx = domain.quantities['xmomentum'].centroid_values
            yy = domain.quantities['ymomentum'].centroid_values
            dd = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values
            #dd_raw=1.0*dd
            dd = (dd)*(dd>1.0e-03)+1.0e-03
            vv = ( (xx/dd)**2 + (yy/dd)**2)**0.5
            vv = vv*(dd>1.0e-03)
            #print 'Peak velocity is: ', vv.max(), vv.argmax()
            #print 'Volume is', sum(dd_raw*domain.areas)


        #print vv.max()

        assert num.all(vv<1.0e-02)


            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb2_domain, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
