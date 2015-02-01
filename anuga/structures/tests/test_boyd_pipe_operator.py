#!/usr/bin/env python


import unittest


from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
from anuga.structures.boyd_pipe_operator import boyd_pipe_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
import numpy

verbose = False
#diameter = width

class Test_boyd_pipe_operator(unittest.TestCase):
    """
	Test the boyd box operator, in particular the discharge_routine!
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    
    def _create_domain(self,d_length,
                            d_width,
                            dx,
                            dy,
                            elevation_0,
                            elevation_1,
                            stage_0,
                            stage_1,
                            xvelocity_0 = 0.0,
                            xvelocity_1 = 0.0,
                            yvelocity_0 = 0.0,
                            yvelocity_1 = 0.0):
        
        points, vertices, boundary = rectangular_cross(int(d_length/dx), int(d_width/dy),
                                                        len1=d_length, len2=d_width)
        domain = Domain(points, vertices, boundary)   
        domain.set_name('Test_Outlet_Inlet')                 # Output name
        domain.set_store()
        domain.set_default_order(2)
        domain.H0 = 0.01
        domain.tight_slope_limiters = 1

        #print 'Size', len(domain)

        #------------------------------------------------------------------------------
        # Setup initial conditions
        #------------------------------------------------------------------------------

        def elevation(x, y):
            """Set up a elevation
            """
            
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = elevation_0
            
            numpy.putmask(z, x > d_length/2, elevation_1)
    
            return z
            
        def stage(x,y):
            """Set up stage
            """
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = stage_0
            
            numpy.putmask(z, x > d_length/2, stage_1)

            return z
        
        def xmom(x,y):
            """Set up xmomentum
            """
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = xvelocity_0*(stage_0-elevation_0)
            
            numpy.putmask(z, x > d_length/2, xvelocity_1*(stage_1-elevation_1) )

            return z
        
        def ymom(x,y):
            """Set up ymomentum
            """
            z = numpy.zeros(x.shape,dtype='d')
            z[:] = yvelocity_0*(stage_0-elevation_0)
            
            numpy.putmask(z, x > d_length/2, yvelocity_1*(stage_1-elevation_1) )

            return z
            
        #print 'Setting Quantities....'
        domain.set_quantity('elevation', elevation)  # Use function for elevation
        domain.set_quantity('stage',  stage)   # Use function for elevation
        domain.set_quantity('xmomentum',  xmom) 
        domain.set_quantity('ymomentum',  ymom) 
        
        return domain

    def test_boyd_non_skew1(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.5
        stage_1 = 10.0
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2
        ##culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 0.50
        expected_v = 0.78
        expected_d = 0.66
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_pipe_operator(domain,
                                    losses=culvert_losses,
                                    diameter=culvert_width,
                                    end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='1.2pipe',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v 
        
    def test_boyd_non_skew2(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 12.2
        stage_1 = 10.0
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2
        ##culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 2.08
        expected_v = 2.13
        expected_d = 0.96
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_pipe_operator(domain,
                                    losses=culvert_losses,
                                    diameter=culvert_width,
                                    end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='1.2pipe',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew2'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v  
        
    def test_boyd_non_skew3(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 15.0
        stage_1 = 10.0
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2
        ##culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 5.59
        expected_v = 4.94
        expected_d = 1.20
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_pipe_operator(domain,
                                    losses=culvert_losses,
                                    diameter=culvert_width,
                                    end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='1.2pipe',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew3'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d

 
        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v 
        
    def test_boyd_non_skew4(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 12.2 #change
        stage_1 = 11.2 #change
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2
        ##culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 2.08
        expected_v = 2.13
        expected_d = 0.96
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_pipe_operator(domain,
                                    losses=culvert_losses,
                                    diameter=culvert_width,
                                    end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='1.2pipe',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew4'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d


        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v 
        
    def test_boyd_non_skew5(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 15.0 #change
        stage_1 = 14.0 #change
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2
        ##culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 3.70
        expected_v = 3.27
        expected_d = 1.20
        

        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_pipe_operator(domain,
                                    losses=culvert_losses,
                                    diameter=culvert_width,
                                    end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='1.2pipe',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew5'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d


        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v  
# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_boyd_pipe_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
