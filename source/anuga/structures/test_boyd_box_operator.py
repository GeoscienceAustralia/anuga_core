#!/usr/bin/env python


import unittest


from anuga.structures.boyd_box_operator import Boyd_box_operator
from anuga.structures.boyd_box_operator import boyd_box_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
import numpy

verbose = False


class Test_boyd_box_operator(unittest.TestCase):
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

    def test_boyd_non_skew(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 6.23
        expected_v = 2.55
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
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x3.6RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
        

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v 
        
        
        
        
        
    
    def test_boyd_skew(self):
        """test_boyd_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0
        
        expected_Q = 6.23
        expected_v = 2.55
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
        
        a = domain_length/2 - culvert_length/2
        b = domain_length/2 + culvert_length/2
        
        el0 = numpy.array([[a, 100.0 - culvert_width/2], [a, 100.0 + culvert_width/2]])
        el1 = numpy.array([[b, 100.0 - culvert_width/2], [b, 100.0 + culvert_width/2]])
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    exchange_lines=[el0, el1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    label='3.6x3.6RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_skew'
            print 'Q: ', Q, 'expected_Q: ', expected_Q

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v         


    def test_boyd_non_skew_enquiry_points(self):
        """test_boyd_skew_enquiry_points
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        stage_0 = 11.0
        stage_1 = 10.0
        elevation_0 = 10.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 3.66
        culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0
        
        expected_Q = 6.23
        expected_v = 2.55
        expected_d = 0.66
        

        # Probably no need to change below here
        
        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1)


        #print 'Defining Structures'
        
        a = domain_length/2 - culvert_length/2
        b = domain_length/2 + culvert_length/2
        
        el0 = numpy.array([[a, 100.0 - culvert_width/2], [a, 100.0 + culvert_width/2]])
        el1 = numpy.array([[b, 100.0 - culvert_width/2], [b, 100.0 + culvert_width/2]])
        
        enquiry_points = (numpy.array([85, 100]), numpy.array([115, 100]))
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    exchange_lines=[el0, el1],
                                    enquiry_points=enquiry_points,
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=culvert_mannings,
                                    label='3.6x3.6RCBC',
                                    verbose=False)


        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_non_skew_enquiry_points'
            print 'Q: ', Q, 'expected_Q: ', expected_Q
            print 'v: ', v, 'expected_v: ', expected_v
            print 'd: ', d, 'expected_d: ', expected_d

        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v         

    def test_boyd_1(self):
        """test_boyd_1
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd     
        """


        g=9.81


        inlet_depth=0.150
        outlet_depth=0.15
        inlet_velocity=1.00
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.6
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        culvert_slope=10.0  # % Downward
        z_in = 10.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 


        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning


        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)


        if verbose:
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST01 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 0.5526, 1.146, 0.1339))
            
        assert numpy.allclose(Q, 0.5526, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, 1.146, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, 0.1339, rtol=1.0e-1) #depth at outlet used to calc v 
   
   
   
   
    def test_boyd_1_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81


        inlet_depth=0.150
        outlet_depth=0.150
        inlet_velocity=1.00
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.6
        culvert_height=1.20
        
        culvert_type='box'

        #sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        culvert_slope=10.0  # % Downward
        z_in = 10.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        
        expected_Q = 0.55
        expected_v = 1.15
        expected_d = 0.13


        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',expected_Q,expected_v, expected_d
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, expected_Q, rtol=5.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=5.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=5.0e-2) #depth at outlet used to calc v 
        
        
        
               
    def test_boyd_2(self):
        """test_boyd_2
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd     
        """


        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=0.500
        outlet_depth=0.700
        inlet_velocity=1.0
        outlet_velocity=0.50
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
        Q_expected = 2.50
        v_expected = 1.897
        d_expected = 0.367
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)



        if verbose:
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST02 Q-v-d',Q,v,d))
            print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 2.508, 1.897, 0.367))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v  


    def test_boyd_2_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=0.500
        outlet_depth=0.700
        inlet_velocity=1.0
        outlet_velocity=0.50
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0

        Q_expected = 2.50
        v_expected = 1.897
        d_expected = 0.367


        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_2'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        
        



    def test_boyd_3(self):
        """test_boyd_3
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd     
        """


        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.800
        outlet_depth=0.80
        inlet_velocity=1.0
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 13.554
        v_expected = 3.329
        d_expected = 1.131
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.2f'%('Delta E = ',delta_total_energy))
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST03 Q-v-d',Q,v,d))
            print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 13.554, 3.329, 1.131))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-2) #depth at outlet used to calc v 


    def test_boyd_3_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.800
        outlet_depth=0.80
        inlet_velocity=1.0
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 13.554
        v_expected = 3.329
        d_expected = 1.131
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        
        




#NOTE FROM HERE DOWN THE UNITS TEST HAVE NOT BEEN AMENDED TO ALLOW VELOCITY COMPONENT TO BE USED. ONLY ABOVE 3 TESTS WORK. PM WILL FIX THE ONES BELOW WHEN THE ABOVE 3 ARE WORKING
    def test_boyd_4(self):
        """test_boyd_4
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.00
        outlet_depth=0.8
        inlet_velocity=1.0
        outlet_velocity=0.5 
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
       
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
        
        
        Q_expected = 6.609
        v_expected = 2.621
        d_expected = 0.70
        
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning

        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)


        if verbose:
            print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.2f'%('Delta E = ',delta_total_energy))
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST04 Q-v-d',Q,v,d))
            print ('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', 6.609, 2.621, 0.70))
            
        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 


    def test_boyd_4_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.00
        outlet_depth=0.8
        inlet_velocity=1.0
        outlet_velocity=0.5 
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
       
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
        
        
        Q_expected = 6.609
        v_expected = 2.621
        d_expected = 0.70
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 

    def test_boyd_5(self):
        """test_boyd_5
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 1.0
        outlet_depth=2.5
        outlet_velocity=0.5
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
    
        Q_expected = 2.961
        v_expected = 0.685
        d_expected = 1.20
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning

        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)


 
        if verbose:
            print ('%s,%.3f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.3f'%('Delta E = ',delta_total_energy))
             
            print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST05 Q-v-d',Q,v,d))
            print ('%s,%.3f,%.3f,%.3f' %('Spreadsheet_Boydcalcs',2.961, 0.685, 1.20))
            
        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        
        
        
        
    def test_boyd_5_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

 
        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 1.0
        outlet_depth=2.5
        outlet_velocity=0.5
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
    
        Q_expected = 2.961
        v_expected = 0.685
        d_expected = 1.20
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 


    def test_boyd_6(self):
        """test_boyd_6
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd     
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.80
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        Q_expected = 15.537
        v_expected = 3.597
        d_expected = 1.20

        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:   
            print ('%s,%.3f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.3f'%('Delta E = ',delta_total_energy))
             
            print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST06 Q-v-d',Q,v,d))
            print ('%s,%.3f,%.3f,%.3f' %('Spreadsheet_Boydcalcs',15.537, 3.597, 1.20))

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        
        
    def test_boyd_6_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.80
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        Q_expected = 15.537
        v_expected = 3.597
        d_expected = 1.20
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 

    def test_boyd_7(self):
        """test_boyd_7
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81


        inlet_depth=0.150
        outlet_depth=0.15
        inlet_velocity=1.00
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.6
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        culvert_slope=1  # % Downward
        z_in = 10.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 


        Q_expected = 0.5526
        v_expected = 1.146
        d_expected = 0.1339
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)
        
#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)
        if verbose:
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST01 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_7_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81


        inlet_depth=0.150
        outlet_depth=0.15
        inlet_velocity=1.00
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.6
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        culvert_slope=1  # % Downward
        z_in = 10.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 


        Q_expected = 0.5526
        v_expected = 1.146
        d_expected = 0.1339
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 


        
    def test_boyd_8(self):
        """test_boyd_8
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=0.500
        outlet_depth=0.700
        inlet_velocity=1.50
        outlet_velocity=0.50
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
        
        Q_expected = 0.224
        v_expected = 0.152
        d_expected = 0.409

        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)

#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)
        if verbose:
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST02 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_8_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=0.500
        outlet_depth=0.700
        inlet_velocity=1.50
        outlet_velocity=0.50
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        
        
        Q_expected = 0.224
        v_expected = 0.152
        d_expected = 0.409
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 


    def test_boyd_9(self):
        """test_boyd_9
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.800
        outlet_depth=0.80
        inlet_velocity=1.0
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 13.554
        v_expected = 3.329
        d_expected = 1.131
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)


#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)

        if verbose:
            print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.2f'%('Delta E = ',delta_total_energy))
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST03 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_9_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.800
        outlet_depth=0.80
        inlet_velocity=1.0
        outlet_velocity=0.5
        
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 13.554
        v_expected = 3.329
        d_expected = 1.131
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 

    def test_boyd_10(self):
        """test_boyd_10
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.00
        outlet_depth=0.8
        inlet_velocity=1.0
        outlet_velocity=0.5 
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
       
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out



        Q_expected = 5.164
        v_expected = 2.047
        d_expected = 0.70
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)

#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)        

        if verbose:
            print ('%s,%.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.2f'%('Delta E = ',delta_total_energy))
            print ('%s,%.2f,%.2f,%.2f' %('ANUGAcalcsTEST04 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_10_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.00
        outlet_depth=0.8
        inlet_velocity=1.0
        outlet_velocity=0.5 
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
       
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out



        Q_expected = 5.164
        v_expected = 2.047
        d_expected = 0.70
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 

    def test_boyd_11(self):
        """test_boyd_11
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.50
        inlet_velocity= 1.0
        outlet_depth=1.3
        outlet_velocity=0.5
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 8.808
        v_expected = 2.039
        d_expected = 1.20
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)

#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)  
        if verbose:      
            print ('%s,%.3f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.3f'%('Delta E = ',delta_total_energy))
        
            print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST05Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_11_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.50
        inlet_velocity= 1.0
        outlet_depth=1.3
        outlet_velocity=0.5
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 8.808
        v_expected = 2.039
        d_expected = 1.20
        
      
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 


    def test_boyd_12(self):
        """test_boyd_12
        
        This tests the Boyd routine with data obtained from ??? by Petar Milevski    
        """
        # FIXME(Ole): This test fails (20 Feb 2009)

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.8
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 13.546
        v_expected = 3.136
        d_expected = 1.20
        
        if verbose:
            print 50*'='
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
            
        Q, v, d, flow_area, case= boyd_box_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)

#         Q, v, d = boyd_generalised_culvert_model(inlet_depth,
#                                                  outlet_depth,
#                                                  inlet_velocity,
#                                                  outlet_velocity,
#                                                  inlet_specific_energy, 
#                                                  delta_total_energy, 
#                                                  g,
#                                                  culvert_length,
#                                                  culvert_width,
#                                                  culvert_height,
#                                                  culvert_type,
#                                                  manning,
#                                                  sum_loss)        
        if verbose:
            print ('%s,%.3f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s,%.3f'%('Delta E = ',delta_total_energy))
        
            print ('%s,%.3f,%.3f,%.3f' %('ANUGAcalcsTEST06 Q-v-d',Q,v,d))
            print('%s,%.2f,%.2f,%.2f' %('Spreadsheet_Boydcalcs', Q_expected, v_expected, d_expected))

        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_boyd_12_operator(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.8
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        Q_expected = 13.546
        v_expected = 3.136
        d_expected = 1.20
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=False)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_1'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        

    def test_boyd_12_operator_invert(self):
        """test_boyd_non_skew
        
        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson, 
        calculation code by MJ Boyd 
        """

        #verbose = True
        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.50
        inlet_velocity= 4.0
        outlet_depth=0.8
        outlet_velocity=4.0
        culvert_length=10.0
        culvert_width=3.60
        culvert_height=1.20
        
        culvert_type='box'
        manning=0.013
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 10.0
        z_out = 10.0-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        invert0 = 20.0
        invert1 = 30.0


        #Q_expected = 13.546
        #v_expected = 3.136
        #d_expected = 1.20
        
        Q_expected = 0.0
        v_expected = 0.0
        d_expected = 0.0
        
        elevation_0 = z_in
        elevation_1 = z_out
        
        stage_0 = elevation_0 + inlet_depth
        stage_1 = elevation_1 + outlet_depth
 

        domain_length = 200.0
        domain_width = 200.0

        #culvert_length = 20.0
        #culvert_width = 3.66
        #culvert_height = 3.66
        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013
        
        culvert_apron = 0.0
        enquiry_gap = 5.0




        domain = self._create_domain(d_length=domain_length,
                                     d_width=domain_width,
                                     dx = 5.0,
                                     dy = 5.0,
                                     elevation_0 = elevation_0,
                                     elevation_1 = elevation_1,
                                     stage_0 = stage_0,
                                     stage_1 = stage_1,
                                     xvelocity_0 = inlet_velocity,
                                     xvelocity_1 = outlet_velocity)
 

        #print 'Defining Structures'
        
        ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        culvert = Boyd_box_operator(domain,
                                    losses=culvert_losses,
                                    width=culvert_width,
                                    end_points=[ep0, ep1],
                                    height=culvert_height,
                                    apron=culvert_apron,
                                    #invert_elevations=[elevation_0,elevation_1],
                                    invert_elevations=[invert0,invert1],
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=False,
                                    use_velocity_head=True,
                                    manning=culvert_mannings,
                                    logging=False,
                                    label='3.6x1.2RCBC',
                                    verbose=verbose)

        #culvert.determine_inflow_outflow()
        
        ( Q, v, d ) = culvert.discharge_routine()
        
        if verbose:
            print 'test_boyd_operator_12_invert'
            print 'expected ',Q_expected,v_expected, d_expected
            print 'calc ',Q,v,d
        

        assert numpy.allclose(Q, Q_expected, rtol=2.0e-2) #inflow
        assert numpy.allclose(v, v_expected, rtol=2.0e-2) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=2.0e-2) #depth at outlet used to calc v 
        
       
# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_boyd_box_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
