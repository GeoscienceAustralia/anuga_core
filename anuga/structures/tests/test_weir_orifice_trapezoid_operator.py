#!/usr/bin/env python


import unittest


from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator
from anuga.structures.weir_orifice_trapezoid_operator import weir_orifice_trapezoid_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
import numpy
import inspect




verbose = False


class Test_weir_orifice_trapezoid_operator(unittest.TestCase):
    """
	Test the weir_orifice_trapezoid_operator, in particular the discharge_routine!
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
    
    #def test_weir_orifice_non_skew(self):

        #stage_0 = 11.0
        #stage_1 = 10.0
        #elevation_0 = 10.0
        #elevation_1 = 10.0

        #domain_length = 200.0
        #domain_width = 200.0

        #culvert_length = 10.0
        #culvert_width = 10.0
        #culvert_height = 3.0
        #culvert_slope=0.1
        #culvert_z1 = 2
        #culvert_z2 = 2
        #culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        #culvert_mannings = 0.015
        
        #culvert_apron = 0.0
        #enquiry_gap = 5.0

        
        #expected_Q = 0.52 
        #expected_v = 1.04
        #expected_d = 0.50
                                    

        #domain = self._create_domain(d_length=domain_length,
                                     #d_width=domain_width,
                                     #dx = 5.0,
                                     #dy = 5.0,
                                     #elevation_0 = elevation_0,
                                     #elevation_1 = elevation_1,
                                     #stage_0 = stage_0,
                                     #stage_1 = stage_1)
 

        ##print 'Defining Structures'
        
        #ep0 = numpy.array([domain_length/2-culvert_length/2, 100.0])
        #ep1 = numpy.array([domain_length/2+culvert_length/2, 100.0])
        
        
        #culvert = Weir_orifice_trapezoid_operator(domain,
                                    #losses=culvert_losses,
                                    #width=culvert_width,
                                    #end_points=[ep0, ep1],
                                    #height=culvert_height,
                                    #culvert_slope=culvert_slope,
                                    #z1=culvert_z1,
                                    #z2=culvert_z2,
                                    #apron=culvert_apron,
                                    #enquiry_gap=enquiry_gap,
                                    #use_momentum_jet=False,
                                    #use_velocity_head=False,
                                    #manning=culvert_mannings,
                                    #logging=False,
                                    #label='10mwidex3mdeep 1V:2H trapezoid',
                                    #verbose=False)

        ##culvert.determine_inflow_outflow()
        
        #( Q, v, d ) = culvert.discharge_routine()
        
        #if verbose:
            #print 50*'='
            #print 'UNITTEST ',inspect.stack()[0][3]
            #print 'Q v d: ', Q, v, d, ' expected_Q v d: ', expected_Q, expected_v, expected_d
        

        #assert numpy.allclose(Q, expected_Q, rtol=1.0e-1) #inflow
        #assert numpy.allclose(v, expected_v, rtol=1.0e-1) #outflow velocity
        #assert numpy.allclose(d, expected_d, rtol=1.0e-1) #depth at outlet used to calc v 
        
      

    def test_weir_orifice_1(self):
        """test_weir_orifice 1
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs    
        """


        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=0.57
        outlet_depth=0.45
        inlet_velocity=1.93 # v jut upstream of culvert
        outlet_velocity=2.03 # v just downstream of culvert
        
        culvert_length=10.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_slope=0.1
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 0.52 
        v_expected = 1.04
        d_expected = 0.05
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3]
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)
        




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %('   ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 
        
    def test_weir_orifice_2(self):
        """test_weir_orifice 2
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs  
        """


        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=1.99
        outlet_depth=1.9
        inlet_velocity=3.59
        outlet_velocity=3.82
        
        culvert_length=10.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_slope=0.01
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 0.09
        v_expected = 0.20
        d_expected = 0.04
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3] 
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %(   'ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 
        
    def test_weir_orifice_3(self):
        """test_weir_orifice 3
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs 
        """


        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=6.15
        outlet_depth=2.02
        inlet_velocity=1.82
        outlet_velocity=8.82
        
        culvert_length=10.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 6.44
        v_expected = 2.33
        d_expected = 0.24
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3]
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %('   ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v  

    def test_weir_orifice_4(self):
        """test_weir_orifice 4
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs  
        """


        g=9.81
        culvert_slope=0.1  # Downward

        inlet_depth=12.0
        outlet_depth=10.0
        inlet_velocity=2.0
        outlet_velocity=2.5
        
        culvert_length=100.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 383.75
        v_expected = 5.09
        d_expected = 3.00
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3]
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %('   ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 
        

    def test_weir_orifice_5(self):
        """test_weir_orifice 5
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs  
        """


        g=9.81
        culvert_slope=1  # Downward

        inlet_depth=12.0
        outlet_depth=10.0
        inlet_velocity=2.0
        outlet_velocity=2.5
        
        culvert_length=100.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 487.48
        v_expected = 6.13
        d_expected = 3.00
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3]
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %('   ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 

    def test_weir_orifice_6(self):
        """test_weir_orifice 6
        
        This tests the test_weir_orifice routine with data obtained from a spreadsheet to do the weir orfice calcs 
        """


        g=9.81
        culvert_slope=10  # Downward

        inlet_depth=12.0
        outlet_depth=10.0
        inlet_velocity=2.0
        outlet_velocity=2.5
        
        culvert_length=100.0
        culvert_width=10.0
        culvert_height=3.0
        culvert_z1=2
        culvert_z2=2
        
        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g 
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out
        

        Q_expected = 1049.86
        v_expected = 8.75
        d_expected = 3.00
        
        if verbose:
            print 50*'='
            print 'UNITTEST ',inspect.stack()[0][3]
            print 'width ',culvert_width
            print 'depth ',culvert_height
            print 'flow_width ',culvert_width
            print 'length ' ,culvert_length
            print 'driving_energy ',inlet_specific_energy
            print 'delta_total_energy ',delta_total_energy
            print 'outlet_enquiry_depth ',outlet_depth
            print 'sum_loss ',sum_loss
            print 'manning ',manning
        
        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(culvert_width, 
                                                    culvert_height, 
                                                    culvert_width,
                                                    culvert_slope,
                                                    culvert_z1,
                                                    culvert_z2, 
                                                    culvert_length, 
                                                    inlet_specific_energy, 
                                                    delta_total_energy, 
                                                    outlet_depth, 
                                                    sum_loss,
                                                    manning)




        if verbose:
            print ('%s %.2f'%('SPEC_E = ',inlet_specific_energy))
            print ('%s %.2f'%('Delta E = ',delta_total_energy))
            print ('%s %.2f,%.2f,%.2f' %('   ANUGAcalcsTEST Q-v-d ',Q,v,d))
            print ('%s %.2f,%.2f,%.2f' %('spreadsheet calculation ', Q_expected, v_expected, d_expected))
            
        assert numpy.allclose(Q, Q_expected, rtol=1.0e-1) #inflow
        assert numpy.allclose(v, v_expected, rtol=1.0e-1) #outflow velocity
        assert numpy.allclose(d, d_expected, rtol=1.0e-1) #depth at outlet used to calc v 
# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_weir_orifice_trapezoid_operator, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
