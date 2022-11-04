#!/usr/bin/env python


import unittest


from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator
from anuga.structures.weir_orifice_trapezoid_operator import weir_orifice_trapezoid_function

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
import numpy
import inspect




verbose =  False


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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        # from Petars spreadsheet
        #Q_expected = 8.791
        #v_expected = 1.116
        #d_expected = 0.692

        Q_expected = 4.95
        v_expected = 1.14
        d_expected = 0.40



        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_height ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_barrels ',culvert_barrels)
            print('culvert_z1 ',culvert_z1)
            print('culvert_z2 ',culvert_z2)
            print('culvert_flow_width ',culvert_width)
            print('culvert_length ' ,culvert_length)
            print('culvert_slope ',culvert_slope)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)
            print('outlet_depth ', outlet_depth)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_velocity ', outlet_velocity)







        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope      = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)



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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        # from petars spreadsheet
        #Q_expected = 8.245
        #v_expected = 0.231
        #d_expected = 2.410

        Q_expected = 0.82
        v_expected = 0.22
        d_expected = 0.35



        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_depth ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_length ' ,culvert_length)
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_depth ', outlet_depth)
            print('outlet_velocity ', outlet_velocity)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('flow_width ',culvert_width)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)

        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope     = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)




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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        # from petars spreadsheet
        #Q_expected = 113.323
        #v_expected = 2.361
        #d_expected = 3.0

        Q_expected = 67.64
        v_expected = 2.36
        d_expected = 2.04


        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_depth ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_length ' ,culvert_length)
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_depth ', outlet_depth)
            print('outlet_velocity ', outlet_velocity)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('flow_width ',culvert_width)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)

        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope     = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)




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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        # values from petars spreadsheet
        #Q_expected = 225.023
        #v_expected = 4.688
        #d_expected = 3.00

        Q_expected = 231.74
        v_expected = 4.83
        d_expected = 3.00


        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_depth ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_length ' ,culvert_length)
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_depth ', outlet_depth)
            print('outlet_velocity ', outlet_velocity)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('flow_width ',culvert_width)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)

        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope     = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)




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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out

        # values from petars spreadsheet
        #Q_expected = 271.275
        #v_expected = 5.652
        #d_expected = 3.00

        Q_expected = 279.38
        v_expected = 5.82
        d_expected = 3.00



        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_depth ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_length ' ,culvert_length)
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_depth ', outlet_depth)
            print('outlet_velocity ', outlet_velocity)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('flow_width ',culvert_width)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)

        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope     = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)




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
        culvert_blockage = 0.0
        culvert_barrels = 1.0

        culvert_type='trapezoid'
        manning=0.015
        sum_loss=1.5

        inlet_specific_energy=inlet_depth + 0.5*inlet_velocity**2/g
        z_in = 0.0
        z_out = z_in-culvert_length*culvert_slope/100
        E_in = z_in+inlet_depth + 0.5*inlet_velocity**2/g
        E_out = z_out+outlet_depth + 0.5*outlet_velocity**2/g
        delta_total_energy = E_in-E_out


        # Values from Petars spreadsheet
        #Q_expected = 420.160
        #v_expected = 11.470
        #d_expected = 3.00

        Q_expected = 419.95
        v_expected = 8.75
        d_expected = 3.00


        if verbose:
            print(50*'=')
            print('UNITTEST ',inspect.stack()[0][3])
            print('culvert_width ',culvert_width)
            print('culvert_depth ',culvert_height)
            print('culvert_blockage ',culvert_blockage)
            print('culvert_length ' ,culvert_length)
            print('inlet_depth ', inlet_depth)
            print('inlet_velocity ', inlet_velocity)
            print('outlet_depth ', outlet_depth)
            print('outlet_velocity ', outlet_velocity)
            print('sum_loss ',sum_loss)
            print('manning ',manning)
            print(' ')
            print('flow_width ',culvert_width)
            print('driving_energy ',inlet_specific_energy)
            print('delta_total_energy ',delta_total_energy)

        Q, v, d, flow_area, case= weir_orifice_trapezoid_function(
                                                    width      = culvert_width,
                                                    depth      = culvert_height,
                                                    blockage   = culvert_blockage,
                                                    barrels    = culvert_barrels,
                                                    z1         = culvert_z1,
                                                    z2         = culvert_z2,
                                                    flow_width = culvert_width,
                                                    length     = culvert_length,
                                                    #culvert_slope     = culvert_slope,
                                                    driving_energy     = inlet_specific_energy,
                                                    delta_total_energy = delta_total_energy,
                                                    outlet_enquiry_depth = outlet_depth,
                                                    sum_loss   = sum_loss,
                                                    manning    = manning)




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
    suite = unittest.makeSuite(Test_weir_orifice_trapezoid_operator, 'test_')
    runner = unittest.TextTestRunner()
    runner.run(suite)
