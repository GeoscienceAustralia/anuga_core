#!/usr/bin/env python

import unittest


from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
from anuga.structures.boyd_pipe_operator import boyd_pipe_function
from anuga.structures.boyd_pipe_operator import circular_critical_depth

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
import numpy
import math

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

        culvert_blockage = 0.0
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0

        # v and d re-baselined when the critical depth became the exact
        # circular-section solution (Q unchanged). The original culvertw
        # values (v=0.78, d=0.66) imply Froude 0.34 at a depth labelled
        # critical, i.e. they encoded the Q/sqrt(g)*D**2.5 transcription.
        expected_Q = 0.50
        expected_v = 1.64
        expected_d = 0.38


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)

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

        culvert_blockage = 0.0
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0

        # v and d re-baselined when the critical depth became the exact
        # circular-section solution (Q unchanged). The original culvertw
        # values (v=2.13, d=0.96) imply Froude 0.68 at a depth labelled
        # critical, i.e. they encoded the Q/sqrt(g)*D**2.5 transcription.
        expected_Q = 2.08
        expected_v = 2.62
        expected_d = 0.79


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew2')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)

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

        culvert_blockage=0.0
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0

        # v and d re-baselined when the critical depth became the exact
        # circular-section solution (Q unchanged). The original culvertw
        # values (v=4.94, d=1.20) had the pipe forced full; the exact
        # critical depth at this discharge is 0.965*D.
        expected_Q = 5.59
        expected_v = 5.00
        expected_d = 1.16


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew3')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)


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

        culvert_blockage = 0.0
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0

        # v and d re-baselined when the critical depth became the exact
        # circular-section solution (Q unchanged). The original culvertw
        # values (v=2.13, d=0.96) imply Froude 0.68 at a depth labelled
        # critical, i.e. they encoded the Q/sqrt(g)*D**2.5 transcription.
        expected_Q = 2.08
        expected_v = 2.62
        expected_d = 0.79


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew4')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)


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

        culvert_blockage = 0.0
        #culvert_barrels = 1.0

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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew5')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)


        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v

    def test_boyd_non_skew6(self):
        """test_boyd_non_skew

        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson,
        calculation code by MJ Boyd
        This tests the blockage code
        """

        stage_0 = 15.0 #change
        stage_1 = 14.0 #change
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2

        culvert_blockage = 0.50
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0


        expected_Q = 1.75
        expected_v = 3.11
        expected_d = 0.85


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew6')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)


        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2) #depth at outlet used to calc v

    def test_boyd_non_skew7(self):
        """test_boyd_non_skew

        This tests the Boyd routine with data obtained from culvertw application 1.1 by IceMindserer  BD Parkinson,
        calculation code by MJ Boyd
        This tests the blockage code
        """

        stage_0 = 15.0 #change
        stage_1 = 14.0 #change
        elevation_0 = 11.0
        elevation_1 = 10.0

        domain_length = 200.0
        domain_width = 200.0

        culvert_length = 20.0
        culvert_width = 1.2

        culvert_blockage = 1.0
        #culvert_barrels = 1.0

        culvert_losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
        culvert_mannings = 0.013

        culvert_apron = 0.0
        enquiry_gap = 5.0


        expected_Q = 0.0
        expected_v = 0.0
        expected_d = 0.0


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
                                    blockage=culvert_blockage,
                                    end_points=[ep0, ep1],
                                    #barrels=culvert_barrels,
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
            print('test_boyd_non_skew7')
            print('Q: ', Q, 'expected_Q: ', expected_Q)
            print('v: ', v, 'expected_v: ', expected_v)
            print('d: ', d, 'expected_d: ', expected_d)


        assert numpy.allclose(Q, expected_Q, rtol=1.0e-2, atol=1.0e-5) #inflow
        assert numpy.allclose(v, expected_v, rtol=1.0e-2, atol=1.0e-5) #outflow velocity
        assert numpy.allclose(d, expected_d, rtol=1.0e-2, atol=1.0e-5) #depth at outlet used to calc v


class Test_boyd_pipe_function_physics(unittest.TestCase):
    """Physics checks on boyd_pipe_function that do not depend on reference
    tables: exact critical depth, dimensional (Froude) similarity, barrel
    additivity and critical-flow consistency."""

    g = 9.8

    @staticmethod
    def _segment(y, D):
        theta = 2.0*math.acos(1.0 - 2.0*y/D)
        area = D*D/8.0*(theta - math.sin(theta))
        top_width = D*math.sin(theta/2.0)
        return area, top_width

    def test_circular_critical_depth_exact(self):
        """Q**2/g == A**3/T at the returned depth, and close to Straub."""
        import anuga
        g = anuga.g
        for D in [0.3, 0.6, 0.9, 1.5, 2.1]:
            for q_star in [0.05, 0.2, 0.5, 0.8]:
                Q = q_star*math.sqrt(g)*D**2.5
                y = circular_critical_depth(Q, D, g)
                self.assertTrue(0.0 < y < D)
                A, T = self._segment(y, D)
                self.assertAlmostEqual(A**3/T/(Q*Q/g), 1.0, places=8)
                # Straub (SI): dc = 1.01/D**0.264 * (Q/sqrt(g))**0.506,
                # stated valid for 0.02 < dc/D < 0.85
                y_straub = 1.01/D**0.264*(Q/math.sqrt(g))**0.506
                if 0.02 < y_straub/D < 0.85:
                    self.assertLess(abs(y/y_straub - 1.0), 0.05)
        self.assertEqual(circular_critical_depth(0.0, 1.0, g), 0.0)

    def _cases(self):
        # (driving_energy, delta_total_energy, outlet_enquiry_depth) for a
        # D = 0.9 m pipe; hits inlet control part-full, outlet control
        # part-full and outlet control with a submerged outlet.
        return [(0.6, 1.0, 0.0),
                (1.2, 0.3, 0.3),
                (1.5, 0.2, 1.2)]

    def test_froude_similarity(self):
        """Scaling all lengths by s (manning by s**(1/6)) must scale
        Q by s**2.5, velocity by s**0.5 and depth by s. Every term in the
        Boyd method is dimensionally consistent, so this holds exactly; a
        non-dimensionless critical-depth group breaks it."""
        D, L, n, losses = 0.9, 20.0, 0.013, 1.5
        for E, dE, tw in self._cases():
            Q0, v0, d0, a0, case0 = boyd_pipe_function(
                0.0, D, 0.0, 1.0, L, E, dE, tw, losses, n)
            for s in [0.5, 2.0]:
                Q1, v1, d1, a1, case1 = boyd_pipe_function(
                    0.0, s*D, 0.0, 1.0, s*L, s*E, s*dE, s*tw, losses,
                    n*s**(1.0/6.0))
                self.assertEqual(case1, case0)
                # barrel_velocity carries an absolute regularisation
                # (velocity_protection), so allow 1e-3 rather than exact.
                self.assertAlmostEqual(Q1/(Q0*s**2.5), 1.0, delta=1.0e-3)
                self.assertAlmostEqual(v1/(v0*s**0.5), 1.0, delta=1.0e-3)
                self.assertAlmostEqual(d1/(d0*s), 1.0, delta=1.0e-3)

    def test_barrels_are_additive(self):
        """N identical barrels carry exactly N times the flow of one, at the
        same velocity and depth."""
        D, L, n, losses = 0.9, 20.0, 0.013, 1.5
        for E, dE, tw in self._cases():
            Q1, v1, d1, a1, _ = boyd_pipe_function(
                0.0, D, 0.0, 1.0, L, E, dE, tw, losses, n)
            Q3, v3, d3, a3, _ = boyd_pipe_function(
                0.0, D, 0.0, 3.0, L, E, dE, tw, losses, n)
            self.assertAlmostEqual(Q3/(3.0*Q1), 1.0, delta=1.0e-3)
            self.assertAlmostEqual(v3/v1, 1.0, delta=1.0e-3)
            self.assertAlmostEqual(d3/d1, 1.0, delta=1.0e-3)

    def test_inlet_control_outlet_is_critical(self):
        """Inlet control with a free outlet: the reported depth and discharge
        are a critical-flow pair (Froude number 1)."""
        import anuga
        D, L, n, losses = 0.9, 20.0, 0.013, 1.5
        Q, v, d, area, case = boyd_pipe_function(
            0.0, D, 0.0, 1.0, L, 0.6, 1.0, 0.0, losses, n)
        self.assertIn('INLET CTRL', case)
        A, T = self._segment(d, D)
        self.assertAlmostEqual(A, area, places=9)
        self.assertAlmostEqual(Q*Q*T/(anuga.g*A**3), 1.0, places=6)


# =========================================================================
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_boyd_pipe_operator)
    runner = unittest.TextTestRunner()
    runner.run(suite)
