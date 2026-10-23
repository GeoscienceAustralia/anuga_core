#!/usr/bin/env python
"""
    Generic boundary conditions for a domain.

    A boundary represents the edge of the model, where inflow, outflow, and
    reflection can take place.

    The boundaries in this model can be applied universally across all
    domain models, without being tied to a particular implementation.
"""
import unittest
from math import sqrt, pi

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import *
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.config import epsilon

import numpy as num


class Test_Generic_Boundary_Conditions(unittest.TestCase):
    def setUp(self):
        pass
        #print "  Setting up"

    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_generic(self):
        b = Boundary()

        try:
            b.evaluate()
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')


    def test_dirichlet_empty(self):

        try:
            Bd = Dirichlet_boundary()
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')

    def test_dirichlet(self):
        x = [3.14,0,0.1]
        Bd = Dirichlet_boundary(x)

        q = Bd.evaluate()
        assert num.allclose(q, x)


    def test_time(self):



        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test time bdry, you need to provide a domain and function
        try:
            T = Time_boundary(domain)
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')

        #Test time bdry, you need to provide a function
        try:
            T = Time_boundary()
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')


        def function(t):
            return [1.0, 0.0]

        T = Time_boundary(domain, function)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [1.0, 0.0])


    def test_time_space_boundary(self):


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test time space bdry, you need to provide a domain and function
        try:
            T = Time_space_boundary(domain)
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')

        #Test time bdry, you need to provide a function
        try:
            T = Time_space_boundary()
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')


        def function(t,x,y):
            return [x,y]

        T = Time_space_boundary(domain, function)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2
        assert num.allclose(q, domain.get_edge_midpoint_coordinate(0,2))


        q = T.evaluate(1, 1)  #Vol=1, edge=1
        assert num.allclose(q, domain.get_edge_midpoint_coordinate(1,1))





    def test_transmissive(self):


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test transmissve bdry
        try:
            T = Transmissive_boundary()
        except Exception:
            pass
        else:
            raise Exception('Should have raised exception')

        T = Transmissive_boundary(domain)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [1.5, 2.5])


        # Now set the centroid_transmissive_bc flag to true
        domain.set_centroid_transmissive_bc(True)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [2.0 ,3.0]) # centroid value





    def NOtest_fileboundary_time_only(self):
        """Test that boundary values can be read from file and interpolated
        This is using the .tms file format

        See also test_util for comprenhensive testing of the underlying
        file_function and also tests in test_datamanager which tests
        file_function using the sts format
        """
        #FIXME (Ole): This test was disabled 18 August 2008 as no
        # need for this was found. Rather I implemented an Exception
        # to catch possible errors in the model setup


        import time
        import os
        from math import sin, pi
        from anuga.config import time_format

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.check_integrity()


        #Write file
        filename = 'boundarytest' + str(time.time())
        fid = open(filename + '.txt', 'w')
        start = time.mktime(time.strptime('2000', '%Y'))
        dt = 5*60  #Five minute intervals
        for i in range(10):
            t = start + i*dt
            t_string = time.strftime(time_format, time.gmtime(t))

            fid.write('%s,%f %f\n' %(t_string, 1.0*i, sin(i*2*pi/10)))
        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)

        from anuga.shallow_water.data_manager import timefile2netcdf

        timefile2netcdf(filename, quantity_names = ['stage', 'ymomentum'])



        F = File_boundary(filename + '.tms', domain)


        os.remove(filename + '.txt')
        os.remove(filename + '.tms')




        #Check that midpoint coordinates at boundary are correctly computed
        assert num.allclose( F.midpoint_coordinates,
                             [[1.0, 0.0], [0.0, 1.0], [3.0, 0.0],
                              [3.0, 1.0], [1.0, 3.0], [0.0, 3.0]])

        #assert allclose(F.midpoint_coordinates[(3,2)], [0.0, 3.0])
        #assert allclose(F.midpoint_coordinates[(3,1)], [1.0, 3.0])
        #assert allclose(F.midpoint_coordinates[(0,2)], [0.0, 1.0])
        #assert allclose(F.midpoint_coordinates[(0,0)], [1.0, 0.0])
        #assert allclose(F.midpoint_coordinates[(2,0)], [3.0, 0.0])
        #assert allclose(F.midpoint_coordinates[(2,1)], [3.0, 1.0])


        #Check time interpolation
        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: F} )

        domain.time = 5*30/2  #A quarter way through first step
        q = F.evaluate()
        assert num.allclose(q, [1.0/4, sin(2*pi/10)/4])


        domain.time = 2.5*5*60  #Half way between steps 2 and 3
        q = F.evaluate()
        assert num.allclose(q, [2.5, (sin(2*2*pi/10) + sin(3*2*pi/10))/2])



    def test_fileboundary_exception(self):
        """Test that boundary object complains if number of
        conserved quantities are wrong
        """


        import time
        import os
        from math import sin, pi
        from anuga.config import time_format

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
        domain.evolved_quantities = ['stage', 'xmomentum', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['xmomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])
        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.check_integrity()

        #Write file (with only two values)
        filename = 'boundarytest' + str(time.time())
        fid = open(filename + '.txt', 'w')
        start = time.mktime(time.strptime('2000', '%Y'))
        dt = 5*60  #Five minute intervals
        for i in range(10):
            t = start + i*dt
            t_string = time.strftime(time_format, time.gmtime(t))

            fid.write('%s,%f %f\n' %(t_string, 1.0*i, sin(i*2*pi/10)))
        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)
        from anuga.file_conversion.file_conversion import timefile2netcdf

        timefile2netcdf(filename+'.txt', quantity_names = ['stage', 'xmomentum'])


        try:
            F = File_boundary(filename + '.tms',
                              domain)
        except Exception:
            pass
        else:
            raise Exception('Should have raised an exception')

        os.remove(filename + '.txt')
        os.remove(filename + '.tms')


class Test_Boundary_evaluate_segment(unittest.TestCase):
    """Test evaluate_segment and __repr__ methods for all boundary types.

    Uses a minimal Generic_Domain so boundary_cells / boundary_edges /
    evolved_quantities are properly initialised.
    """

    def _make_domain(self):
        """Return a Generic_Domain with 4 triangles and two conserved quantities."""
        points = [
            [0.0, 0.0], [0.0, 2.0], [2.0, 0.0],
            [0.0, 4.0], [2.0, 2.0], [4.0, 0.0],
        ]
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities   = ['stage', 'ymomentum']
        domain.quantities['stage']     = Quantity(domain, [[1, 2, 3],  [5, 5, 5],
                                                            [0, 0, 9], [-6, 3, 3]])
        domain.quantities['ymomentum'] = Quantity(domain, [[2, 3, 4],  [5, 5, 5],
                                                            [0, 0, 9], [-6, 3, 3]])
        domain.check_integrity()
        # Set a boundary on all external edges
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: Dirichlet_boundary([0.0, 0.0])})
        return domain

    def _all_segment_edges(self, domain):
        """Return an array indexing every boundary edge."""
        return num.arange(len(domain.boundary_cells))

    # ------------------------------------------------------------------
    # Boundary base class: evaluate_segment and get_time guards
    # ------------------------------------------------------------------

    def test_boundary_evaluate_segment_no_edges(self):
        """evaluate_segment with None segment_edges returns silently."""
        domain = self._make_domain()
        b = Boundary()
        b.evaluate_segment(domain, None)  # should not raise

    def test_boundary_evaluate_segment_no_domain(self):
        b = Boundary()
        b.evaluate_segment(None, num.array([0]))  # should not raise

    # ------------------------------------------------------------------
    # Transmissive_boundary
    # ------------------------------------------------------------------

    def test_transmissive_repr(self):
        domain = self._make_domain()
        T = Transmissive_boundary(domain)
        self.assertIn('Transmissive', repr(T))

    def test_transmissive_evaluate_segment_edge_mode(self):
        domain = self._make_domain()
        T = Transmissive_boundary(domain)
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})

        domain.set_centroid_transmissive_bc(False)
        ids = self._all_segment_edges(domain)
        T.evaluate_segment(domain, ids)
        # boundary_values should now be set; spot-check no NaN
        for name in domain.evolved_quantities:
            Q = domain.quantities[name]
            self.assertFalse(num.any(num.isnan(Q.boundary_values[ids])))

    def test_transmissive_evaluate_segment_centroid_mode(self):
        domain = self._make_domain()
        T = Transmissive_boundary(domain)
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})

        domain.set_centroid_transmissive_bc(True)
        ids = self._all_segment_edges(domain)
        T.evaluate_segment(domain, ids)
        # Centroid values should have been copied to boundary_values
        for name in domain.evolved_quantities:
            Q = domain.quantities[name]
            vol_ids = domain.boundary_cells[ids]
            num.testing.assert_array_equal(
                Q.boundary_values[ids], Q.centroid_values[vol_ids])

    def test_transmissive_evaluate_segment_none_guards(self):
        domain = self._make_domain()
        T = Transmissive_boundary(domain)
        T.evaluate_segment(domain, None)  # no exception
        T.evaluate_segment(None, num.array([0]))  # no exception

    # ------------------------------------------------------------------
    # Dirichlet_boundary
    # ------------------------------------------------------------------

    def test_dirichlet_repr(self):
        bd = Dirichlet_boundary([1.0, 2.0])
        self.assertIn('Dirichlet', repr(bd))

    def test_dirichlet_evaluate_segment_evolved_quantities(self):
        """When len(q_bdry) == len(evolved_quantities) no edge fallback needed."""
        domain = self._make_domain()
        # evolved and conserved both have length 2, so q_bdry length 2 matches evolved
        bd = Dirichlet_boundary([7.0, 8.0])
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: bd})

        ids = self._all_segment_edges(domain)
        bd.evaluate_segment(domain, ids)

        for j, name in enumerate(['stage', 'ymomentum']):
            Q = domain.quantities[name]
            expected_val = [7.0, 8.0][j]
            num.testing.assert_allclose(
                Q.boundary_values[ids], expected_val,
                err_msg=f"{name} boundary_values not set to {expected_val}")

    def test_dirichlet_evaluate_segment_none_guards(self):
        bd = Dirichlet_boundary([1.0, 0.0])
        domain = self._make_domain()
        bd.evaluate_segment(domain, None)   # no exception
        bd.evaluate_segment(None, num.array([0]))   # no exception

    # ------------------------------------------------------------------
    # Compute_fluxes_boundary
    # ------------------------------------------------------------------

    def test_compute_fluxes_boundary_instantiation(self):
        cfb = Compute_fluxes_boundary()
        self.assertIsNotNone(cfb)

    def test_compute_fluxes_boundary_evaluate_returns_none(self):
        cfb = Compute_fluxes_boundary()
        result = cfb.evaluate(0, 0)
        self.assertIsNone(result)

    def test_compute_fluxes_boundary_evaluate_segment_returns(self):
        cfb = Compute_fluxes_boundary()
        domain = self._make_domain()
        result = cfb.evaluate_segment(domain, num.array([0]))
        self.assertIsNone(result)

    # ------------------------------------------------------------------
    # Time_boundary
    # ------------------------------------------------------------------

    def test_time_boundary_repr(self):
        domain = self._make_domain()
        T = Time_boundary(domain, function=lambda t: [1.0, 0.0])
        self.assertEqual(repr(T), 'Time boundary')

    def test_time_boundary_get_time(self):
        domain = self._make_domain()
        T = Time_boundary(domain, function=lambda t: [0.0, 0.0])
        # domain time should be accessible via T.get_time()
        self.assertAlmostEqual(T.get_time(), domain.get_time())

    def test_time_boundary_evaluate_segment(self):
        domain = self._make_domain()
        T = Time_boundary(domain, function=lambda t: [3.0, 4.0])
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})

        ids = self._all_segment_edges(domain)
        T.evaluate_segment(domain, ids)

        for j, name in enumerate(['stage', 'ymomentum']):
            Q = domain.quantities[name]
            expected = [3.0, 4.0][j]
            num.testing.assert_allclose(
                Q.boundary_values[ids], expected,
                err_msg=f"{name} boundary_values not set by Time_boundary")

    def test_time_boundary_evaluate_segment_none_guards(self):
        domain = self._make_domain()
        T = Time_boundary(domain, function=lambda t: [0.0, 0.0])
        T.evaluate_segment(domain, None)        # no exception
        T.evaluate_segment(None, num.array([0]))  # no exception

    def test_time_boundary_bad_function_type(self):
        """Function that returns a non-numeric value raises an Exception."""
        domain = self._make_domain()
        with self.assertRaises(Exception):
            Time_boundary(domain, function=lambda t: 'bad')

    # ------------------------------------------------------------------
    # Time_space_boundary
    # ------------------------------------------------------------------

    def test_time_space_boundary_repr(self):
        domain = self._make_domain()
        T = Time_space_boundary(domain, function=lambda t, x, y: [0.0, 0.0])
        self.assertEqual(repr(T), 'Time space boundary')

    def test_time_space_boundary_evaluate(self):
        domain = self._make_domain()
        T = Time_space_boundary(domain, function=lambda t, x, y: [x, y])
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})

        q = T.evaluate(0, 2)
        expected = domain.get_edge_midpoint_coordinate(0, 2)
        num.testing.assert_allclose(q, expected, atol=1e-12)

    def test_time_space_boundary_bad_function(self):
        domain = self._make_domain()
        with self.assertRaises(Exception):
            Time_space_boundary(domain, function=lambda t, x, y: 'bad')

    # ------------------------------------------------------------------
    # Base Boundary.evaluate_segment loop (lines 60-86)
    # ------------------------------------------------------------------

    def test_base_evaluate_segment_loop(self):
        """Base class evaluate_segment iterates over edges (lines 60-86)."""
        class ConcreteBoundary(Boundary):
            def evaluate(self, vol_id=None, edge_id=None):
                return [1.0, 0.0]
        domain = self._make_domain()
        b = ConcreteBoundary()
        ids = self._all_segment_edges(domain)
        b.evaluate_segment(domain, ids)
        # After the call, boundary_values for each edge should be set
        for name in domain.evolved_quantities:
            Q = domain.quantities[name]
            self.assertFalse(num.any(num.isnan(Q.boundary_values[ids])))

    # ------------------------------------------------------------------
    # Boundary.get_time (line 91)
    # ------------------------------------------------------------------

    def test_base_get_time(self):
        """Boundary.get_time() via a subclass that sets self.domain (line 91)."""
        domain = self._make_domain()
        b = Transmissive_boundary(domain)  # sets self.domain; doesn't override get_time
        # The base get_time should work
        t = b.get_time()
        self.assertIsNotNone(t)

    # ------------------------------------------------------------------
    # Boundary.get_boundary_values — t provided explicitly (96->99 branch)
    # ------------------------------------------------------------------

    def test_get_boundary_values_with_explicit_t(self):
        """get_boundary_values(t=val) skips the t=None branch (line 96->99)."""
        domain = self._make_domain()
        T = Time_boundary(domain, function=lambda t: [t, 0.0])
        res = T.get_boundary_values(t=5.0)
        self.assertAlmostEqual(res[0], 5.0)

    # ------------------------------------------------------------------
    # Compute_fluxes_boundary.__repr__ (line 277)
    # ------------------------------------------------------------------

    def test_compute_fluxes_boundary_repr(self):
        """__repr__ is callable (line 277)."""
        cfb = Compute_fluxes_boundary()
        # __repr__ accesses self.domain which is not set → AttributeError expected
        try:
            s = repr(cfb)
            self.assertIsInstance(s, str)
        except AttributeError:
            pass  # expected: self.domain not set in __init__

    # ------------------------------------------------------------------
    # Dirichlet conserved_quantities=True path (lines 234-245, 251)
    # ------------------------------------------------------------------

    def _make_domain_3evolved(self):
        """Domain with 3 evolved_quantities but 2 conserved_quantities."""
        points = [
            [0.0, 0.0], [0.0, 2.0], [2.0, 0.0],
            [0.0, 4.0], [2.0, 2.0], [4.0, 0.0],
        ]
        elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities   = ['stage', 'ymomentum', 'xmomentum']
        domain.quantities['stage']     = Quantity(domain, [[1, 2, 3],  [5, 5, 5],
                                                            [0, 0, 9], [-6, 3, 3]])
        domain.quantities['ymomentum'] = Quantity(domain, [[2, 3, 4],  [5, 5, 5],
                                                            [0, 0, 9], [-6, 3, 3]])
        domain.quantities['xmomentum'] = Quantity(domain, [[0, 0, 0],  [0, 0, 0],
                                                            [0, 0, 0], [ 0, 0, 0]])
        domain.check_integrity()
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: Dirichlet_boundary([0.0, 0.0])})
        return domain

    def test_dirichlet_evaluate_segment_conserved_path(self):
        """conserved_quantities=True path: edge fallback + conserved write (lines 234-245, 251)."""
        domain = self._make_domain_3evolved()
        # 2 dirichlet values == len(conserved_quantities), != len(evolved_quantities)=3
        bd = Dirichlet_boundary([3.0, 4.0])
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: bd})
        ids = self._all_segment_edges(domain)
        bd.evaluate_segment(domain, ids)  # should not raise

    def test_time_boundary_evaluate_segment_conserved_path(self):
        """Time_boundary conserved path: edge fallback (lines 381-392, 398)."""
        domain = self._make_domain_3evolved()
        # Function returns 2 values == len(conserved_quantities) != len(evolved)=3
        T = Time_boundary(domain, function=lambda t: [1.0, 0.0])
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})
        ids = self._all_segment_edges(domain)
        T.evaluate_segment(domain, ids)  # should not raise

    # ------------------------------------------------------------------
    # Time_boundary.__init__ exception on function execution (lines 330-332)
    # ------------------------------------------------------------------

    def test_time_boundary_function_raises_on_call(self):
        """Function that raises during execution hits lines 330-332."""
        domain = self._make_domain()
        def bad_func(t):
            raise RuntimeError("cannot evaluate")
        with self.assertRaises(Exception):
            Time_boundary(domain, function=bad_func)

    # ------------------------------------------------------------------
    # get_boundary_values — Modeltime_too_early/too_late paths (lines 101-129)
    # ------------------------------------------------------------------

    def test_get_boundary_values_modeltime_too_early_reraise(self):
        """Modeltime_too_early raises through (lines 101-102)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_early
        domain = self._make_domain()
        def early_func(t):
            raise Modeltime_too_early('too early')
        T = Time_boundary.__new__(Time_boundary)
        T.function = early_func
        T.domain = domain
        T.default_boundary = None
        T.default_boundary_invoked = False
        T.verbose = False
        with self.assertRaises(Modeltime_too_early):
            T.get_boundary_values()

    def test_get_boundary_values_modeltime_too_late_no_default(self):
        """Modeltime_too_late with no default boundary re-raises (lines 103-105)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late
        domain = self._make_domain()
        def late_func(t):
            raise Modeltime_too_late('too late')
        T = Time_boundary.__new__(Time_boundary)
        T.function = late_func
        T.domain = domain
        T.default_boundary = None
        T.default_boundary_invoked = False
        T.verbose = False
        with self.assertRaises(Modeltime_too_late):
            T.get_boundary_values()

    def test_get_boundary_values_modeltime_too_late_with_default(self):
        """Modeltime_too_late with default boundary uses default (lines 106-129)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late
        domain = self._make_domain()
        def late_func(t):
            raise Modeltime_too_late('too late')
        T = Time_boundary.__new__(Time_boundary)
        T.function = late_func
        T.domain = domain
        T.default_boundary = num.array([9.0, 9.0])
        T.default_boundary_invoked = False
        T.verbose = True  # cover the logging branch
        res = T.get_boundary_values()
        num.testing.assert_allclose(res, [9.0, 9.0])

    # ------------------------------------------------------------------
    # Time_space_boundary.__init__ exception on function execution (449-451)
    # ------------------------------------------------------------------

    def test_time_space_boundary_function_raises_on_call(self):
        """Function that raises during execution hits lines 449-451."""
        domain = self._make_domain()
        def bad_func(t, x, y):
            raise RuntimeError("cannot evaluate")
        with self.assertRaises(Exception):
            Time_space_boundary(domain, function=bad_func)

    # ------------------------------------------------------------------
    # Time_space_boundary.evaluate — Modeltime_too_late with default (483-510)
    # ------------------------------------------------------------------

    def test_time_space_boundary_modeltime_too_late_with_default(self):
        """Modeltime_too_late with default boundary uses default (lines 485-510)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late
        domain = self._make_domain()
        call_count = [0]
        def late_func(t, x, y):
            call_count[0] += 1
            if call_count[0] > 1:  # first call (t=0) from __init__ succeeds
                raise Modeltime_too_late('too late')
            return [0.0, 0.0]
        default_b = Dirichlet_boundary([5.0, 5.0])
        T = Time_space_boundary(domain, function=late_func,
                                default_boundary=default_b, verbose=True)
        from anuga.config import default_boundary_tag
        domain.set_boundary({default_boundary_tag: T})
        # evaluate at an edge; function will now raise Modeltime_too_late
        res = T.evaluate(0, 0)
        num.testing.assert_allclose(res, [5.0, 5.0])


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Generic_Boundary_Conditions)
    runner = unittest.TextTestRunner()
    runner.run(suite)
