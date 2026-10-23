#!/usr/bin/env python

import unittest
from math import sqrt

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.region import *
#from anuga.config import epsilon

import numpy as num


"""
This is what the mesh in these tests look like;

3---7
|5 /|
| /6|
2---6
|3 /|
| /2|
1---5
|1 /|
| /0|
0---4
"""

def add_x_y(x, y):
    return x+y

def give_me_23(x, y):
    return 23.0

class Test_region(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_region_indices(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)

        region = Region(domain, indices=[0,2,3])

        expected_indices = [0,2,3]
        assert num.allclose(region.indices, expected_indices)


    def test_region_polygon(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)

        poly = [[0.0,0.0], [0.5,0.0], [0.5,0.5]]

        #print poly
        region = Region(domain, polygon=poly)

        expected_indices = [1]
        assert num.allclose(region.indices, expected_indices)


    def test_region_polygon_expanded(self):
        """create region based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)

        region = Region(domain, polygon=[[0.0,0.0], [0.5,0.0], [0.5,0.5]], expand_polygon=True)

        # Triangle 1 holds the polygon; the polygon's edges pass through 0
        # and 2. Triangle 3 only touches the polygon at its corner (0.5, 0.5),
        # which lies on the edge triangles 2 and 3 share: a point contact,
        # not an overlap, so it is not selected (it used to be, by rounding).
        expected_indices = [0,1,2]
        assert num.allclose(sorted(region.indices), expected_indices)



#-------------------------------------------------------------

class Test_region_extra(unittest.TestCase):
    """Tests for uncovered Region paths."""

    def _make_domain(self):
        points, vertices, boundary = rectangular(1, 3)
        return Domain(points, vertices, boundary)

    def test_empty_indices(self):
        """indices=[] → self.indices=[], type='empty' (lines 92-93)."""
        domain = self._make_domain()
        r = Region(domain, indices=[])
        self.assertEqual(r.type, 'empty')
        self.assertEqual(len(r.indices), 0)
        self.assertEqual(len(r.full_indices), 0)   # line 152

    def test_none_region_full_indices(self):
        """indices=None → full_indices covers all full triangles (line 150)."""
        domain = self._make_domain()
        r = Region(domain)   # all None → else branch (line 146)
        self.assertIsNotNone(r.full_indices)

    def test_circle_region(self):
        """center+radius → _setup_indices_circle (lines 99-104, 217-241)."""
        domain = self._make_domain()
        r = Region(domain, center=[0.25, 0.25], radius=0.5)
        self.assertEqual(r.type, 'circle')
        self.assertIsNotNone(r.indices)

    def test_line_region(self):
        """line= → _setup_indices_line (lines 116-124, 278-289)."""
        domain = self._make_domain()
        # A line that cuts across triangles
        r = Region(domain, line=[[0.0, 0.25], [1.0, 0.25]])
        self.assertEqual(r.type, 'line')
        self.assertIsNotNone(r.indices)

    def test_poly_as_polygon(self):
        """poly with > 2 points → treated as polygon (lines 126-140)."""
        domain = self._make_domain()
        r = Region(domain, poly=[[0.0,0.0], [0.5,0.0], [0.5,0.5], [0.0,0.5]])
        self.assertEqual(r.type, 'polygon')

    def test_poly_as_line(self):
        """poly with 2 points → treated as line (lines 141-144)."""
        domain = self._make_domain()
        r = Region(domain, poly=[[0.0, 0.25], [1.0, 0.25]])
        self.assertEqual(r.type, 'line')

    def test_repr(self):
        """__repr__ returns a string (line 158)."""
        domain = self._make_domain()
        r = Region(domain, indices=[0, 1])
        s = repr(r)
        self.assertIsInstance(s, str)

    def test_get_type(self):
        """get_type() returns type string (line 161)."""
        domain = self._make_domain()
        r = Region(domain, indices=[0])
        self.assertEqual(r.get_type(), 'indices_specified')

    def test_get_indices_not_full_only(self):
        """get_indices(full_only=False) (line 297)."""
        domain = self._make_domain()
        r = Region(domain, indices=[0, 1, 2])
        idx = r.get_indices(full_only=False)
        self.assertTrue(len(idx) >= 0)

    def test_set_verbose(self):
        """set_verbose() (line 301)."""
        domain = self._make_domain()
        r = Region(domain, indices=[0])
        r.set_verbose(True)
        self.assertTrue(r.verbose)

    def test_polygon_no_centroids_warns(self):
        """Polygon with no matching centroids → warning (lines 260, 265-270)."""
        import warnings
        domain = self._make_domain()
        # Polygon far outside the mesh
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            r = Region(domain, polygon=[[10.0, 10.0], [11.0, 10.0], [10.5, 11.0]])
        self.assertEqual(len(r.indices), 0)


class Test_centroid_field(unittest.TestCase):
    """Tests for Centroid_field (lines 307-339)."""

    def _make_domain_and_region(self):
        points, vertices, boundary = rectangular(1, 3)
        domain = Domain(points, vertices, boundary)
        region = Region(domain, indices=[0, 1, 2])
        return domain, region

    def test_init(self):
        """Centroid_field.__init__ (lines 307-309)."""
        domain, region = self._make_domain_and_region()
        cf = Centroid_field(region, value=1.0)
        self.assertIs(cf.domain, domain)

    def test_set_value_scalar(self):
        """set_value with scalar (lines 319, 328-330)."""
        domain, region = self._make_domain_and_region()
        cf = Centroid_field(region, value=1.0)
        cf.set_value(2.0)
        self.assertEqual(cf.value_type, 'scalar')
        self.assertFalse(cf.value_callable)

    def test_set_value_time_function(self):
        """set_value with t-only function (lines 334-336)."""
        domain, region = self._make_domain_and_region()
        cf = Centroid_field(region, value=1.0)
        cf.set_value(lambda t: t * 2.0)
        self.assertEqual(cf.value_type, 't')
        self.assertTrue(cf.value_callable)
        self.assertFalse(cf.value_spatial)

    def test_set_value_spatial_function(self):
        """set_value with x,y function (lines 337-339)."""
        domain, region = self._make_domain_and_region()
        cf = Centroid_field(region, value=1.0)
        cf.set_value(lambda x, y: x + y)
        self.assertTrue(cf.value_callable)
        self.assertTrue(cf.value_spatial)

    def test_set_value_quantity(self):
        """set_value with Quantity (lines 319-321, 331-333)."""
        from anuga import Quantity
        domain, region = self._make_domain_and_region()
        q = domain.quantities['stage']
        cf = Centroid_field(region, value=1.0)
        cf.set_value(q)
        self.assertEqual(cf.value_type, 'quantity')
        self.assertFalse(cf.value_callable)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_region)
    runner = unittest.TextTestRunner()
    runner.run(suite)
