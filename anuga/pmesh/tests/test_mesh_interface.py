#!/usr/bin/env python


import tempfile
import unittest
import os
import numpy
import pytest

from anuga.pmesh.mesh import import_mesh_from_file
from anuga.pmesh.mesh_interface import create_pmesh_from_regions
from anuga.pmesh.mesh_interface import _create_pmesh_from_regions, PolygonError
from anuga.geometry.polygon import polygon_outlines_cross
from anuga.load_mesh.loadASCII import *
from anuga.geometry.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference,DEFAULT_ZONE


class TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_create_mesh_from_regions(self):
        x=-500
        y=-1000
        mesh_geo = geo_reference=Geo_reference(56, x, y)

        # These are the absolute values
        polygon_absolute = [[0,0], [100,0], [100,100], [0,100]]

        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner2_polygon_absolute)

        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]

        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.get_user_segments()
        self.assertTrue(len(segs)==12, 'FAILED!')
        self.assertTrue(len(m.userVertices)==12, 'FAILED!')
        self.assertTrue(segs[0].tag=='walls', 'FAILED!')
        self.assertTrue(segs[1].tag=='walls', 'FAILED!')
        self.assertTrue(segs[2].tag=='bom', 'FAILED!')
        self.assertTrue(segs[3].tag=='bom', 'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[0]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        msg = ('Expected point (%s,%s) to be inside polygon %s'
               % (str(poly_point.x+x), str(poly_point.y+y),
                  str(polygon_absolute)))
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          polygon_absolute, closed=False),
                        msg)

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[1]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[2]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner2_polygon_absolute,
                                          closed=False),
                        'FAILED!')


    def test_create_mesh_from_regions_with_caching(self):
        x=-500
        y=-1000
        mesh_geo = geo_reference=Geo_reference(56, x, y)

        # These are the absolute values
        polygon_absolute = [[0,0], [100,0], [100,100], [0,100]]

        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                             change_points_geo_ref(inner2_polygon_absolute)

        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]

        interior_holes = None

        # Clear cache first
        from anuga.caching import cache

        cache(_create_pmesh_from_regions,
              (polygon, boundary_tags),
              {'minimum_triangle_angle': 28.0,
               'maximum_triangle_area': 10000000,
               'interior_regions': interior_regions,
               'interior_holes': interior_holes,
               'poly_geo_reference': geo_ref_poly,
               'mesh_geo_reference': mesh_geo,
               'verbose': False},
              verbose=False,
              clear=1)

        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     maximum_triangle_area=10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo,
                                     verbose=False,
                                     use_cache=True)


        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.get_user_segments()
        self.assertTrue(len(segs)==12, 'FAILED!')
        self.assertTrue(len(m.userVertices)==12, 'FAILED!')
        self.assertTrue(segs[0].tag=='walls', 'FAILED!')
        self.assertTrue(segs[1].tag=='walls', 'FAILED!')
        self.assertTrue(segs[2].tag=='bom', 'FAILED!')
        self.assertTrue(segs[3].tag=='bom', 'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[0]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[1]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_pmesh_from_regions
        # a black box)
        poly_point = m.get_regions()[2]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner2_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Now create m using cached values
        m_cache = create_pmesh_from_regions(polygon,
                                           boundary_tags,
                                           10000000,
                                           interior_regions=interior_regions,
                                           poly_geo_reference=geo_ref_poly,
                                           mesh_geo_reference=mesh_geo,
                                           verbose=False,
                                           use_cache=True)


    def test_create_mesh_from_regions2(self):
        # These are the absolute values
        min_x = -10
        min_y = -88
        polygon_absolute = [[min_x,min_y], [1000,100], [1000,1000], [100,1000]]

        x_p = -10
        y_p = -40
        zone = 23
        geo_ref_poly = Geo_reference(zone, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner2_polygon_absolute)

        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly)

        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.get_user_segments()
        self.assertTrue(len(segs)==12, 'FAILED!')
        self.assertTrue(len(m.userVertices)==12, 'FAILED!')
        self.assertTrue(segs[0].tag=='walls', 'FAILED!')
        self.assertTrue(segs[1].tag=='walls', 'FAILED!')
        self.assertTrue(segs[2].tag=='bom', 'FAILED!')
        self.assertTrue(segs[3].tag=='bom', 'FAILED!')
        self.assertTrue(m.geo_reference.get_zone()==zone, 'FAILED!')
        self.assertTrue(m.geo_reference.get_xllcorner()==min_x, 'FAILED!')
        self.assertTrue(m.geo_reference.get_yllcorner()==min_y, 'FAILED!')


    def test_create_mesh_from_regions3(self):
        # These are the absolute values
        min_x = -10
        min_y = -88
        polygon = [[min_x,min_y], [1000,100], [1000,1000], [100,1000]]


        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner2_polygon_absolute)

        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        # Test the mesh instance
        self.assertTrue(len(m.regions) == 3, 'FAILED!')
        segs = m.get_user_segments()
        self.assertTrue(len(segs) == 12, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 12, 'FAILED!')
        self.assertTrue(segs[0].tag == 'walls', 'FAILED!')
        self.assertTrue(segs[1].tag == 'walls', 'FAILED!')
        self.assertTrue(segs[2].tag == 'bom', 'FAILED!')
        self.assertTrue(segs[3].tag == 'bom', 'FAILED!')
        self.assertTrue(m.geo_reference.get_zone() == DEFAULT_ZONE, 'FAILED!')
        self.assertTrue(m.geo_reference.get_xllcorner() == min_x, 'FAILED!')
        self.assertTrue(m.geo_reference.get_yllcorner() == min_y, 'FAILED!')


    def test_create_mesh_from_regions4(self):
        fd, file_name = tempfile.mkstemp('.tsh')
        os.close(fd)

        # These are the absolute values
        density_outer = 1000
        min_outer = 0
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer], [max_outer,min_outer],
                         [max_outer,max_outer], [min_outer,max_outer]]

        density_inner1 = 10000000
        inner_buffer = 100
        min_inner1 = min_outer + inner_buffer
        max_inner1 = max_outer - inner_buffer
        inner1_polygon = [[min_inner1,min_inner1], [max_inner1,min_inner1],
                          [max_inner1,max_inner1], [min_inner1,max_inner1]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        interior_regions = [(inner1_polygon, density_inner1)]
        create_pmesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name, verbose=False)

        m = import_mesh_from_file(file_name)

        self.assertTrue(len(m.get_triangulation()) <= 900,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.get_triangulation()) >= 200,
                        'Test mesh interface failed!')

        create_pmesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = import_mesh_from_file(file_name)

        self.assertTrue(len(m.get_triangulation()) <= 100,
                        'Test mesh interface failed!')

        os.remove(file_name)


    def test_create_mesh_from_regions5(self):
        fd, file_name = tempfile.mkstemp('.tsh')
        os.close(fd)

        # These are the absolute values
        density_outer = 10000000
        min_outer = 0
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer], [max_outer,min_outer],
                         [max_outer,max_outer], [min_outer,max_outer]]

        density_inner1 = 1000
        inner_buffer = 100
        min_inner1 = min_outer + inner_buffer
        max_inner1 = max_outer - inner_buffer
        inner1_polygon = [[min_inner1,min_inner1], [max_inner1,min_inner1],
                          [max_inner1,max_inner1], [min_inner1,max_inner1]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        interior_regions = [(inner1_polygon, density_inner1)]
        create_pmesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = import_mesh_from_file(file_name)
        self.assertTrue(len(m.get_triangulation()) <= 2000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.get_triangulation()) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)


    def test_create_mesh_from_regions6(self):
        fd, file_name = tempfile.mkstemp('.tsh')
        os.close(fd)

        # These are the absolute values
        density_outer = 1000
        min_outer = 0
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer], [max_outer,min_outer],
                         [max_outer,max_outer], [min_outer,max_outer]]

        delta = 10
        density_inner1 = 1000
        min_inner1 = min_outer + delta
        max_inner1 = max_outer - delta
        inner1_polygon = [[min_inner1,min_inner1], [max_inner1,min_inner1],
                          [max_inner1,max_inner1], [min_inner1,max_inner1]]

        density_inner2 = 10000000
        min_inner2 = min_outer +  2*delta
        max_inner2 = max_outer -  2*delta
        inner2_polygon = [[min_inner2,min_inner2], [max_inner2,min_inner2],
                          [max_inner2,max_inner2], [min_inner2,max_inner2]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        interior_regions = [(inner1_polygon, density_inner1),
                            (inner2_polygon, density_inner2)]
        create_pmesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = import_mesh_from_file(file_name)
        self.assertTrue(len(m.get_triangulation()) <= 2000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.get_triangulation()) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)


    def test_create_mesh_from_regions7(self):
        fd, file_name = tempfile.mkstemp('.tsh')
        os.close(fd)

        # These are the absolute values
        density_outer = 1001
        min_outer = 0
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer], [max_outer,min_outer],
                         [max_outer,max_outer], [min_outer,max_outer]]

        delta = 10
        density_inner1 = 100000000
        min_inner1 = min_outer + delta
        max_inner1 = max_outer - delta
        inner1_polygon = [[min_inner1,min_inner1], [max_inner1,min_inner1],
                          [max_inner1,max_inner1], [min_inner1,max_inner1]]

        density_inner2 = 1000
        min_inner2 = min_outer +  2*delta
        max_inner2 = max_outer -  2*delta
        inner2_polygon = [[min_inner2,min_inner2], [max_inner2,min_inner2],
                          [max_inner2,max_inner2], [min_inner2,max_inner2]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        # Note the list order is important
        # The last region added will be the region triangle uses,
        # if two regions points are in the same bounded area.
        interior_regions = [(inner2_polygon, density_inner2),
                            (inner1_polygon, density_inner1)]
        create_pmesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = import_mesh_from_file(file_name)
        self.assertTrue(len(m.get_triangulation()) <= 3000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.get_triangulation()) >= 2000,
                        'Test mesh interface failed!')

        os.remove(file_name)


    def test_create_mesh_from_regions_interior_regions(self):
        '''Test that create_pmesh_from_regions fails when an interior
        region is outside bounding polygon.
        '''

        # These are the absolute values
        min_x = 10
        min_y = 88
        polygon = [[min_x,min_y], [1000,100], [1000,1000], [100,1000]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400], [900,500], [800,600]]

        interior_regions = [(inner_polygon, 5)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        # This one sticks outside bounding polygon - should fail
        inner_polygon = [[800,400], [900,500], [800,600], [200, 995]]
        inner_polygon1 = [[800,400], [1100,500], [800,600]]
        interior_regions = [[inner_polygon, 50], [inner_polygon1, 50]]

        try:
            m = create_pmesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions,
                                         verbose=False)
        except Exception:
            pass
        else:
            msg = 'Interior polygon sticking outside bounding polygon should '
            msg += 'cause an Exception to be raised'
            raise Exception(msg)


    def test_create_mesh_from_regions_interior_regions1(self):
        '''Test that create_pmesh_from_regions fails
        when an interior region is outside bounding polygon.
        '''

        # These are the values
        d0 = [310000, 7690000]
        d1 = [280000, 7690000]
        d2 = [270000, 7645000]
        d3 = [240000, 7625000]
        d4 = [270000, 7580000]
        d5 = [300000, 7590000]
        d6 = [340000, 7610000]
        poly_all = [d0, d1, d2, d3, d4, d5, d6]

        i0 = [304000, 7607000]
        i1 = [302000, 7605000]
        i2 = [304000, 7603000]
        i3 = [307000, 7602000]
        i4 = [309000, 7603000]
#        i4 = [310000, 7580000]
        i5 = [307000, 7606000]
        poly_onslow = [i0, i1, i2, i3, i4, i5]

        # Thevenard Island
        j0 = [294000, 7629000]
        j1 = [285000, 7625000]
        j2 = [294000, 7621000]
        j3 = [299000, 7625000]
        poly_thevenard = [j0, j1, j2, j3]

        # med res around onslow
        l0 = [300000, 7610000]
        l1 = [285000, 7600000]
        l2 = [300000, 7597500]
        l3 = [310000, 7770000] # this one is outside
#        l3 = [310000, 7630000] # this one is NOT outside
        l4 = [315000, 7610000]
        poly_coast = [l0, l1, l2, l3, l4]

        # general coast and local area to onslow region
        m0 = [270000, 7581000]
        m1 = [300000, 7591000]
        m2 = [339000, 7610000]
        m3 = [330000, 7630000]
        m4 = [290000, 7640000]
        m5 = [260000, 7600000]
        poly_region = [m0, m1, m2, m3, m4, m5]

        # This one sticks outside bounding polygon - should fail
        interior_regions = [[poly_onslow, 50000], [poly_region, 50000],
                            [poly_coast, 100000], [poly_thevenard, 100000]]
        boundary_tags = {'walls': [0,1], 'bom': [2], 'exterior': [3,4,5,6]}

        try:
            m = create_pmesh_from_regions(poly_all,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions,
                                         verbose=False)
        except Exception:
            pass
        else:
            msg = 'Interior polygon sticking outside bounding polygon should '
            msg += 'cause an Exception to be raised'
            raise Exception(msg)


    def FIXMEtest_create_mesh_with_multiply_tagged_segments(self):
        '''Test that create_pmesh_from_regions fails when
        segments are listed repeatedly in boundary_tags.
        '''

        # These are the absolute values
        min_x = 10
        min_y = 88
        polygon = [[min_x,min_y], [1000,100], [1000,1000], [100,1000]]
        boundary_tags = {'walls': [0,1], 'bom': [1,2]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400], [900,500], [800,600]]
        interior_regions = [(inner_polygon, 5)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     verbose=False)

        # This one sticks outside bounding polygon - should fail
        inner_polygon = [[800,400], [900,500], [800,600]]
        interior_regions = [(inner_polygon, 5)]

        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)
        try:
            m = create_pmesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except Exception:
            pass
        else:
            msg = 'Tags are listed repeatedly, but create mesh from regions '
            msg += 'does not cause an Exception to be raised'
            raise Exception(msg)


    def test_create_mesh_with_segments_out_of_bounds(self):
        """Test that create_pmesh_from_regions fails when a segment is out of bounds.
        """

        # These are the absolute values
        min_x = 10
        min_y = 88
        polygon = [[min_x,min_y],[1000,100],[1000,1000],[100,1000]]


        boundary_tags = {'walls': [0,1], 'bom':[2,3], 'out': [5]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400],[900,500],[800,600]]

        interior_regions = [(inner_polygon, 5)]


        try:
            m = create_pmesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except Exception:
            pass
        else:
            msg = 'Tags are listed repeatedly, but create mesh from regions '
            msg += 'does not cause an Exception to be raised'
            raise Exception(msg)


    def test_create_mesh_with_breaklines(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     breaklines=[[[50,50],[2000,2000]]])

        self.assertTrue(len(m.regions) == 1, 'FAILED!')
        segs = m.get_user_segments()
        self.assertTrue(len(segs) == 5, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 6, 'FAILED!')


    def test_create_mesh_with_interior_holes(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]

        interior_poly1 = [[101,101],[200,200], [101,200]]
        interior_poly2 = [[300,300],[500,500], [400,200]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}


        # This should work with one hole
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1])

        self.assertTrue(len(m.get_user_segments()) == 7, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 7, 'FAILED!')


        # This should work with two holes
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1, interior_poly2])
        #print len(m.get_user_segments())
        #print len(m.userVertices)

        self.assertTrue(len(m.get_user_segments()) == 10, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 10, 'FAILED!')

        #-------------------------------------
        # Error testing
        #-------------------------------------


        # try passing just one polygon, not a list of polygons, should throw exception
        try:
            m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=interior_poly1)
        except Exception:
            pass
        else:
            msg = 'Passing a single polygon should have raised an error '
            raise Exception(msg)



        # interior polygon outside bounding polygon, should throw exception
        interior_poly3  = [[50,50],[500,500], [400,200]]
        try:
            m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly3])
        except Exception:
            pass
        else:
            msg = 'Passing a single polygon should have raised an error '
            raise Exception(msg)


    def test_create_mesh_with_interior_holes_and_tags(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]

        interior_poly1 = [[101,101],[200,200], [101,200]]
        interior_poly2 = [[300,300],[500,500], [400,200]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}


        # This should work with one hole
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1],
                                     hole_tags=[{'edge0' : [0], 'edge1': [1], 'edge2': [2]}])

        self.assertTrue(len(m.get_user_segments()) == 7, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 7, 'FAILED!')

        m.generate_mesh()

        tags_list = m.get_mesh_segment_tags()

        assert len([x for x in tags_list if x == 'edge0']) == 7
        assert len([x for x in tags_list if x == 'edge1']) == 6
        assert len([x for x in tags_list if x == 'edge2']) == 54


    def test_create_mesh_from_regions_with_duplicate_verts(self):
        # These are the absolute values
        polygon_absolute = [[0.0, 0.0],
                            [0, 4.0],
                            [4.0, 4.0],
                            [4.0, 0.0],
                            [4.0, 0.0]]
        x_p = -10
        y_p = -40
        zone = 23
        geo_ref_poly = Geo_reference(zone, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)
        boundary_tags = {'50': [0],
                         '40': [1],
                         '30': [2],
                         'no where seg': [3],
                         '20': [4]}
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     poly_geo_reference=geo_ref_poly,
                                     verbose=False)

        fileName = 'badmesh.tsh'
        #m.export_mesh_file(fileName)

    def concept_create_mesh_from_regions_with_ungenerate(self):
        x=0
        y=0
        mesh_geo = geo_reference=Geo_reference(56, x, y)

        # These are the absolute values
        polygon_absolute = [[0,0], [100,0], [100,100], [0,100]]
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls': [0,1], 'bom': [2]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner2_polygon_absolute)

        max_area = 10000000
        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        m.export_mesh_file('a_test_mesh_iknterface.tsh')

        fd, fileName = tempfile.mkstemp('.txt')
        os.close(fd)
        file = open(fileName, 'w')
        file.write('         1       ??      ??\n\
       90.0       90.0\n\
       81.0       90.0\n\
       81.0       81.0\n\
       90.0       81.0\n\
       90.0       90.0\n\
END\n\
         2      ?? ??\n\
       10.0       80.0\n\
       10.0       90.0\n\
       20.0       90.0\n\
       10.0       80.0\n\
END\n\
END\n')
        file.close()

        m.import_ungenerate_file(fileName, tag='wall')
        os.remove(fileName)
        m.generate_mesh(maximum_triangle_area=max_area, verbose=False)
        m.export_mesh_file('b_test_mesh_iknterface.tsh')

    def concept_ungenerateII(self):
        from anuga import Domain, Reflective_boundary, Dirichlet_boundary

        x=0
        y=0
        mesh_geo = geo_reference=Geo_reference(56, x, y)

        # These are the absolute values
        polygon_absolute = [[0,0], [100,0], [100,100], [0,100]]
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'wall': [0,1,3], 'wave': [2]}

        inner1_polygon_absolute = [[10,10], [20,10], [20,20], [10,20]]
        inner1_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30], [40,30], [40,40], [30,40]]
        inner2_polygon = geo_ref_poly.\
                            change_points_geo_ref(inner2_polygon_absolute)

        max_area = 1
        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        m.export_mesh_file('a_test_mesh_iknterface.tsh')

        fd, fileName = tempfile.mkstemp('.txt')
        os.close(fd)
        file = open(fileName, 'w')
        file.write('         1       ??      ??\n\
       90.0       90.0\n\
       81.0       90.0\n\
       81.0       81.0\n\
       90.0       81.0\n\
       90.0       90.0\n\
END\n\
         2      ?? ??\n\
       10.0       80.0\n\
       10.0       90.0\n\
       20.0       90.0\n\
       10.0       80.0\n\
END\n\
END\n')
        file.close()

        m.import_ungenerate_file(fileName)      #, tag='wall')
        os.remove(fileName)
        m.generate_mesh(maximum_triangle_area=max_area, verbose=False)
        mesh_filename = 'bento_b.tsh'
        m.export_mesh_file(mesh_filename)

        domain = Domain(mesh_filename, use_cache = False)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([3, 0, 0])
        domain.set_boundary({'wall': Br, 'wave': Bd})
        yieldstep = 0.1
        finaltime = 10
        for t in domain.evolve(yieldstep, finaltime):
            domain.write_time()

    def concept_ungenerateIII(self):
        from anuga import Domain, Reflective_boundary, \
                            Dirichlet_boundary
        from anuga.pmesh.mesh_interface import create_pmesh_from_regions

        # These are the absolute values
        polygon = [[0,0], [100,0], [100,100], [0,100]]

        boundary_tags = {'wall': [0,1,3], 'wave': [2]}
        inner1_polygon = [[10,10], [20,10], [20,20], [10,20]]
        inner2_polygon = [[30,30], [40,30], [40,40], [30,40]]

        max_area = 1
        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions)

        fd, fileName = tempfile.mkstemp('.txt')
        os.close(fd)
        file = open(fileName, 'w')
        file.write('         1       ??      ??\n\
       90.0       90.0\n\
       81.0       90.0\n\
       81.0       81.0\n\
       90.0       81.0\n\
       90.0       90.0\n\
END\n\
         2      ?? ??\n\
       10.0       80.0\n\
       10.0       90.0\n\
       20.0       90.0\n\
       10.0       80.0\n\
END\n\
END\n')
        file.close()

        m.import_ungenerate_file(fileName)
        os.remove(fileName)
        m.generate_mesh(maximum_triangle_area=max_area, verbose=False)
        mesh_filename = 'mesh.tsh'
        m.export_mesh_file(mesh_filename)

        domain = Domain(mesh_filename, use_cache=False)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([3, 0, 0])
        domain.set_boundary({'wall': Br, 'wave': Bd})
        yieldstep = 0.1
        finaltime = 10
        for t in domain.evolve(yieldstep, finaltime):
            domain.write_time()

    def test_create_mesh_from_regions_check_segs(self):
        '''Test that create_pmesh_from_regions fails
        when an interior region is outside bounding polygon.
        '''

        # These are the absolute values
        min_x = 10
        min_y = 88
        polygon = [[min_x,min_y], [1000,100], [1000,1000], [100,1000]]
        boundary_tags = {'walls': [0,1,3], 'bom': [2]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400], [900,500], [800,600]]
        interior_regions = [(inner_polygon, 5)]
        m = create_pmesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        boundary_tags = {'walls': [0,1,3,4], 'bom': [2]}
        try:
            m = create_pmesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except Exception:
            pass
        else:
            msg = 'Segment out of bounds not caught '
            raise Exception(msg)

class TestCrossingOutlines:
    """Interior region and hole outlines must not cross (issue #497)"""

    bounding = [[0., 0.], [20., 0.], [20., 10.], [0., 10.]]
    town = [[8., 3.], [12., 3.], [12., 7.], [8., 7.]]
    crossing = [[10., 5.], [15., 5.], [15., 8.], [10., 8.]]
    nested = [[9., 4.], [11., 4.], [11., 6.], [9., 6.]]
    separate = [[14., 3.], [16., 3.], [16., 7.], [14., 7.]]
    shares_edge = [[12., 3.], [15., 3.], [15., 7.], [12., 7.]]

    def mesh(self, regions=None, holes=None):
        return _create_pmesh_from_regions(self.bounding,
                                          {'all': [0, 1, 2, 3]},
                                          maximum_triangle_area=1.0,
                                          interior_regions=regions,
                                          interior_holes=holes,
                                          verbose=False)

    def test_crossing_regions_raise_in_either_order(self):
        for a, b in ((self.town, self.crossing), (self.crossing, self.town)):
            with pytest.raises(PolygonError, match='size_field_from_regions'):
                self.mesh(regions=[(a, 0.1), (b, 0.5)])

    def test_message_names_the_two_polygons(self):
        with pytest.raises(PolygonError, match=r'interior_regions\[0\] and interior_regions\[2\]'):
            self.mesh(regions=[(self.town, 0.1), (self.separate, 0.1), (self.crossing, 0.5)])

    def test_crossing_with_no_vertex_inside_the_other_raises(self):
        wide = [[6., 4.], [14., 4.], [14., 6.], [6., 6.]]
        tall = [[9., 1.], [11., 1.], [11., 9.], [9., 9.]]
        with pytest.raises(PolygonError):
            self.mesh(regions=[(wide, 0.1), (tall, 0.5)])

    def test_edge_through_a_vertex_raises(self):
        # Two edges enter town exactly through its corners (12, 3) and (12, 7)
        through = [[13., 2.], [10., 5.], [13., 8.]]
        assert polygon_outlines_cross(self.town, through)
        assert polygon_outlines_cross(through, self.town)

    def test_inside_and_touching_at_corners_does_not_cross(self):
        diamond = [[10., 1.], [14., 5.], [10., 9.], [6., 5.]]
        assert not polygon_outlines_cross(self.town, diamond)

    def test_crossing_holes_raise(self):
        with pytest.raises(PolygonError, match='interior_holes'):
            self.mesh(holes=[self.town, self.crossing])

    def test_hole_crossing_a_region_raises(self):
        with pytest.raises(PolygonError, match='interior_regions.*interior_holes'):
            self.mesh(regions=[(self.town, 0.1)], holes=[self.crossing])

    def test_polygons_that_do_not_cross_are_accepted(self):
        self.mesh(regions=[(self.town, 0.1), (self.nested, 0.05)])
        self.mesh(regions=[(self.nested, 0.05), (self.town, 0.1)])
        self.mesh(regions=[(self.town, 0.1), (self.separate, 0.1)])
        self.mesh(regions=[(self.town, 0.1), (self.shares_edge, 0.1)])
        self.mesh(regions=[(self.town, 0.1)], holes=[self.separate])
        self.mesh(regions=[(self.town, 0.1)], holes=[self.nested])
        self.mesh(regions=[(self.nested, 0.1)], holes=[self.town])

    def test_touching_outlines_do_not_cross(self):
        corner = [[12., 7.], [15., 7.], [15., 9.], [12., 9.]]
        part_of_edge = [[12., 4.], [15., 4.], [15., 6.], [12., 6.]]
        for other in (self.shares_edge, corner, part_of_edge):
            assert not polygon_outlines_cross(self.town, other)
            assert not polygon_outlines_cross(other, self.town)

    def test_large_utm_coordinates(self):
        shift = numpy.array([500000., 6200000.])
        town = numpy.array(self.town) + shift
        assert polygon_outlines_cross(town, numpy.array(self.crossing) + shift)
        assert not polygon_outlines_cross(town, numpy.array(self.shares_edge) + shift)
        assert not polygon_outlines_cross(town, numpy.array(self.nested) + shift)


################################################################################


class TestNestedSeedPoints:
    """A region's or hole's seed point must not fall inside a polygon nested
    in it, or the nested polygon takes its area (issue #504)"""

    bounding = [[0., 0.], [20., 0.], [20., 10.], [0., 10.]]
    tags = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}
    outer = [[8., 3.], [12., 3.], [12., 7.], [8., 7.]]
    # Includes the first vertex of outer
    inner = [[8., 3.], [10., 3.], [10., 5.], [8., 5.]]
    # Includes the first vertex of the bounding polygon
    corner = [[0., 0.], [3., 0.], [3., 3.], [0., 3.]]
    small_hole = [[8., 3.], [9., 3.], [9., 4.], [8., 4.]]

    def domain(self, shift=(0., 0.), **kwargs):
        import anuga
        def move(polygon):
            return [[x + shift[0], y + shift[1]] for x, y in polygon]
        for key in ('interior_regions',):
            if key in kwargs:
                kwargs[key] = [(move(p), r) for p, r in kwargs[key]]
        if 'interior_holes' in kwargs:
            kwargs['interior_holes'] = [move(p) for p in kwargs['interior_holes']]
        domain = anuga.create_domain_from_regions(move(self.bounding), self.tags,
                                                  verbose=False, **kwargs)
        domain.test_shift = numpy.array(shift)
        return domain

    def largest(self, domain, polygon, without=()):
        """Largest triangle with its centroid in polygon but not in without"""
        from anuga.geometry.polygon import inside_polygon
        c = domain.get_centroid_coordinates(absolute=True)
        shift = domain.test_shift
        inside = numpy.zeros(len(c), bool)
        inside[inside_polygon(c, numpy.array(polygon) + shift)] = True
        for other in without:
            inside[inside_polygon(c, numpy.array(other) + shift)] = False
        assert inside.any()
        return domain.areas[inside].max()

    def test_inner_region_at_outer_first_vertex(self):
        for regions in ([(self.outer, 0.05), (self.inner, 0.01)],
                        [(self.inner, 0.01), (self.outer, 0.05)]):
            domain = self.domain(maximum_triangle_area=0.2,
                                 interior_regions=regions)
            assert self.largest(domain, self.outer, [self.inner]) <= 0.05
            assert self.largest(domain, self.inner) <= 0.01

    def test_region_at_bounding_first_vertex(self):
        domain = self.domain(maximum_triangle_area=0.2,
                             interior_regions=[(self.corner, 5.0)])
        assert self.largest(domain, self.bounding, [self.corner]) <= 0.2

    def test_hole_at_region_first_vertex(self):
        domain = self.domain(maximum_triangle_area=0.2,
                             interior_regions=[(self.outer, 0.05)],
                             interior_holes=[self.small_hole])
        assert self.largest(domain, self.outer) <= 0.05
        assert numpy.allclose(domain.areas.sum(), 199.0)

    def test_hole_at_bounding_first_vertex(self):
        domain = self.domain(maximum_triangle_area=0.2,
                             interior_holes=[self.corner])
        assert domain.areas.max() <= 0.2
        assert numpy.allclose(domain.areas.sum(), 191.0)

    def test_utm_coordinates(self):
        domain = self.domain(shift=(300000., 6180000.),
                             maximum_triangle_area=0.2,
                             interior_regions=[(self.outer, 0.05),
                                               (self.inner, 0.01)])
        assert self.largest(domain, self.outer, [self.inner]) <= 0.05
        assert self.largest(domain, self.bounding, [self.outer]) <= 0.2

    def test_geo_referenced_polygons(self):
        # Polygons given relative to poly_geo_reference
        geo = Geo_reference(zone=56, xllcorner=300000., yllcorner=6180000.)
        m = _create_pmesh_from_regions(self.bounding, self.tags,
                                       maximum_triangle_area=0.2,
                                       interior_regions=[(self.outer, 0.05),
                                                         (self.inner, 0.01)],
                                       poly_geo_reference=geo,
                                       verbose=False)
        seeds = [(r.x, r.y, r.get_max_area()) for r in m.regions
                 if r.get_max_area() is not None]
        assert sorted(area for x, y, area in seeds) == [0.01, 0.05, 0.2]
        outer = numpy.array(self.outer)
        inner = numpy.array(self.inner)
        for x, y, area in seeds:
            # Seeds are relative to the mesh geo_reference
            [[x, y]] = m.geo_reference.get_absolute([[x, y]])
            x -= 300000.
            y -= 6180000.
            if area == 0.05:
                assert is_inside_polygon([x, y], outer, closed=False)
                assert not is_inside_polygon([x, y], inner)


class TestBreaklineAreaLimit(unittest.TestCase):
    """#511: breaklines from boundary to boundary used to drop the area limit
    beyond them."""

    @staticmethod
    def _cell_areas(m):
        pts = numpy.array(m.get_mesh_vertices())
        P = pts[numpy.array(m.get_triangulation())]
        area = 0.5 * numpy.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
                               - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
        return P.mean(axis=1), area

    _BOX = [[0.0, 0.0], [100.0, 0.0], [100.0, 30.0], [0.0, 30.0]]
    _TAGS = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}

    def test_breaklines_across_domain_keep_maximum_area(self):
        breaklines = [[[0.0, 10.0], [100.0, 10.0]], [[0.0, 20.0], [100.0, 20.0]]]
        m = create_pmesh_from_regions(self._BOX, self._TAGS, maximum_triangle_area=1.0,
                                      breaklines=breaklines)
        m.generate_mesh(verbose=False)
        c, area = self._cell_areas(m)
        for lo, hi in ((0, 10), (10, 20), (20, 30)):
            strip = (c[:, 1] > lo) & (c[:, 1] < hi)
            assert strip.sum() > 0
            assert area[strip].max() <= 1.0 + 1e-9, (lo, hi, area[strip].max())

    def test_coarse_interior_region_stays_coarse(self):
        # The global limit must not refine a region coarser than the default
        coarse = [[[50.0, 0.0], [100.0, 0.0], [100.0, 30.0], [50.0, 30.0]], 25.0]
        m = create_pmesh_from_regions(self._BOX, self._TAGS, maximum_triangle_area=1.0,
                                      interior_regions=[coarse])
        m.generate_mesh(verbose=False)
        c, area = self._cell_areas(m)
        assert area[c[:, 0] < 50].max() <= 1.0 + 1e-9
        assert 5.0 < area[c[:, 0] > 50].max() <= 25.0 + 1e-9
        assert m.default_maximum_triangle_area == 25.0

    def test_area_limit_in_exponent_notation(self):
        # str(1e20) is '1e+20', and Triangle read 'a1e+20' as a limit of 1
        # (the weir_1 validation mesh went from 2 219 to 16 129 triangles).
        breaklines = [[[50.0, 0.0], [50.0, 30.0]]]
        m = create_pmesh_from_regions(self._BOX, self._TAGS, maximum_triangle_area=1.0e20,
                                      breaklines=breaklines,
                                      regionPtArea=[[25.0, 15.0, 50.0], [75.0, 15.0, 50.0]])
        m.generate_mesh(verbose=False)
        c, area = self._cell_areas(m)
        assert 10.0 < area.max() <= 50.0 + 1e-9, area.max()

    def test_breaklines_across_domain_via_create_domain(self):
        import anuga
        domain = anuga.create_domain_from_regions(
            self._BOX, self._TAGS, maximum_triangle_area=1.0,
            breaklines=[[[0.0, 10.0], [100.0, 10.0]]])
        assert domain.areas.max() <= 1.0 + 1e-9


class TestCreateDomainMeshFilename:
    """mesh_filename writes the mesh the domain is built on (#506)."""

    _BOX = [[0.0, 0.0], [100.0, 0.0], [100.0, 30.0], [0.0, 30.0]]
    _TAGS = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}

    def _check_file_matches(self, domain, filename):
        import anuga
        assert os.path.exists(filename)
        from_file = anuga.create_domain_from_file(filename)
        assert from_file.number_of_triangles == domain.number_of_triangles
        numpy.testing.assert_allclose(
            from_file.get_vertex_coordinates(absolute=True),
            domain.get_vertex_coordinates(absolute=True))
        assert (sorted(from_file.get_boundary_tags())
                == sorted(domain.get_boundary_tags()))

    @pytest.mark.parametrize('ext', ['.msh', '.tsh'])
    def test_regions_writes_mesh(self, tmp_path, ext):
        import anuga
        filename = str(tmp_path / ('mesh' + ext))
        inner = [[40.0, 10.0], [60.0, 10.0], [60.0, 20.0], [40.0, 20.0]]
        domain = anuga.create_domain_from_regions(
            self._BOX, self._TAGS, maximum_triangle_area=20.0,
            interior_regions=[[inner, 2.0]], mesh_filename=filename)
        self._check_file_matches(domain, filename)

    def test_regions_writes_relaxed_mesh(self, tmp_path):
        # The file is written after relaxation, so it is the domain's mesh.
        import anuga
        filename = str(tmp_path / 'mesh.msh')
        hole = [[40.0, 10.0], [60.0, 10.0], [60.0, 20.0], [40.0, 20.0]]
        domain = anuga.create_domain_from_regions(
            self._BOX, self._TAGS, maximum_triangle_area=20.0,
            interior_holes=[hole], relax_passes=5, mesh_filename=filename)
        self._check_file_matches(domain, filename)

    def test_regions_writes_nothing_by_default(self, tmp_path, monkeypatch):
        import anuga
        monkeypatch.chdir(tmp_path)
        anuga.create_domain_from_regions(self._BOX, self._TAGS,
                                         maximum_triangle_area=20.0)
        assert os.listdir(tmp_path) == []

    def test_dem_writes_mesh(self, tmp_path):
        import anuga
        filename = str(tmp_path / 'mesh.msh')
        x = numpy.linspace(0.0, 100.0, 51)
        y = numpy.linspace(0.0, 30.0, 16)
        Z = 0.01 * x[:, None] + 0.05 * numpy.abs(y[None, :] - 15.0)
        domain = anuga.create_domain_from_dem(
            (x, y, Z), maximum_triangle_area=20.0, mesh_filename=filename)
        self._check_file_matches(domain, filename)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCase)
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)

