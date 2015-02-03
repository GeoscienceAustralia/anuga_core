#!/usr/bin/env python

import tempfile
import unittest
import os

from anuga.pmesh.mesh import importMeshFromFile
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.pmesh.mesh_interface import _create_mesh_from_regions
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

        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==12, 'FAILED!')
        self.assertTrue(len(m.userVertices)==12, 'FAILED!')
        self.assertTrue(segs[0].tag=='walls', 'FAILED!')
        self.assertTrue(segs[1].tag=='walls', 'FAILED!')
        self.assertTrue(segs[2].tag=='bom', 'FAILED!')
        self.assertTrue(segs[3].tag=='bom', 'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[0]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        msg = ('Expected point (%s,%s) to be inside polygon %s'
               % (str(poly_point.x+x), str(poly_point.y+y),
                  str(polygon_absolute)))
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          polygon_absolute, closed=False),
                        msg)

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[1]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[2]

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

        cache(_create_mesh_from_regions,
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

        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     maximum_triangle_area=10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo,
                                     verbose=False,
                                     use_cache=True)


        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs)==12, 'FAILED!')
        self.assertTrue(len(m.userVertices)==12, 'FAILED!')
        self.assertTrue(segs[0].tag=='walls', 'FAILED!')
        self.assertTrue(segs[1].tag=='walls', 'FAILED!')
        self.assertTrue(segs[2].tag=='bom', 'FAILED!')
        self.assertTrue(segs[3].tag=='bom', 'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[0]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[1]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[2]

        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.assertTrue(is_inside_polygon([poly_point.x+x, poly_point.y+y],
                                          inner2_polygon_absolute,
                                          closed=False),
                        'FAILED!')

        # Now create m using cached values
        m_cache = create_mesh_from_regions(polygon,
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
        zone = 808
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly)

        # Test the mesh instance
        self.assertTrue(len(m.regions)==3, 'FAILED!')
        segs = m.getUserSegments()
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        # Test the mesh instance
        self.assertTrue(len(m.regions) == 3, 'FAILED!')
        segs = m.getUserSegments()
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
        file_name = tempfile.mktemp('.tsh')

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
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name, verbose=False)

        m = importMeshFromFile(file_name)

        self.assertTrue(len(m.getTriangulation()) <= 900,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.getTriangulation()) >= 200,
                        'Test mesh interface failed!')

        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = importMeshFromFile(file_name)

        self.assertTrue(len(m.getTriangulation()) <= 100,
                        'Test mesh interface failed!')

        os.remove(file_name)

    def test_create_mesh_from_regions5(self):
        file_name = tempfile.mktemp('.tsh')

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
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = importMeshFromFile(file_name)
        self.assertTrue(len(m.getTriangulation()) <= 2000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.getTriangulation()) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)

    def test_create_mesh_from_regions6(self):
        file_name = tempfile.mktemp('.tsh')

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
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = importMeshFromFile(file_name)
        self.assertTrue(len(m.getTriangulation()) <= 2000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.getTriangulation()) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)

    def test_create_mesh_from_regions7(self):
        file_name = tempfile.mktemp('.tsh')

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
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)

        m = importMeshFromFile(file_name)
        self.assertTrue(len(m.getTriangulation()) <= 3000,
                        'Test mesh interface failed!')
        self.assertTrue(len(m.getTriangulation()) >= 2000,
                        'Test mesh interface failed!')

        os.remove(file_name)

    def test_create_mesh_from_regions_interior_regions(self):
        '''Test that create_mesh_from_regions fails when an interior
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        # This one sticks outside bounding polygon - should fail
        inner_polygon = [[800,400], [900,500], [800,600], [200, 995]]
        inner_polygon1 = [[800,400], [1100,500], [800,600]]
        interior_regions = [[inner_polygon, 50], [inner_polygon1, 50]]

        try:
            m = create_mesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions,
                                         verbose=False)
        except:
            pass
        else:
            msg = 'Interior polygon sticking outside bounding polygon should '
            msg += 'cause an Exception to be raised'
            raise Exception, msg

    def test_create_mesh_from_regions_interior_regions1(self):
        '''Test that create_mesh_from_regions fails
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
        boundary_tags = {'walls': [0,1], 'bom': [2]}

        try:
            m = create_mesh_from_regions(poly_all,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions,
                                         verbose=False)
        except:
            pass
        else:
            msg = 'Interior polygon sticking outside bounding polygon should '
            msg += 'cause an Exception to be raised'
            raise Exception, msg

    def FIXMEtest_create_mesh_with_multiply_tagged_segments(self):
        '''Test that create_mesh_from_regions fails when
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     verbose=False)

        # This one sticks outside bounding polygon - should fail
        inner_polygon = [[800,400], [900,500], [800,600]]
        interior_regions = [(inner_polygon, 5)]

        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)
        try:
            m = create_mesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except:
            pass
        else:
            msg = 'Tags are listed repeatedly, but create mesh from regions '
            msg += 'does not cause an Exception to be raised'
            raise Exception, msg

            
    def test_create_mesh_with_segments_out_of_bounds(self):
        """Test that create_mesh_from_regions fails when a segment is out of bounds.
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
            m = create_mesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except:
            pass
        else:
            msg = 'Tags are listed repeatedly, but create mesh from regions '
            msg += 'does not cause an Exception to be raised'
            raise Exception, msg
            
    def test_create_mesh_with_breaklines(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}

        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     breaklines=[[[50,50],[2000,2000]]])

        self.assertTrue(len(m.regions) == 1, 'FAILED!')
        segs = m.getUserSegments()
        self.assertTrue(len(segs) == 5, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 6, 'FAILED!')


    def test_create_mesh_with_interior_holes(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]
        
        interior_poly1 = [[101,101],[200,200], [101,200]]
        interior_poly2 = [[300,300],[500,500], [400,200]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}


        # This should work with one hole
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1])

        self.assertTrue(len(m.getUserSegments()) == 7, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 7, 'FAILED!')


        # This should work with two holes
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1, interior_poly2])
        #print len(m.getUserSegments())
        #print len(m.userVertices)
        
        self.assertTrue(len(m.getUserSegments()) == 10, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 10, 'FAILED!')

        #-------------------------------------
        # Error testing
        #-------------------------------------

        
        # try passing just one polygon, not a list of polygons, should throw exception
        try:
            m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=interior_poly1)
        except:
            pass
        else:
            msg = 'Passing a single polygon should have raised an error '
            raise Exception, msg



        # interior polygon outside bounding polygon, should throw exception
        interior_poly3  = [[50,50],[500,500], [400,200]]
        try:
            m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly3])
        except:
            pass
        else:
            msg = 'Passing a single polygon should have raised an error '
            raise Exception, msg
        
        
    def test_create_mesh_with_interior_holes_and_tags(self):
        # These are the absolute values
        polygon = [[100,100], [1000,100], [1000,1000], [100,1000]]
        
        interior_poly1 = [[101,101],[200,200], [101,200]]
        interior_poly2 = [[300,300],[500,500], [400,200]]

        boundary_tags = {'walls': [0,1], 'bom': [2,3]}


        # This should work with one hole
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_holes=[interior_poly1],
                                     hole_tags=[{'edge0' : [0], 'edge1': [1], 'edge2': [2]}])

        self.assertTrue(len(m.getUserSegments()) == 7, 'FAILED!')
        self.assertTrue(len(m.userVertices) == 7, 'FAILED!')

        m.generate_mesh()
        
        tags_list = m.getMeshSegmentTags()
        
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
        zone = 808
        geo_ref_poly = Geo_reference(zone, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)
        boundary_tags = {'50': [0],
                         '40': [1],
                         '30': [2],
                         'no where seg': [3],
                         '20': [4]}
        m = create_mesh_from_regions(polygon,
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        m.export_mesh_file('a_test_mesh_iknterface.tsh')

        fileName = tempfile.mktemp('.txt')
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        m.export_mesh_file('a_test_mesh_iknterface.tsh')

        fileName = tempfile.mktemp('.txt')
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
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        # These are the absolute values
        polygon = [[0,0], [100,0], [100,100], [0,100]]

        boundary_tags = {'wall': [0,1,3], 'wave': [2]}
        inner1_polygon = [[10,10], [20,10], [20,20], [10,20]]
        inner2_polygon = [[30,30], [40,30], [40,40], [30,40]]

        max_area = 1
        interior_regions = [(inner1_polygon, 5), (inner2_polygon, 10)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     max_area,
                                     interior_regions=interior_regions)

        fileName = tempfile.mktemp('.txt')
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
        '''Test that create_mesh_from_regions fails
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
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)

        boundary_tags = {'walls': [0,1,3,4], 'bom': [2]}
        try:
            m = create_mesh_from_regions(polygon,
                                         boundary_tags,
                                         10000000,
                                         interior_regions=interior_regions)
        except:
            pass
        else:
            msg = 'Segment out of bounds not caught '
            raise Exception, msg

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)

