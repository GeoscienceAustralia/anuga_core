#!/usr/bin/env python
#
import tempfile
import unittest
import os
from anuga.pmesh.mesh import importMeshFromFile
from anuga.pmesh.mesh_interface import create_mesh_from_regions

# This not work on DSG's desktop
#from anuga.pmesh.mesh_interface import _create_mesh_from_regions
#But this does.  maybe due to me using python 2.3?
from mesh_interface import _create_mesh_from_regions
from load_mesh.loadASCII import *
from anuga.utilities.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference,DEFAULT_ZONE

class TestCase(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_create_mesh_from_regions(self):
        x=-500
        y=-1000
        mesh_geo = geo_reference=Geo_reference(56,x,y)
        
        # These are the absolute values
        polygon_absolute = [[0,0],[100,0],[100,100],[0,100]]
        
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        inner1_polygon_absolute = [[10,10],[20,10],[20,20],[10,20]]
        inner1_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30],[40,30],[40,40],[30,40]]
        inner2_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner2_polygon_absolute)
        
        interior_regions = [(inner1_polygon, 5),(inner2_polygon, 0.2)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly,
                                     mesh_geo_reference=mesh_geo)

        # Test the mesh instance
        self.failUnless(len(m.regions)==3,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==12,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==12,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='walls',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='walls',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[0]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          polygon_absolute, closed = False),
                        'FAILED!')
        
        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[1]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed = False),
                        'FAILED!')
        
        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[2]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          inner2_polygon_absolute,
                                          closed = False),
                        'FAILED!')


    def test_create_mesh_from_regions_with_caching(self):
        x=-500
        y=-1000
        mesh_geo = geo_reference=Geo_reference(56,x,y)
        
        # These are the absolute values
        polygon_absolute = [[0,0],[100,0],[100,100],[0,100]]
        
        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        inner1_polygon_absolute = [[10,10],[20,10],[20,20],[10,20]]
        inner1_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30],[40,30],[40,40],[30,40]]
        inner2_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner2_polygon_absolute)
        
        interior_regions = [(inner1_polygon, 5),(inner2_polygon, 0.2)]

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
        self.failUnless(len(m.regions)==3,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==12,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==12,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='walls',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='walls',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')

        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[0]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          polygon_absolute, closed = False),
                        'FAILED!')
        
        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[1]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          inner1_polygon_absolute,
                                          closed = False),
                        'FAILED!')
        
        # Assuming the order of the region points is known.
        # (This isn't true, if you consider create_mesh_from_regions
        # a black box)
        poly_point = m.getRegions()[2]
        
        #print "poly_point", poly_point
        #print "polygon_absolute",polygon_absolute
         
        # poly_point values are relative to the mesh geo-ref
        # make them absolute
        self.failUnless(is_inside_polygon([poly_point.x+x,poly_point.y+y],
                                          inner2_polygon_absolute,
                                          closed = False),
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
        polygon_absolute = [[min_x,min_y],[1000,100],[1000,1000],[100,1000]]
        
        x_p = -10
        y_p = -40
        zone = 808
        geo_ref_poly = Geo_reference(zone, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        inner1_polygon_absolute = [[10,10],[20,10],[20,20],[10,20]]
        inner1_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30],[40,30],[40,40],[30,40]]
        inner2_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner2_polygon_absolute)
        
        interior_regions = [(inner1_polygon, 5),(inner2_polygon, 0.2)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,
                                     poly_geo_reference=geo_ref_poly)
        

        # Test the mesh instance
        self.failUnless(len(m.regions)==3,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==12,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==12,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='walls',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='walls',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')
        
        self.failUnless(m.geo_reference.get_zone()==zone,
                        'FAILED!')
        self.failUnless(m.geo_reference.get_xllcorner()==min_x,
                        'FAILED!')
        self.failUnless(m.geo_reference.get_yllcorner()==min_y,
                        'FAILED!')

        
    def test_create_mesh_from_regions3(self):

        # These are the absolute values
        min_x = -10 
        min_y = -88
        polygon = [[min_x,min_y],[1000,100],[1000,1000],[100,1000]]
        

        x_p = -10
        y_p = -40
        geo_ref_poly = Geo_reference(56, x_p, y_p)
        
        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        inner1_polygon_absolute = [[10,10],[20,10],[20,20],[10,20]]
        inner1_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner1_polygon_absolute)

        inner2_polygon_absolute = [[30,30],[40,30],[40,40],[30,40]]
        inner2_polygon = geo_ref_poly. \
                         change_points_geo_ref(inner2_polygon_absolute)
        
        interior_regions = [(inner1_polygon, 5),(inner2_polygon, 0.2)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)
        

        # Test the mesh instance
        self.failUnless(len(m.regions)==3,
                        'FAILED!')
        segs = m.getUserSegments()
        self.failUnless(len(segs)==12,
                        'FAILED!')
        self.failUnless(len(m.userVertices)==12,
                        'FAILED!') 
        self.failUnless(segs[0].tag=='walls',
                        'FAILED!')  
        self.failUnless(segs[1].tag=='walls',
                        'FAILED!') 
         
        self.failUnless(segs[2].tag=='bom',
                        'FAILED!') 
        self.failUnless(segs[3].tag=='',
                        'FAILED!')
        
        self.failUnless(m.geo_reference.get_zone()==DEFAULT_ZONE,
                        'FAILED!')
        self.failUnless(m.geo_reference.get_xllcorner()==min_x,
                        'FAILED!')
        self.failUnless(m.geo_reference.get_yllcorner()==min_y,
                        'FAILED!')

    def test_create_mesh_from_regions4(self):

        file_name = tempfile.mktemp(".tsh")
        
        # These are the absolute values
        density_outer = 1000
        min_outer = 0 
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer],[max_outer,min_outer],
                   [max_outer,max_outer],[min_outer,max_outer]]
        
        density_inner1 = 10000000
        inner_buffer = 100
        min_inner1 = min_outer + inner_buffer
        max_inner1 = max_outer - inner_buffer
        inner1_polygon = [[min_inner1,min_inner1],[max_inner1,min_inner1],
                   [max_inner1,max_inner1],[min_inner1,max_inner1]]
      
        
        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        interior_regions = [(inner1_polygon, density_inner1)]
        create_mesh_from_regions(polygon_outer
                                     , boundary_tags
                                     , density_outer
                                     , interior_regions=interior_regions
                                     ,filename=file_name
                                     #,verbose=True
                                     ,verbose=False
                                     )
        
        m = importMeshFromFile(file_name)
        
        #print "file_name",file_name 
        self.failUnless(len(m.meshTriangles) <= 900,
                        'Test mesh interface failed!')
        self.failUnless(len(m.meshTriangles) >= 200,
                        'Test mesh interface failed!')
        
        create_mesh_from_regions(polygon_outer
                                     , boundary_tags
                                     , interior_regions=interior_regions
                                     ,filename=file_name
                                     #,verbose=True
                                     ,verbose=False
                                     )
        
        m = importMeshFromFile(file_name)
        
        #print "len(m.meshTriangles)",len(m.meshTriangles) 
        self.failUnless(len(m.meshTriangles) <= 100,
                        'Test mesh interface failed!')

        os.remove(file_name)
        
    def test_create_mesh_from_regions5(self):

        file_name = tempfile.mktemp(".tsh")
        
        # These are the absolute values
        density_outer = 10000000 
        min_outer = 0 
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer],[max_outer,min_outer],
                   [max_outer,max_outer],[min_outer,max_outer]]
        
        density_inner1 = 1000
        inner_buffer = 100
        min_inner1 = min_outer + inner_buffer
        max_inner1 = max_outer - inner_buffer
        inner1_polygon = [[min_inner1,min_inner1],[max_inner1,min_inner1],
                   [max_inner1,max_inner1],[min_inner1,max_inner1]]
      
        
        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        interior_regions = [(inner1_polygon, density_inner1)]
        create_mesh_from_regions(polygon_outer
                                     , boundary_tags
                                     , density_outer
                                     , interior_regions=interior_regions
                                     ,filename=file_name
                                     #,verbose=True
                                     ,verbose=False
                                     )
        
        m = importMeshFromFile(file_name)
        #print "file_name",file_name
        #print "len(m.meshTriangles",len(m.meshTriangles) 
        self.failUnless(len(m.meshTriangles) <= 2000, 
                        'Test mesh interface failed!')
 
        self.failUnless(len(m.meshTriangles) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)
        
    def test_create_mesh_from_regions6(self):

        file_name = tempfile.mktemp(".tsh")
        
        # These are the absolute values
        density_outer = 1000
        min_outer = 0 
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer],[max_outer,min_outer],
                         [max_outer,max_outer],[min_outer,max_outer]]

        delta = 10
        density_inner1 = 1000
        min_inner1 = min_outer + delta
        max_inner1 = max_outer - delta
        inner1_polygon = [[min_inner1,min_inner1],[max_inner1,min_inner1],
                          [max_inner1,max_inner1],[min_inner1,max_inner1]]
      
        
        density_inner2 = 10000000 
        min_inner2 = min_outer +  2*delta
        max_inner2 = max_outer -  2*delta
        inner2_polygon = [[min_inner2,min_inner2],[max_inner2,min_inner2],
                          [max_inner2,max_inner2],[min_inner2,max_inner2]]
        
        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        interior_regions = [(inner1_polygon, density_inner1),
                            (inner2_polygon, density_inner2)]
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)
        
        m = importMeshFromFile(file_name)
        #print "file_name",file_name
        #print "len(m.meshTriangles",len(m.meshTriangles) 
        self.failUnless(len(m.meshTriangles) <= 2000, 
                        'Test mesh interface failed!')
 
        self.failUnless(len(m.meshTriangles) >= 900,
                        'Test mesh interface failed!')

        os.remove(file_name)
        
    def test_create_mesh_from_regions7(self):

        file_name = tempfile.mktemp(".tsh")
        
        # These are the absolute values
        density_outer = 1001
        min_outer = 0 
        max_outer = 1000
        polygon_outer = [[min_outer,min_outer],[max_outer,min_outer],
                   [max_outer,max_outer],[min_outer,max_outer]]

        delta = 10
        density_inner1 = 100000000
        min_inner1 = min_outer + delta
        max_inner1 = max_outer - delta
        inner1_polygon = [[min_inner1,min_inner1],[max_inner1,min_inner1],
                   [max_inner1,max_inner1],[min_inner1,max_inner1]]
      
        
        density_inner2 = 1000 
        min_inner2 = min_outer +  2*delta
        max_inner2 = max_outer -  2*delta
        inner2_polygon = [[min_inner2,min_inner2],[max_inner2,min_inner2],
                   [max_inner2,max_inner2],[min_inner2,max_inner2]]
        
        boundary_tags = {'walls':[0,1],'bom':[2]}

        #Note the list order is important
        # The last region added will be the region triangle uses,
        # if two regions points are in the same bounded area.
        interior_regions = [(inner2_polygon, density_inner2),(inner1_polygon, density_inner1)]
        create_mesh_from_regions(polygon_outer,
                                 boundary_tags,
                                 density_outer,
                                 interior_regions=interior_regions,
                                 filename=file_name,
                                 verbose=False)
        
        m = importMeshFromFile(file_name)
        #print "file_name",file_name
        #print "len(m.meshTriangles",len(m.meshTriangles) 
        self.failUnless(len(m.meshTriangles) <= 3000, 
                        'Test mesh interface failed!')
 
        self.failUnless(len(m.meshTriangles) >= 2000,
                        'Test mesh interface failed!')

        os.remove(file_name)

        
    def test_create_mesh_from_regions_interior_regions(self):
        """Test that create_mesh_from_regions fails when an interior region is
         outside bounding polygon.       """
        

        # These are the absolute values
        min_x = 10 
        min_y = 88
        polygon = [[min_x,min_y],[1000,100],[1000,1000],[100,1000]]

        
        boundary_tags = {'walls':[0,1],'bom':[2]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400],[900,500],[800,600]]
        
        interior_regions = [(inner_polygon, 5)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions)


        # This one sticks outside bounding polygon - should fail
        inner_polygon = [[800,400],[1100,500],[800,600]]
        interior_regions = [(inner_polygon, 5)]


        
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
            raise msg



    def FIXME_test_create_mesh_with_multiply_tagged_segments(self):
        """Test that create_mesh_from_regions fails when
        segments are listed repeatedly in boundary_tags.
        """
        
        
        

        # These are the absolute values
        min_x = 10 
        min_y = 88
        polygon = [[min_x,min_y],[1000,100],[1000,1000],[100,1000]]

        
        boundary_tags = {'walls':[0,1],'bom':[1,2]}

        # This one is inside bounding polygon - should pass
        inner_polygon = [[800,400],[900,500],[800,600]]
        
        interior_regions = [(inner_polygon, 5)]
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     interior_regions=interior_regions,verbose=False)


        # This one sticks outside bounding polygon - should fail
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
            raise msg

        


    def test_create_mesh_from_regions_with_duplicate_verts(self):

        # These are the absolute values
        
        polygon_absolute = [[0.0,0.0],
                            [0,4.0],
                            [4.0,4.0],
                            [4.0,0.0],
                            [4.0,0.0]]
        
        x_p = -10
        y_p = -40
        zone = 808
        geo_ref_poly = Geo_reference(zone, x_p, y_p)
        polygon = geo_ref_poly.change_points_geo_ref(polygon_absolute)

        boundary_tags = {'50':[0],
                         '40':[1],
                         '30':[2],
                         'no where seg':[3],
                         '20':[4]
                         }
       
        m = create_mesh_from_regions(polygon,
                                     boundary_tags,
                                     10000000,
                                     poly_geo_reference=geo_ref_poly,verbose=False)
        

        fileName = "badmesh.tsh"
        #m.export_mesh_file(fileName)
        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(TestCase,'test')
    #suite = unittest.makeSuite(meshTestCase,'test_asciiFile')
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
    
