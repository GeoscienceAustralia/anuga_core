#!/usr/bin/env python

import os
import unittest
import anuga
import anuga.utilities.plot_utils as util
import numpy as np

flow_algorithms = ['DE0', 'DE1', '1_5', '2_0', 'tsunami']
verbose = False

class Test_plot_utils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def create_domain(self, InitialOceanStage, InitialLandStage, flowAlg, verbose):
        """
         Make the domain and set the flow algorithm for a test. Produces an sww
         that we can use for testing
        """
        boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                   boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 20.,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'test_plot_utils.msh',
                                   interior_regions =[ ],
                                   verbose=False)

        domain=anuga.create_domain_from_file('test_plot_utils.msh')
        
        os.remove('test_plot_utils.msh')

        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm(flowAlg)
        domain.set_name('test_plot_utils')

        domain.set_store_vertices_uniquely()
       
        def topography(x,y):
            return -x/150. 

        def stagefun(x,y):
            stg=InitialOceanStage*(x>=50.) + InitialLandStage*(x<50.)
            return stg 

        domain.set_flow_algorithm(flowAlg)
        #domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('elevation',topography)     
        domain.set_quantity('friction',0.03)             
        #domain.set_quantity('stage', stagefun,location='centroids')            
        domain.set_quantity('stage', stagefun)            
       
        if(verbose):
            if(domain.store_centroids):
                print '   Centroids stored'
            else:
                print '    Centroids estimated from vertices'

        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

        for t in domain.evolve(yieldstep=0.2,finaltime=1.0):
            pass

        return 

    def everything_equal(self, p1, slice1, p2, slice2):
        """
            Convenience function to check that all variables in two 'get_output' objects
            are the same.
            
        """
        # Array variables
        assert(all(p1.stage[slice1,:]==p2.stage[slice2,:]))
        assert(all(p1.vel[slice1,:]==p2.vel[slice2,:]))
        assert(all(p1.height[slice1,:]==p2.height[slice2,:]))
        assert(all(p1.xmom[slice1,:]==p2.xmom[slice2,:]))
        assert(all(p1.ymom[slice1,:]==p2.ymom[slice2,:]))
        assert(all(p1.xvel[slice1,:]==p2.xvel[slice2,:]))
        assert(all(p1.yvel[slice1,:]==p2.yvel[slice2,:]))

        # Constant or 1D variables
        assert(all(p1.elev==p2.elev))
        assert(all(p1.x==p2.x))
        assert(all(p1.y==p2.y))
        assert(all(p1.friction==p2.friction))
        assert(p1.xllcorner==p2.xllcorner)
        assert(p1.yllcorner==p2.yllcorner)
        assert((p1.time[slice1]==p2.time[slice2]))
        return 

    def basic_checks(self):
        """
            Check that dimensions are as required, and that
            if we extract centroids by passing the file name
            we get the same result as if we pass the get_output object
        """

        p=util.get_output('test_plot_utils.sww')
       
        # Check that dimesions are ok
        l_time=len(p.time) 
        l_space=len(p.x)

        assert(len(p.y)==l_space)
        assert(p.stage.shape==(l_time,l_space)) 
        assert(p.vel.shape==(l_time,l_space)) 
        assert(p.xmom.shape==(l_time,l_space)) 
        assert(p.ymom.shape==(l_time,l_space)) 
        assert(p.xvel.shape==(l_time,l_space)) 
        assert(p.yvel.shape==(l_time,l_space)) 
        assert(p.elev.shape==(l_space,))
        assert(p.friction.shape==(l_space,))
       
        p2=util.get_centroids(p,velocity_extrapolation=True) 
        l_space=len(p.vols) #  len(vols in vertex quantities) = len(centroids)
        assert(p2.stage.shape==(l_time,l_space)) 
        assert(p2.vel.shape==(l_time,l_space)) 
        assert(p2.xmom.shape==(l_time,l_space)) 
        assert(p2.ymom.shape==(l_time,l_space)) 
        assert(p2.xvel.shape==(l_time,l_space)) 
        assert(p2.yvel.shape==(l_time,l_space)) 
        assert(p2.elev.shape==(l_space,))
        assert(p2.friction.shape==(l_space,))

        # Read centroids as a file, and check that it is the same as above for a couple of time-slices
        p3=util.get_centroids('test_plot_utils.sww',velocity_extrapolation=True)
        self.everything_equal(p3, 2, p2, 2)
        self.everything_equal(p3, 4, p2, 4)

    def velExtrap_timeSlices_check(self, ve):
        """

            Check that time-slices are behaving as expected. Assumes sww has been made

            Is called by a test function for various flow algorithms
        """
        p=util.get_output('test_plot_utils.sww')
        p2=util.get_centroids(p,velocity_extrapolation=ve) 
       
        # Check that dimesions are ok
        l_time=len(p.time) 
        l_space=len(p.x)

        assert(p.timeSlices==range(l_time))
        assert(p2.timeSlices==range(l_time))

        # Try getting some time-slices, and checking all is as intended
        p_12=util.get_output('test_plot_utils.sww', timeSlices=[0, 3])
        pc_12=util.get_centroids(p_12,velocity_extrapolation=ve)

        assert(p_12.timeSlices==[0,3])
        assert(pc_12.timeSlices==[0,3])

        self.everything_equal(p_12, 0, p, 0)
        self.everything_equal(p_12, 1, p, 3)
        self.everything_equal(pc_12, 0, p2, 0)
        self.everything_equal(pc_12, 1, p2, 3)
        

        # Try getting some time-slices, and checking all is as intended
        p_12a=util.get_output('test_plot_utils.sww', timeSlices=3)
        pc_12a=util.get_centroids(p_12a,velocity_extrapolation=ve)

        #print p_12a.timeSlices
        #print pc_12a.timeSlices
        assert(p_12a.timeSlices==[3])
        assert(pc_12a.timeSlices==[3])

        self.everything_equal(p_12a, 0, p, 3)
        self.everything_equal(pc_12a, 0, p2, 3)
       

        # Try getting some time-slices, and checking all is as intended
        p_12b=util.get_output('test_plot_utils.sww')
        pc_12b=util.get_centroids(p_12b,velocity_extrapolation=ve, timeSlices=3)

        #print p_12b.timeSlices
        #print pc_12b.timeSlices
        assert(p_12b.timeSlices==[0, 1, 2, 3, 4, 5])
        assert(pc_12b.timeSlices==[3])

        self.everything_equal(p_12b, 0, p, 0)
        self.everything_equal(p_12b, 5, p, 5)
        self.everything_equal(pc_12b, 0, p2, 3)
       

        # Check we can get the 'last' time, and it is correct 
        p_l=util.get_output('test_plot_utils.sww', timeSlices='last')
        pc_l=util.get_centroids(p_l,velocity_extrapolation=ve)
        l_time=len(p.time) 
        self.everything_equal(p_l, 0, p, l_time-1)
        self.everything_equal(pc_l, 0, p2, l_time-1)
        
        assert(p_l.timeSlices==[l_time-1])
        assert(pc_l.timeSlices==[l_time-1])
        
        # Check that we can get the 'max' time
        p_m=util.get_output('test_plot_utils.sww', timeSlices='max')
        pc_m=util.get_centroids(p_m,velocity_extrapolation=ve)
        
        assert(p_m.time==p.time.max())
        assert(pc_m.time==p2.time.max())
        assert(p_m.timeSlices=='max')
        assert(pc_m.timeSlices=='max')
        assert(all(p_m.stage[0,:]==p.stage.max(axis=0)))
        assert(all(pc_m.stage[0,:]==p2.stage.max(axis=0)))
        assert(all(p_m.vel[0,:]==p.vel.max(axis=0)))
        assert(all(pc_m.vel[0,:]==p2.vel.max(axis=0)))
        assert(all(p_m.height[0,:]==p.height.max(axis=0)))
        assert(all(pc_m.height[0,:]==p2.height.max(axis=0)))
        # Somewhat lazy test of variables where the sign is important
        assert(all(np.abs(p_m.xmom)[0,:]==np.abs(p.xmom).max(axis=0)))
        assert(all(np.abs(pc_m.xmom)[0,:]==np.abs(p2.xmom).max(axis=0)))
        assert(all(np.abs(p_m.ymom)[0,:]==np.abs(p.ymom).max(axis=0)))
        assert(all(np.abs(pc_m.ymom)[0,:]==np.abs(p2.ymom).max(axis=0)))
        assert(all(np.abs(p_m.xvel)[0,:]==np.abs(p.xvel).max(axis=0)))
        assert(all(np.abs(pc_m.xvel)[0,:]==np.abs(p2.xvel).max(axis=0)))
        assert(all(np.abs(p_m.yvel)[0,:]==np.abs(p.yvel).max(axis=0)))
        assert(all(np.abs(pc_m.yvel)[0,:]==np.abs(p2.yvel).max(axis=0)))

        return

    def quantity_consistency_check(self, p, slice1):
        """
            Check logical things such as that e.g. xmom = xvel*depth, etc
        """
        assert(np.allclose(p.xmom[slice1,:],p.xvel[slice1,:]*p.height[slice1,:]))
        assert(np.allclose(p.ymom[slice1,:],p.yvel[slice1,:]*p.height[slice1,:]))
        assert(np.allclose(p.stage[slice1,:], p.height[slice1,:]+p.elev, atol=1.0e-07))
        assert(np.allclose(p.vel[slice1,:], (p.xvel[slice1,:]**2. + p.yvel[slice1,:]**2)**0.5 ))
        return
        

    def test_basic(self):
        """
            Check dimensions of important quantities
        """
        # Run flowAlg_test_basic for a range of flow algorithms
        for flowAlg in flow_algorithms:
            if(verbose):
                print flowAlg
            self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg=flowAlg, verbose=verbose)
            self.basic_checks()

        os.remove('test_plot_utils.sww')
        

    def test_timeslices(self):
        """
            Check that outputs from timeslice-subsets agree with bulk outputs
        """
        for flowAlg in flow_algorithms:
            if(verbose):
                print flowAlg
            self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg=flowAlg, verbose=verbose)
            # Test time-slices with velocity_extrapolation=True
            self.velExtrap_timeSlices_check(ve=True)
            self.velExtrap_timeSlices_check(ve=False)
            
        os.remove('test_plot_utils.sww')

    def test_quantity_consistency(self):
        """
            Check that the expected relations exist between stage, elevation,
                height, xmom,xvel,ymom,yvel
        """
        for flowAlg in flow_algorithms:
            if(verbose):
                print flowAlg
            self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg=flowAlg, verbose=verbose)
            pc=util.get_centroids('test_plot_utils.sww')
            l=len(pc.time)-1
            self.quantity_consistency_check(pc,l)
            
        os.remove('test_plot_utils.sww')
      

    def test_near_points(self):
        """
            Check that the near-points function is working
        """ 
        self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg='DE0', verbose=verbose)
        p=util.get_output('test_plot_utils.sww')
        pc=util.get_centroids(p,velocity_extrapolation=True)
        # First check -- get points along y==50
        nt=util.near_transect(pc, [20., 50.], [80., 50.], tol=10.)
        assert(all(abs(pc.y[nt[0]]-50.)<10.))
        assert(all(nt[1]>=0.))
        assert(all(nt[1]<=60.))
        assert(np.allclose(pc.x[nt[0]]-20., nt[1]))
        
        # Next check -- get points along x==50
        nt=util.near_transect(pc, [50., 20.], [50., 80.], tol=10.)
        assert(all(abs(pc.x[nt[0]]-50.)<10.))
        assert(all(nt[1]>=0.))
        assert(all(nt[1]<=60.))
        assert(np.allclose(pc.y[nt[0]]-20., nt[1]))
        
        # Next check -- get points along x==y
        nt=util.near_transect(pc, [20., 20.], [80., 80.], tol=10.)
        assert(all(nt[1]>=0.))
        # Length of line is 60*sqrt(2)
        assert(all(nt[1]<=60.*2**0.5))
        # Coords must be within 10*sqrt(2) of each other
        assert(all(abs( pc.x[nt[0]]-pc.y[nt[0]]) < 10.*2**0.5))
        # The dot product of the points along the line is equal to nt[1]
        dt_Prd=( (pc.x[nt[0]]-20.)/2.**0.5 + (pc.y[nt[0]]-20.)/2.**0.5)
        assert(np.allclose(dt_Prd , nt[1]))
        
        # Next check -- get points along x==2*y + 5
        nt=util.near_transect(pc, [25., 10.], [85., 40.], tol=10.)
        assert(all(nt[1]>=0.))
        # Length of line is sqrt(60^2+30^3)
        assert(all(nt[1]<=(60.**2+30.**2)**0.5))
        # The dot product of the points along the line is equal to nt[1]
        # Unit vector along line is (1,0.5)/ll
        ll=(1.**2+0.5**2)**0.5
        dt_Prd=( (pc.x[nt[0]]-25.)/ll + (pc.y[nt[0]]-10.)*0.5/ll)
        assert(np.allclose(dt_Prd , nt[1]))
        
        os.remove('test_plot_utils.sww')

    def test_triangle_areas(self):
        """
            Check that triangle areas is working as it should
        """
        self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg='DE0', verbose=verbose)
        p=util.get_output('test_plot_utils.sww')
        pc=util.get_centroids(p,velocity_extrapolation=True)

        # Check that subsetting works        
        ta=util.triangle_areas(p)
        ta_1=util.triangle_areas(p,range(10))
        assert(all(ta[range(10)]==ta_1))

        # Independently compute an example and check it
        x0=p.x[p.vols[0][0]]
        y0=p.y[p.vols[0][0]]
        x1=p.x[p.vols[0][1]]
        y1=p.y[p.vols[0][1]]
        x2=p.x[p.vols[0][2]]
        y2=p.y[p.vols[0][2]]
        
        # Use 0.5 a * b
        len_a = ((x1-x0)**2 + (y1-y0)**2)**0.5
        vec_01 = np.array([ x1-x0, y1-y0])
        vec_01 = vec_01/((vec_01**2).sum())**0.5
        vec_01_perp=np.array([vec_01[1], -vec_01[0]])
        len_b = (x2-x0)*vec_01_perp[0] + (y2-y0)*vec_01_perp[1]
        assert(np.allclose(abs(0.5*len_a*len_b),ta[0]))
        
        os.remove('test_plot_utils.sww')

    def test_water_volume(self):
        """ Check that water volume is right
            We assume triangle areas are computed ok, but note they are tested above
        """
        self.create_domain(InitialOceanStage=1., InitialLandStage=0., flowAlg='DE0', verbose=verbose)
        p=util.get_output('test_plot_utils.sww')
        pc=util.get_centroids(p,velocity_extrapolation=True)

        # Check that subsetting works        
        ta=util.triangle_areas(p)

        # Independently computed water volume
        wVol_2=(ta*pc.height[2,:]).sum()

        wv=util.water_volume(p,pc)
        assert(np.allclose(wVol_2, wv[2]))
        
        os.remove('test_plot_utils.sww')

    def test_Make_Geotif(self):
        # VERY BASIC TEST
        #
        # Simply create some data and grid it
        #
        # If nothing fails, that's good
        #
        # Pick a domain that makes sense in EPSG:32756
        x=np.linspace(307000., 308000., 100)     
        y=np.linspace(6193000., 6194000., 100)
        myCellSize=5.0
        xG,yG=np.meshgrid(x, y)
        xG=xG.flatten()
        yG=yG.flatten()
        # Surface is z=x+y
        fakeZ=xG-min(xG)+yG -min(yG)
        dataToGrid=np.vstack([xG,yG,fakeZ]).transpose()
        util.Make_Geotif(dataToGrid, output_quantities=['TestData'],
                         EPSG_CODE=32756, output_dir='.', CellSize=myCellSize)

       
        # Use gdal to check that at least the data extent is ok
        import gdal
        raster=gdal.Open('PointData_TestData.tif')
        rasterGeoTrans=raster.GetGeoTransform()
        assert(np.allclose(x.min()-myCellSize/2.0, rasterGeoTrans[0]))
        assert(np.allclose(y.max()+myCellSize/2.0, rasterGeoTrans[3]))
        #release data file
        raster = None
        # Delete tif made with Make_Geotif
        os.remove('PointData_TestData.tif')

    def test_Make_Geotif_with_knn(self):
        # VERY BASIC TEST using knn+inverse distance interpolation to make the grid
        #
        # Simply create some data and grid it
        #
        # If nothing fails, that's good
        #
        # Pick a domain that makes sense in EPSG:32756
        x=np.linspace(307000., 308000., 100)     
        y=np.linspace(6193000., 6194000., 100)
        myCellSize=5.0
        xG,yG=np.meshgrid(x, y)
        xG=xG.flatten()
        yG=yG.flatten()
        # Surface is z=x+y
        fakeZ=xG-min(xG)+yG -min(yG)
        dataToGrid=np.vstack([xG,yG,fakeZ]).transpose()
        util.Make_Geotif(dataToGrid, output_quantities=['TestData'],
                         EPSG_CODE=32756, output_dir='.', CellSize=myCellSize,
                         k_nearest_neighbours=4)

       
        # Use gdal to check that at least the data extent is ok
        import gdal
        raster=gdal.Open('PointData_TestData.tif')
        rasterGeoTrans=raster.GetGeoTransform()
        assert(np.allclose(x.min()-myCellSize/2.0, rasterGeoTrans[0]))
        assert(np.allclose(y.max()+myCellSize/2.0, rasterGeoTrans[3]))
        #release data file
        raster = None
        # Delete tif made with Make_Geotif
        os.remove('PointData_TestData.tif')
        

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_plot_utils, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

