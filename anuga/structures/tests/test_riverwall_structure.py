#!/usr/bin/env python


import unittest
import os

#from anuga.structures.riverwall import Boyd_box_operator
#from anuga.structures.riverwall import boyd_box_function

#from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
import anuga
import numpy
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.utilities import plot_utils as util
from anuga.config import g

## Generic data for scenario

verbose = False
boundaryPolygon=[ [0., 0.], [0., 100.], [100.0, 100.0], [100.0, 0.0]]
wallLoc=50.
# The boundary polygon + riverwall breaks the mesh into multiple regions
# Must define the resolution in these areas with an xy point + maximum area
# Otherwise triangle.c gets confused
regionPtAreas=[ [99., 99., 20.0*20.0*0.5],
                [1., 1., 20.0*20.0*0.5]]

class Test_riverwall_structure(unittest.TestCase):
    """
	Test the riverwall structure
    """

    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('test_riverwall.sww')
        except:
            pass

        try:
            os.remove('testRiverwall.msh')
        except:
            pass
    
    def create_domain_DE0(self, wallHeight, InitialOceanStage, InitialLandStage, riverWall=None, riverWall_Par=None):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values

        if(riverWall is None):
            riverWall={ 'centralWall':
                            [ [wallLoc, 0.0, wallHeight],
                              [wallLoc, 100.0, wallHeight]] 
                      }
        if(riverWall_Par is None):
            riverWall_Par={'centralWall':{'Qfactor':1.0}}
        # Make the domain
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'testRiverwall.msh',
                                   interior_regions =[ ], #[ [higherResPolygon, 1.*1.*0.5],
                                                          #  [midResPolygon, 3.0*3.0*0.5]],
                                   breaklines=list(riverWall.values()),
                                   use_cache=False,
                                   verbose=verbose,
                                   regionPtArea=regionPtAreas)

        domain=anuga.create_domain_from_file('testRiverwall.msh')

        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm('DE0')
        domain.set_name('test_riverwall')

        domain.set_store_vertices_uniquely()
       
        def topography(x,y):
            return -x/150. 

        def stagefun(x,y):
            stg=InitialOceanStage*(x>=50.) + InitialLandStage*(x<50.)
            return stg 

        # NOTE: Setting quantities at centroids is important for exactness of tests
        domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('friction',0.03)             
        domain.set_quantity('stage', stagefun,location='centroids')            
        
        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=verbose) 

        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

        return domain
 
    def create_domain_DE1(self, wallHeight, InitialOceanStage, InitialLandStage):
        # Riverwall = list of lists, each with a set of x,y,z (and optional QFactor) values
        riverWall={ 'centralWall':
                        [ [wallLoc, 0.0, wallHeight],
                          [wallLoc, 100.0, wallHeight]] 
                  }
        riverWall_Par={'centralWall':{'Qfactor':1.0}}
        # Make the domain
        anuga.create_mesh_from_regions(boundaryPolygon, 
                                 boundary_tags={'left': [0],
                                                'top': [1],
                                                'right': [2],
                                                'bottom': [3]},
                                   maximum_triangle_area = 200.,
                                   minimum_triangle_angle = 28.0,
                                   filename = 'testRiverwall.msh',
                                   interior_regions =[ ], #[ [higherResPolygon, 1.*1.*0.5],
                                                          #  [midResPolygon, 3.0*3.0*0.5]],
                                   breaklines=list(riverWall.values()),
                                   use_cache=False,
                                   verbose=verbose,
                                   regionPtArea=regionPtAreas)

        domain=anuga.create_domain_from_file('testRiverwall.msh')

        # 05/05/2014 -- riverwalls only work with DE0 and DE1
        domain.set_flow_algorithm('DE1')
        domain.set_name('test_riverwall')

        domain.set_store_vertices_uniquely()
       
        def topography(x,y):
            return -x/150. 

        def stagefun(x,y):
            stg=InitialOceanStage*(x>=50.) + InitialLandStage*(x<50.)
            return stg 

        # NOTE: Setting quantities at centroids is important for exactness of tests
        domain.set_quantity('elevation',topography,location='centroids')     
        domain.set_quantity('friction',0.03)             
        domain.set_quantity('stage', stagefun,location='centroids')            
        
        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=verbose) 

        # Boundary conditions
        Br=anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

        return domain
           
    def test_noflux_riverwall_DE1(self):
        """test_noflux_riverwall
            
            Tests that the riverwall blocks water when the stage is < wall height 
            
        """
        wallHeight=-0.2
        InitialOceanStage=-0.3
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 

        # Run the model for a few seconds, and check that no water has flowed past the riverwall 
        for t in domain.evolve(yieldstep=10.0,finaltime=10.0):
            if(verbose): 
                print(domain.timestepping_statistics())

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', 0.)

        assert numpy.allclose(landVol,0., atol=1.0e-12)
        
    def test_simpleflux_riverwall_DE1(self):
        """test_simpleflux_riverwall
            
            Tests that the riverwall flux (when dry on one edge) is
            2/3*H*sqrt(2/3*g*H)
            
        """
        wallHeight=2.
        InitialOceanStage=2.50
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        H=InitialOceanStage-wallHeight # Upstream head
        dt=ft
        theoretical_flux_vol=dt*L*2./3.*H*(2./3.*g*H)**0.5

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()

        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)
    
    def test_simpleflux_riverwall_with_Qfactor_DE1(self):
        """test_simpleflux_riverwall with Qfactor != 1.0
            
            Tests that the riverwall flux (when dry on one edge) is
            2/3*H*sqrt(2/3*g*H)*Qfactor
            
        """
        wallHeight=2.
        InitialOceanStage=2.50
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 
        # Redefine the riverwall, with 2x the discharge Qfactor
        riverWall={ 'centralWall':
                        [ [wallLoc, 0.0, wallHeight],
                          [wallLoc, 100.0, wallHeight]] 
                  }
        Qmult=2.0
        riverWall_Par={'centralWall':{'Qfactor':Qmult}}
        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=verbose)

        #import pdb
        #pdb.set_trace()

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        H=InitialOceanStage-wallHeight # Upstream head
        dt=ft
        theoretical_flux_vol=Qmult*dt*L*2./3.*H*(2./3.*g*H)**0.5

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()

        #import pdb
        #pdb.set_trace()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)

   
    def test_submergeflux_riverwall_DE1(self):
        """test_submergedflux_riverwall
            
            Tests that the riverwall flux is
            Q1*(1-Q2/Q1)**(0.385)

            Where Q1, Q2 =2/3*H*sqrt(2/3*g*H)
             with H computed from the smallest (Q2) / largest (Q1) head side
        """
        wallHeight=20.
        InitialOceanStage=20.50
        InitialLandStage=20.44
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        InitialHeight=(domain.quantities['stage'].centroid_values[landInds]-domain.quantities['elevation'].centroid_values[landInds])
        InitialLandVol=InitialHeight*(InitialHeight>0.)*domain.areas[landInds]            
        InitialLandVol=InitialLandVol.sum()

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        dt=ft
        H1=max(InitialOceanStage-wallHeight,0.) # Upstream head
        H2=max(InitialLandStage-wallHeight,0.) # Downstream head

        Q1=2./3.*H1*(2./3.*g*H1)**0.5
        Q2=2./3.*H2*(2./3.*g*H2)**0.5

        theoretical_flux_vol=dt*L*Q1*(1.-Q2/Q1)**0.385

        # Compute volume of water landward of riverwall
        FinalLandVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        FinalLandVol=FinalLandVol.sum()

        landVol=FinalLandVol-InitialLandVol

        #import pdb
        #pdb.set_trace()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)
 
    def test_noflux_riverwall_DE0(self):
        """test_noflux_riverwall
            
            Tests that the riverwall blocks water when the stage is < wall height 
            
        """
        wallHeight=-0.2
        InitialOceanStage=-0.3
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE0(wallHeight,InitialOceanStage, InitialLandStage) 

        # Run the model for a few seconds, and check that no water has flowed past the riverwall 
        for t in domain.evolve(yieldstep=10.0,finaltime=10.0):
            if(verbose): 
                print(domain.timestepping_statistics())

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', 0.)

        assert numpy.allclose(landVol,0., atol=1.0e-12)
        
    def test_simpleflux_riverwall_DE0(self):
        """test_simpleflux_riverwall
            
            Tests that the riverwall flux (when dry on one edge) is
            2/3*H*sqrt(2/3*g*H)
            
        """
        wallHeight=2.
        InitialOceanStage=2.50
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE0(wallHeight,InitialOceanStage, InitialLandStage) 

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        H=InitialOceanStage-wallHeight # Upstream head
        dt=ft
        theoretical_flux_vol=dt*L*2./3.*H*(2./3.*g*H)**0.5

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()

        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)
    
    def test_simpleflux_riverwall_with_Qfactor_DE0(self):
        """test_simpleflux_riverwall with Qfactor != 1.0
            
            Tests that the riverwall flux (when dry on one edge) is
            2/3*H*sqrt(2/3*g*H)*Qfactor
            
        """
        wallHeight=2.
        InitialOceanStage=2.50
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE0(wallHeight,InitialOceanStage, InitialLandStage) 
        # Redefine the riverwall, with 2x the discharge Qfactor
        riverWall={ 'centralWall':
                        [ [wallLoc, 0.0, wallHeight],
                          [wallLoc, 100.0, wallHeight]] 
                  }
        Qmult=2.0
        riverWall_Par={'centralWall':{'Qfactor':Qmult}}
        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=verbose)

        #import pdb
        #pdb.set_trace()

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        H=InitialOceanStage-wallHeight # Upstream head
        dt=ft
        theoretical_flux_vol=Qmult*dt*L*2./3.*H*(2./3.*g*H)**0.5

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()

        #import pdb
        #pdb.set_trace()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)

   
    def test_submergeflux_riverwall_DE0(self):
        """test_submergedflux_riverwall
            
            Tests that the riverwall flux is
            Q1*(1-Q2/Q1)**(0.385)

            Where Q1, Q2 =2/3*H*sqrt(2/3*g*H)
             with H computed from the smallest (Q2) / largest (Q1) head side
        """
        wallHeight=20.
        InitialOceanStage=20.50
        InitialLandStage=20.44
       
        domain=self.create_domain_DE0(wallHeight,InitialOceanStage, InitialLandStage) 

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        InitialHeight=(domain.quantities['stage'].centroid_values[landInds]-domain.quantities['elevation'].centroid_values[landInds])
        InitialLandVol=InitialHeight*(InitialHeight>0.)*domain.areas[landInds]            
        InitialLandVol=InitialLandVol.sum()

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        dt=ft
        H1=max(InitialOceanStage-wallHeight,0.) # Upstream head
        H2=max(InitialLandStage-wallHeight,0.) # Downstream head

        Q1=2./3.*H1*(2./3.*g*H1)**0.5
        Q2=2./3.*H2*(2./3.*g*H2)**0.5

        theoretical_flux_vol=dt*L*Q1*(1.-Q2/Q1)**0.385

        # Compute volume of water landward of riverwall
        FinalLandVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        FinalLandVol=FinalLandVol.sum()

        landVol=FinalLandVol-InitialLandVol

        #import pdb
        #pdb.set_trace()
        
        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)
    
    def test_riverwall_includes_specified_points_in_domain(self): 
        """
            Check that all domain points that should be on the riverwall
            actually are, and that there are no 'non-riverwall' points on the riverwall
        """
        wallHeight=-0.2
        InitialOceanStage=-0.3
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 

        edgeInds_on_wall=domain.riverwallData.riverwall_edges.tolist()

        riverWall_x_coord=domain.edge_coordinates[edgeInds_on_wall,0]-wallLoc

        assert(numpy.allclose(riverWall_x_coord,0.))
   
        # Now check that all the other domain edge coordinates are not on the wall 
        # Note the threshold requires a sufficiently coarse mesh
        notriverWall_x_coord=numpy.delete(domain.edge_coordinates[:,0], edgeInds_on_wall) 
        assert(min(abs(notriverWall_x_coord-wallLoc))>1.0e-01)

    def test_is_vertex_on_boundary(self):
        """
            Check that is_vertex_on_boundary is working as expected
        """
        wallHeight=-0.2
        InitialOceanStage=-0.3
        InitialLandStage=-999999.
       
        domain=self.create_domain_DE1(wallHeight,InitialOceanStage, InitialLandStage) 

        allVertices=numpy.array(list(range(len(domain.vertex_coordinates))))
        boundaryFlag=domain.riverwallData.is_vertex_on_boundary(allVertices)
        boundaryVerts=boundaryFlag.nonzero()[0].tolist()

        # Check that all boundary vertices are on the boundary 
        check2=(domain.vertex_coordinates[boundaryVerts,0]==0.)+\
               (domain.vertex_coordinates[boundaryVerts,0]==100.)+\
               (domain.vertex_coordinates[boundaryVerts,1]==100.)+\
               (domain.vertex_coordinates[boundaryVerts,1]==0.)

        assert(all(check2>0.))

        # Check that all non-boundary vertices are not
        nonboundaryVerts=(boundaryFlag==0).nonzero()[0].tolist()
        check2=(domain.vertex_coordinates[nonboundaryVerts,0]==0.)+\
               (domain.vertex_coordinates[nonboundaryVerts,0]==100.)+\
               (domain.vertex_coordinates[nonboundaryVerts,1]==100.)+\
               (domain.vertex_coordinates[nonboundaryVerts,1]==0.)

        assert(all(check2==0))

    def test_multiple_riverwalls(self):
        """
            Testcase with multiple riverwalls -- check all is working as required

         Idea -- add other riverwalls with different Qfactor / height.
                 Set them up to have no hydraulic effect, but
                 so that we are likely to catch bugs if the code is not right
        """
        wallHeight=2.
        InitialOceanStage=2.50
        InitialLandStage=-999999.
        
       
        riverWall={ 'awall1':
                        [ [wallLoc+20., 0.0, -9999.],
                          [wallLoc+20., 100.0, -9999.]],  
                    'centralWall':
                        [ [wallLoc, 0.0, wallHeight],
                          [wallLoc, 100.0, wallHeight]] ,
                    'awall2':
                        [ [wallLoc-20., 0.0, 30.],
                          [wallLoc-20., 100.0, 30.]],  
                  }

        newQfac=2.0
        riverWall_Par={'centralWall':{'Qfactor':newQfac}, 'awall1':{'Qfactor':100.}, 'awall2':{'Qfactor':0.}}

        domain=self.create_domain_DE0(wallHeight,InitialOceanStage, InitialLandStage, riverWall=riverWall,riverWall_Par=riverWall_Par) 
        

        domain.riverwallData.create_riverwalls(riverWall,riverWall_Par,verbose=verbose) 

        # Run the model for a few fractions of a second
        # Any longer, and the evolution of stage starts causing
        # significant changes to the flow
        yst=1.0e-04
        ft=1.0e-03
        for t in domain.evolve(yieldstep=yst,finaltime=ft):
            if(verbose):
                print(domain.timestepping_statistics())

        # Compare with theoretical result
        L= 100. # Length of riverwall
        H=InitialOceanStage-wallHeight # Upstream head
        dt=ft
        theoretical_flux_vol=newQfac*dt*L*2./3.*H*(2./3.*g*H)**0.5

        # Indices landward of the riverwall
        landInds=(domain.centroid_coordinates[:,0]<50.).nonzero()[0]
        # Compute volume of water landward of riverwall
        landVol=domain.quantities['height'].centroid_values[landInds]*domain.areas[landInds]            
        landVol=landVol.sum()

        if(verbose):
            print('Land Vol: ', landVol, 'theoretical vol: ', theoretical_flux_vol)

        assert numpy.allclose(landVol,theoretical_flux_vol, rtol=1.0e-03)

# =========================================================================
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_riverwall_structure, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
