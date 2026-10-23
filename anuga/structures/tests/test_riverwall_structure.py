#!/usr/bin/env python


import unittest
import os
import pytest

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
        except OSError:
            pass

        try:
            os.remove('testRiverwall.msh')
        except OSError:
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
        anuga.create_pmesh_from_regions(boundaryPolygon,
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
        anuga.create_pmesh_from_regions(boundaryPolygon,
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

    # -----------------------------------------------------------------------
    # Throughflow (Cd_through) tests
    # -----------------------------------------------------------------------

    def test_throughflow_parameter_stored(self):
        """Cd_through is stored in hydraulic_properties table at column index 5."""
        wallHeight = 2.0
        Cd = 0.5
        riverWall = {'centralWall': [[wallLoc, 0.0, wallHeight],
                                     [wallLoc, 100.0, wallHeight]]}
        riverWall_Par = {'centralWall': {'Cd_through': Cd}}
        domain = self.create_domain_DE0(wallHeight, 1.0, 0.5,
                                        riverWall=riverWall,
                                        riverWall_Par=riverWall_Par)
        hp = domain.riverwallData.hydraulic_properties
        names = domain.riverwallData.hydraulic_variable_names
        assert 'Cd_through' in names, "Cd_through missing from hydraulic_variable_names"
        col = list(names).index('Cd_through')
        assert col == 5, f"Cd_through should be column 5, got {col}"
        assert numpy.allclose(hp[0, col], Cd), \
            f"Expected Cd_through={Cd}, got {hp[0, col]}"

    def test_throughflow_default_zero(self):
        """Cd_through defaults to 0 — behaviour matches impermeable wall."""
        wallHeight = 2.0
        InitialOceanStage = 1.0
        InitialLandStage = -999999.

        # With default Cd_through=0 the wall is impermeable below the crest.
        # Stage below the wall on both sides: no overtopping, no throughflow → land stays dry.
        domain = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage)

        for t in domain.evolve(yieldstep=10.0, finaltime=10.0):
            pass

        landInds = (domain.centroid_coordinates[:, 0] < 50.).nonzero()[0]
        landVol = (domain.quantities['height'].centroid_values[landInds]
                   * domain.areas[landInds]).sum()
        assert landVol < 1.0e-6, \
            f"With Cd_through=0 land should stay dry, got volume={landVol}"

    def test_throughflow_dry_downstream(self):
        """With Cd_through>0 and dry downstream side, water flows through the wall."""
        wallHeight = 2.0
        InitialOceanStage = 1.0    # below crest, so no overtopping
        InitialLandStage = -999999.

        riverWall = {'centralWall': [[wallLoc, 0.0, wallHeight],
                                     [wallLoc, 100.0, wallHeight]]}
        riverWall_Par = {'centralWall': {'Cd_through': 0.5}}

        domain = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage,
                                        riverWall=riverWall, riverWall_Par=riverWall_Par)

        for t in domain.evolve(yieldstep=10.0, finaltime=10.0):
            pass

        landInds = (domain.centroid_coordinates[:, 0] < 50.).nonzero()[0]
        landVol = (domain.quantities['height'].centroid_values[landInds]
                   * domain.areas[landInds]).sum()
        assert landVol > 1.0e-4, \
            f"With Cd_through=0.5 and dry downstream, expected water flow, got volume={landVol}"

    def test_throughflow_more_cd_more_flow(self):
        """Higher Cd_through produces more throughflow volume."""
        wallHeight = 2.0
        InitialOceanStage = 1.0
        InitialLandStage = -999999.

        def run(Cd):
            riverWall = {'centralWall': [[wallLoc, 0.0, wallHeight],
                                         [wallLoc, 100.0, wallHeight]]}
            riverWall_Par = {'centralWall': {'Cd_through': Cd}}
            domain = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage,
                                            riverWall=riverWall, riverWall_Par=riverWall_Par)
            for t in domain.evolve(yieldstep=0.001, finaltime=0.001):
                pass
            landInds = (domain.centroid_coordinates[:, 0] < 50.).nonzero()[0]
            return (domain.quantities['height'].centroid_values[landInds]
                    * domain.areas[landInds]).sum()

        vol_lo = run(0.1)
        vol_hi = run(0.5)
        assert vol_hi > vol_lo, \
            f"Higher Cd_through should produce more throughflow: Cd=0.1 → {vol_lo}, Cd=0.5 → {vol_hi}"

    def test_throughflow_additive_to_overtopping(self):
        """Throughflow adds to overtopping: total flux > overtopping-only flux."""
        wallHeight = 0.5
        InitialOceanStage = 1.0   # above crest: overtopping occurs
        InitialLandStage = -999999.

        yst = 0.001
        ft  = 0.001

        # No throughflow
        domain_no = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage)
        for t in domain_no.evolve(yieldstep=yst, finaltime=ft):
            pass
        landInds = (domain_no.centroid_coordinates[:, 0] < 50.).nonzero()[0]
        vol_no = (domain_no.quantities['height'].centroid_values[landInds]
                  * domain_no.areas[landInds]).sum()

        # With throughflow
        riverWall = {'centralWall': [[wallLoc, 0.0, wallHeight],
                                     [wallLoc, 100.0, wallHeight]]}
        riverWall_Par = {'centralWall': {'Cd_through': 0.5}}
        domain_cd = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage,
                                           riverWall=riverWall, riverWall_Par=riverWall_Par)
        for t in domain_cd.evolve(yieldstep=yst, finaltime=ft):
            pass
        vol_cd = (domain_cd.quantities['height'].centroid_values[landInds]
                  * domain_cd.areas[landInds]).sum()

        assert vol_cd > vol_no, \
            f"Throughflow+overtopping should exceed overtopping alone: {vol_cd} vs {vol_no}"

    def test_throughflow_backward_compatible(self):
        """Existing riverwall tests unaffected: Cd_through=0 leaves hydraulic_properties unchanged."""
        wallHeight = 2.0
        InitialOceanStage = 2.5
        InitialLandStage = -999999.

        domain = self.create_domain_DE0(wallHeight, InitialOceanStage, InitialLandStage)
        hp = domain.riverwallData.hydraulic_properties
        names = domain.riverwallData.hydraulic_variable_names
        col = list(names).index('Cd_through')
        # Default must be exactly 0.0
        assert hp[0, col] == 0.0, f"Default Cd_through must be 0.0, got {hp[0, col]}"
        # Should now have 6 columns
        assert hp.shape[1] == 6, f"Expected 6 hydraulic property columns, got {hp.shape[1]}"


class Test_riverwall_notebook(unittest.TestCase):
    """Tests using create_domain_from_regions with breaklines, mirroring the
    docs/source/examples/notebook_create_domain_with_riverwalls.ipynb scenario.

    This exercises the full mesh-from-regions path (as opposed to the
    create_pmesh_from_regions + create_domain_from_file path used in
    Test_riverwall_structure).
    """

    def _make_notebook_domain(self):
        """Reproduce the notebook's 3-wall domain and return it (not evolved)."""
        bounding_polygon = [[0.0, 0.0], [20.0, 0.0], [20.0, 10.0], [0.0, 10.0]]
        boundary_tags = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}
        riverWalls = {
            'wall1': [[5.0,  0.0,  0.5], [5.0,  4.0,  0.5]],
            'wall2': [[15.0, 0.0, -0.5], [15.0, 4.0, -0.5]],
            'wall3': [[10.0, 10.0, 0.0], [10.0, 6.0,  0.0]],
        }
        domain = anuga.create_domain_from_regions(
            bounding_polygon, boundary_tags,
            maximum_triangle_area=0.5,
            breaklines=list(riverWalls.values()))
        domain.set_store(False)
        domain.set_quantity('elevation', lambda x, y: -x / 10, location='centroids')
        domain.set_quantity('friction', 0.01, location='centroids')
        domain.set_quantity('stage', expression='elevation', location='centroids')
        Bi = anuga.Dirichlet_boundary([0.4, 0, 0])
        Bo = anuga.Dirichlet_boundary([-2, 0, 0])
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
        domain.create_riverwalls(riverWalls, verbose=False)
        return domain, riverWalls

    def test_create_domain_from_regions_with_breaklines(self):
        """create_domain_from_regions accepts riverWall dict values as breaklines
        and create_riverwalls registers all three named walls."""
        domain, riverWalls = self._make_notebook_domain()
        # All three walls must be registered
        self.assertEqual(len(domain.riverwallData.names), 3)
        for name in riverWalls:
            self.assertIn(name, domain.riverwallData.names)

    def test_riverwall_crest_heights_stored(self):
        """Crest height (z coordinate) is stored correctly for each wall segment."""
        domain, riverWalls = self._make_notebook_domain()
        # wall1 crest is 0.5 m — all edges on wall1 should have elevation ~0.5
        # wall2 crest is -0.5 m — submerged wall
        # Just check that the riverwall_elevation array is non-empty and finite
        rw_elev = domain.riverwallData.riverwall_elevation
        self.assertGreater(len(rw_elev), 0)
        self.assertTrue(numpy.all(numpy.isfinite(rw_elev)))

    def test_riverwall_edges_lie_on_breaklines(self):
        """All edges flagged as riverwall edges must lie on a breakline (x ≈ 5, 15, or 10)."""
        domain, _ = self._make_notebook_domain()
        edge_inds = domain.riverwallData.riverwall_edges
        self.assertGreater(len(edge_inds), 0)
        edge_x = domain.edge_coordinates[edge_inds, 0]
        # Each edge midpoint must be close to one of the three wall x-positions
        on_wall = (numpy.abs(edge_x - 5.0) < 0.1) | \
                  (numpy.abs(edge_x - 15.0) < 0.1) | \
                  (numpy.abs(edge_x - 10.0) < 0.1)
        self.assertTrue(numpy.all(on_wall),
                        f"Some riverwall edges not near x=5,10,15: {edge_x[~on_wall]}")

    @pytest.mark.slow
    def test_impermeable_wall_blocks_sub_crest_flow(self):
        """Single wall (crest 0.5 m), upstream Dirichlet stage 0.4 m.
        No overtopping should occur so the downstream half stays essentially dry."""
        import tempfile
        import os
        bounding_polygon = [[0, 0], [20, 0], [20, 10], [0, 10]]
        boundary_tags = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}
        riverWalls = {'levee': [[10.0, 0.0, 0.5], [10.0, 10.0, 0.5]]}

        domain = anuga.create_domain_from_regions(
            bounding_polygon, boundary_tags,
            maximum_triangle_area=0.5,
            breaklines=list(riverWalls.values()))
        domain.set_store(False)
        domain.set_quantity('elevation', 0.0, location='centroids')
        domain.set_quantity('friction',  0.01, location='centroids')
        domain.set_quantity('stage',     0.0, location='centroids')
        Bi = anuga.Dirichlet_boundary([0.4, 0.0, 0.0])   # stage < crest
        Bo = anuga.Dirichlet_boundary([-2.0, 0.0, 0.0])
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
        domain.create_riverwalls(riverWalls, verbose=False)

        for t in domain.evolve(yieldstep=5, duration=20):
            pass

        x = domain.centroid_coordinates[:, 0]
        depth = numpy.maximum(
            domain.quantities['stage'].centroid_values -
            domain.quantities['elevation'].centroid_values, 0.0)
        mean_downstream_depth = float(depth[x > 10].mean())
        self.assertLess(mean_downstream_depth, 0.01,
                        f"Downstream should be dry; got mean depth={mean_downstream_depth:.4f} m")

    @pytest.mark.slow
    def test_wall_overtopping_above_crest(self):
        """Single wall (crest 0.5 m), upstream Dirichlet stage 0.8 m.
        Water overtops so the downstream half should become wet."""
        bounding_polygon = [[0, 0], [20, 0], [20, 10], [0, 10]]
        boundary_tags = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}
        riverWalls = {'levee': [[10.0, 0.0, 0.5], [10.0, 10.0, 0.5]]}

        domain = anuga.create_domain_from_regions(
            bounding_polygon, boundary_tags,
            maximum_triangle_area=0.5,
            breaklines=list(riverWalls.values()))
        domain.set_store(False)
        domain.set_quantity('elevation', 0.0, location='centroids')
        domain.set_quantity('friction',  0.01, location='centroids')
        domain.set_quantity('stage',     0.0, location='centroids')
        Bi = anuga.Dirichlet_boundary([0.8, 0.0, 0.0])   # stage > crest
        Bo = anuga.Dirichlet_boundary([-2.0, 0.0, 0.0])
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
        domain.create_riverwalls(riverWalls, verbose=False)

        for t in domain.evolve(yieldstep=5, duration=30):
            pass

        x = domain.centroid_coordinates[:, 0]
        depth = numpy.maximum(
            domain.quantities['stage'].centroid_values -
            domain.quantities['elevation'].centroid_values, 0.0)
        mean_downstream_depth = float(depth[x > 10].mean())
        self.assertGreater(mean_downstream_depth, 0.05,
                           f"Downstream should be wet after overtopping; got depth={mean_downstream_depth:.4f} m")


# =========================================================================
# Tests for the RiverWall runtime interface
# (get/set elevation, get/set hydraulic parameters, get edge coordinates)
# =========================================================================

class Test_riverwall_interface(unittest.TestCase):
    """Unit tests for RiverWall.get/set elevation & hydraulic-parameter methods."""

    def setUp(self):
        """Build a minimal two-wall domain for interface testing."""
        bounding_polygon = [[0., 0.], [0., 20.], [20., 20.], [20., 0.]]
        boundary_tags = {'left': [0], 'top': [1], 'right': [2], 'bottom': [3]}
        riverWalls = {
            'wallA': [[5.0,  0.0, 1.0], [5.0,  20.0, 1.0]],
            'wallB': [[10.0, 0.0, 2.0], [10.0, 20.0, 2.0]],
        }
        anuga.create_pmesh_from_regions(
            bounding_polygon, boundary_tags,
            maximum_triangle_area=4.0,
            breaklines=list(riverWalls.values()),
            use_cache=False, verbose=False,
            filename='test_rw_interface.msh')
        self.domain = anuga.create_domain_from_file('test_rw_interface.msh')
        self.domain.set_flow_algorithm('DE0')
        self.domain.set_name('test_rw_interface')
        self.domain.set_store(False)
        self.domain.set_quantity('elevation', 0.0, location='centroids')
        self.domain.set_quantity('friction', 0.03, location='centroids')
        self.domain.set_quantity('stage', 0.0, location='centroids')
        Br = anuga.Reflective_boundary(self.domain)
        self.domain.set_boundary(
            {'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        self.domain.riverwallData.create_riverwalls(
            riverWalls,
            {'wallA': {'Qfactor': 0.8, 's1': 0.7, 's2': 0.85,
                       'h1': 0.9, 'h2': 1.4, 'Cd_through': 0.1},
             'wallB': {'Qfactor': 1.2}},
            verbose=False)
        self.rw = self.domain.riverwallData

    def tearDown(self):
        for f in ('test_rw_interface.msh', 'test_rw_interface.sww'):
            try:
                os.remove(f)
            except OSError:
                pass

    # --- get_wall_names ---

    def test_get_wall_names_returns_both(self):
        names = self.rw.get_wall_names()
        self.assertIn('wallA', names)
        self.assertIn('wallB', names)
        self.assertEqual(len(names), 2)

    def test_get_wall_names_returns_copy(self):
        names = self.rw.get_wall_names()
        names.append('extra')
        self.assertEqual(len(self.rw.get_wall_names()), 2)

    # --- _name_to_index / _param_to_col error paths ---

    def test_name_to_index_bad_name_raises(self):
        with self.assertRaises(KeyError):
            self.rw._name_to_index('nonexistent')

    def test_param_to_col_bad_param_raises(self):
        with self.assertRaises(KeyError):
            self.rw._param_to_col('bad_param')

    # --- get_edge_coordinates ---

    def test_get_edge_coordinates_shape(self):
        xy = self.rw.get_edge_coordinates('wallA')
        self.assertEqual(xy.ndim, 2)
        self.assertEqual(xy.shape[1], 2)
        self.assertGreater(xy.shape[0], 0)

    def test_get_edge_coordinates_near_wall(self):
        """All returned x-coords should be close to wallA x=5."""
        xy = self.rw.get_edge_coordinates('wallA')
        numpy.testing.assert_allclose(xy[:, 0], 5.0, atol=0.5)

    # --- get_elevation ---

    def test_get_elevation_returns_correct_values(self):
        elev = self.rw.get_elevation('wallA')
        numpy.testing.assert_allclose(elev, 1.0, atol=1e-10)

    def test_get_elevation_returns_copy(self):
        elev = self.rw.get_elevation('wallA')
        elev[:] = 999.0
        elev2 = self.rw.get_elevation('wallA')
        numpy.testing.assert_allclose(elev2, 1.0, atol=1e-10)

    # --- set_elevation (scalar) ---

    def test_set_elevation_scalar_updates_all_edges(self):
        self.rw.set_elevation('wallA', 3.5)
        numpy.testing.assert_allclose(
            self.rw.get_elevation('wallA'), 3.5, atol=1e-10)

    def test_set_elevation_scalar_does_not_affect_other_wall(self):
        self.rw.set_elevation('wallA', 3.5)
        numpy.testing.assert_allclose(
            self.rw.get_elevation('wallB'), 2.0, atol=1e-10)

    # --- set_elevation (array) ---

    def test_set_elevation_array_updates_edges(self):
        n = len(self.rw.get_elevation('wallA'))
        new_vals = numpy.linspace(1.0, 2.0, n)
        self.rw.set_elevation('wallA', new_vals)
        numpy.testing.assert_allclose(
            self.rw.get_elevation('wallA'), new_vals, atol=1e-10)

    def test_set_elevation_wrong_array_length_raises(self):
        with self.assertRaises(ValueError):
            self.rw.set_elevation('wallA', [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])

    # --- set_elevation_offset ---

    def test_set_elevation_offset_adds_to_current(self):
        before = self.rw.get_elevation('wallA').copy()
        self.rw.set_elevation_offset('wallA', 0.5)
        numpy.testing.assert_allclose(
            self.rw.get_elevation('wallA'), before + 0.5, atol=1e-10)

    def test_set_elevation_offset_negative(self):
        before = self.rw.get_elevation('wallB').copy()
        self.rw.set_elevation_offset('wallB', -1.0)
        numpy.testing.assert_allclose(
            self.rw.get_elevation('wallB'), before - 1.0, atol=1e-10)

    # --- get_hydraulic_parameter ---

    def test_get_hydraulic_parameter_qfactor(self):
        val = self.rw.get_hydraulic_parameter('wallA', 'Qfactor')
        self.assertAlmostEqual(val, 0.8)

    def test_get_hydraulic_parameter_default(self):
        val = self.rw.get_hydraulic_parameter('wallB', 'Qfactor')
        self.assertAlmostEqual(val, 1.2)

    def test_get_hydraulic_parameter_cd_through(self):
        val = self.rw.get_hydraulic_parameter('wallA', 'Cd_through')
        self.assertAlmostEqual(val, 0.1)

    # --- set_hydraulic_parameter ---

    def test_set_hydraulic_parameter_updates_value(self):
        self.rw.set_hydraulic_parameter('wallA', 'Qfactor', 0.5)
        self.assertAlmostEqual(
            self.rw.get_hydraulic_parameter('wallA', 'Qfactor'), 0.5)

    def test_set_hydraulic_parameter_does_not_affect_other_wall(self):
        self.rw.set_hydraulic_parameter('wallA', 'Qfactor', 0.5)
        self.assertAlmostEqual(
            self.rw.get_hydraulic_parameter('wallB', 'Qfactor'), 1.2)

    def test_set_hydraulic_parameter_all_params(self):
        for param, val in [('Qfactor', 0.9), ('s1', 0.6), ('s2', 0.8),
                           ('h1', 0.5), ('h2', 1.0), ('Cd_through', 0.3)]:
            self.rw.set_hydraulic_parameter('wallB', param, val)
            self.assertAlmostEqual(
                self.rw.get_hydraulic_parameter('wallB', param), val,
                msg=f"param={param}")

    # --- runtime use: set elevation inside a yield loop ---

    def test_elevation_change_mid_evolve(self):
        """Changing wall elevation mid-simulation should not crash."""
        self.domain.set_quantity('stage', 0.5, location='centroids')
        Br = anuga.Reflective_boundary(self.domain)
        self.domain.set_boundary(
            {'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        for t in self.domain.evolve(yieldstep=5.0, finaltime=10.0):
            if t == 5.0:
                self.rw.set_elevation('wallA', 0.3)
        elev = self.rw.get_elevation('wallA')
        numpy.testing.assert_allclose(elev, 0.3, atol=1e-10)


# =========================================================================
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_riverwall_structure)
    runner = unittest.TextTestRunner()
    runner.run(suite)
