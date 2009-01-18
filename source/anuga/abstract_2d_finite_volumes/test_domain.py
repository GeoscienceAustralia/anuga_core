#!/usr/bin/env python

import unittest
from math import sqrt

from anuga.abstract_2d_finite_volumes.domain import *
from anuga.config import epsilon

import Numeric as num


def add_to_verts(tag, elements, domain):
    if tag == "mound":
        domain.test = "Mound"


def set_bottom_friction(tag, elements, domain):
    if tag == "bottom":
        #print 'bottom - indices',elements
        domain.set_quantity('friction', 0.09, indices = elements)

def set_top_friction(tag, elements, domain):
    if tag == "top":
        #print 'top - indices',elements
        domain.set_quantity('friction', 1., indices = elements)


def set_all_friction(tag, elements, domain):
    if tag == 'all':
        new_values = domain.get_quantity('friction').get_values(indices = elements) + 10.0

        domain.set_quantity('friction', new_values, indices = elements)


class Test_Domain(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_simple(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
        other_quantities = ['elevation', 'friction']

        domain = Domain(points, vertices, None,
                        conserved_quantities, other_quantities)
        domain.check_integrity()

        for name in conserved_quantities + other_quantities:
            assert domain.quantities.has_key(name)


        assert domain.get_conserved_quantities(0, edge=1) == 0.


    def test_conserved_quantities(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'])


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        domain.set_quantity('xmomentum', [[1,2,3], [5,5,5],
                                          [0,0,9], [-6, 3, 3]])

        domain.check_integrity()

        #Centroids
        q = domain.get_conserved_quantities(0)
        assert num.allclose(q, [2., 2., 0.])

        q = domain.get_conserved_quantities(1)
        assert num.allclose(q, [5., 5., 0.])

        q = domain.get_conserved_quantities(2)
        assert num.allclose(q, [3., 3., 0.])

        q = domain.get_conserved_quantities(3)
        assert num.allclose(q, [0., 0., 0.])


        #Edges
        q = domain.get_conserved_quantities(0, edge=0)
        assert num.allclose(q, [2.5, 2.5, 0.])
        q = domain.get_conserved_quantities(0, edge=1)
        assert num.allclose(q, [2., 2., 0.])
        q = domain.get_conserved_quantities(0, edge=2)
        assert num.allclose(q, [1.5, 1.5, 0.])

        for i in range(3):
            q = domain.get_conserved_quantities(1, edge=i)
            assert num.allclose(q, [5, 5, 0.])


        q = domain.get_conserved_quantities(2, edge=0)
        assert num.allclose(q, [4.5, 4.5, 0.])
        q = domain.get_conserved_quantities(2, edge=1)
        assert num.allclose(q, [4.5, 4.5, 0.])
        q = domain.get_conserved_quantities(2, edge=2)
        assert num.allclose(q, [0., 0., 0.])


        q = domain.get_conserved_quantities(3, edge=0)
        assert num.allclose(q, [3., 3., 0.])
        q = domain.get_conserved_quantities(3, edge=1)
        assert num.allclose(q, [-1.5, -1.5, 0.])
        q = domain.get_conserved_quantities(3, edge=2)
        assert num.allclose(q, [-1.5, -1.5, 0.])



    def test_create_quantity_from_expression(self):
        """Quantity created from other quantities using arbitrary expression

        """


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'],
                        other_quantities = ['elevation', 'friction'])


        domain.set_quantity('elevation', -1)


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        domain.set_quantity('xmomentum', [[1,2,3], [5,5,5],
                                          [0,0,9], [-6, 3, 3]])

        domain.set_quantity('ymomentum', [[3,3,3], [4,2,1],
                                          [2,4,-1], [1, 0, 1]])

        domain.check_integrity()



        expression = 'stage - elevation'
        Q = domain.create_quantity_from_expression(expression)

        assert num.allclose(Q.vertex_values, [[2,3,4], [6,6,6],
                                              [1,1,10], [-5, 4, 4]])

        expression = '(xmomentum*xmomentum + ymomentum*ymomentum)**0.5'
        Q = domain.create_quantity_from_expression(expression)

        X = domain.quantities['xmomentum'].vertex_values
        Y = domain.quantities['ymomentum'].vertex_values

        assert num.allclose(Q.vertex_values, (X**2 + Y**2)**0.5)



    def test_set_quanitities_to_be_monitored(self):
        """test_set_quanitities_to_be_monitored
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'],
                        other_quantities = ['elevation', 'friction', 'depth'])


        assert domain.quantities_to_be_monitored is None
        domain.set_quantities_to_be_monitored(['stage', 'stage-elevation'])
        assert len(domain.quantities_to_be_monitored) == 2
        assert domain.quantities_to_be_monitored.has_key('stage')
        assert domain.quantities_to_be_monitored.has_key('stage-elevation')
        for key in domain.quantities_to_be_monitored['stage'].keys():
            assert domain.quantities_to_be_monitored['stage'][key] is None


        # Check that invalid requests are dealt with
        try:
            domain.set_quantities_to_be_monitored(['yyyyy'])        
        except:
            pass
        else:
            msg = 'Should have caught illegal quantity'
            raise Exception, msg

        try:
            domain.set_quantities_to_be_monitored(['stage-xx'])        
        except NameError:
            pass
        else:
            msg = 'Should have caught illegal quantity'
            raise Exception, msg

        try:
            domain.set_quantities_to_be_monitored('stage', 'stage-elevation')
        except:
            pass
        else:
            msg = 'Should have caught too many arguments'
            raise Exception, msg

        try:
            domain.set_quantities_to_be_monitored('stage', 'blablabla')
        except:
            pass
        else:
            msg = 'Should have caught polygon as a string'
            raise Exception, msg        



        # Now try with a polygon restriction
        domain.set_quantities_to_be_monitored('xmomentum',
                                              polygon=[[1,1], [1,3], [3,3], [3,1]],
                                              time_interval = [0,3])
        assert domain.monitor_indices[0] == 1
        assert domain.monitor_time_interval[0] == 0
        assert domain.monitor_time_interval[1] == 3        
        

    def test_set_quantity_from_expression(self):
        """Quantity set using arbitrary expression

        """


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'],
                        other_quantities = ['elevation', 'friction', 'depth'])


        domain.set_quantity('elevation', -1)


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        domain.set_quantity('xmomentum', [[1,2,3], [5,5,5],
                                          [0,0,9], [-6, 3, 3]])

        domain.set_quantity('ymomentum', [[3,3,3], [4,2,1],
                                          [2,4,-1], [1, 0, 1]])




        domain.set_quantity('depth', expression = 'stage - elevation')

        domain.check_integrity()




        Q = domain.quantities['depth']

        assert num.allclose(Q.vertex_values, [[2,3,4], [6,6,6],
                                              [1,1,10], [-5, 4, 4]])



                                      
    def test_add_quantity(self):
        """Test that quantities already set can be added to using
        add_quantity

        """


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'],
                        other_quantities = ['elevation', 'friction', 'depth'])


        A = num.array([[1,2,3], [5,5,-5], [0,0,9], [-6,3,3]], num.Float)
        B = num.array([[2,4,4], [3,2,1], [6,-3,4], [4,5,-1]], num.Float)
        
        #print A
        #print B
        #print A+B
        
        # Shorthands
        stage = domain.quantities['stage']
        elevation = domain.quantities['elevation']
        depth = domain.quantities['depth']
        
        # Go testing
        domain.set_quantity('elevation', A)
        domain.add_quantity('elevation', B)
        assert num.allclose(elevation.vertex_values, A+B)
        
        domain.add_quantity('elevation', 4)
        assert num.allclose(elevation.vertex_values, A+B+4)        
        
        
        # Test using expression
        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])        
        domain.set_quantity('depth', 1.0)                                     
        domain.add_quantity('depth', expression = 'stage - elevation')        
        assert num.allclose(depth.vertex_values, stage.vertex_values-elevation.vertex_values+1)
                
        
        # Check self referential expression
        reference = 2*stage.vertex_values - depth.vertex_values 
        domain.add_quantity('stage', expression = 'stage - depth')                
        assert num.allclose(stage.vertex_values, reference)        
                                      

        # Test using a function
        def f(x, y):
            return x+y
            
        domain.set_quantity('elevation', f)            
        domain.set_quantity('stage', 5.0)
        domain.set_quantity('depth', expression = 'stage - elevation')
        
        domain.add_quantity('depth', f)
        assert num.allclose(stage.vertex_values, depth.vertex_values)                
         
            
        
        
                                      
                                      
    def test_setting_timestepping_method(self):
        """test_set_quanitities_to_be_monitored
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe, daf, dae
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]


        domain = Domain(points, vertices, boundary=None,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'],
                        other_quantities = ['elevation', 'friction', 'depth'])


        domain.timestepping_method = None


        # Check that invalid requests are dealt with
        try:
            domain.set_timestepping_method('eee')        
        except:
            pass
        else:
            msg = 'Should have caught illegal method'
            raise Exception, msg


        #Should have no trouble with euler, rk2 or rk3
        domain.set_timestepping_method('euler')
        domain.set_timestepping_method('rk2')
        domain.set_timestepping_method('rk3')

        #test get timestepping method
        assert domain.get_timestepping_method() == 'rk3'



    def test_boundary_indices(self):

        from anuga.config import default_boundary_tag


        a = [0.0, 0.5]
        b = [0.0, 0.0]
        c = [0.5, 0.5]

        points = [a, b, c]
        vertices = [ [0,1,2] ]
        domain = Domain(points, vertices)

        domain.set_boundary( {default_boundary_tag: Dirichlet_boundary([5,2,1])} )


        domain.check_integrity()

        assert num.allclose(domain.neighbours, [[-1,-2,-3]])



    def test_boundary_conditions(self):

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'First',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Second'}


        domain = Domain(points, vertices, boundary,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'])
        domain.check_integrity()



        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])


        domain.set_boundary( {'First': Dirichlet_boundary([5,2,1]),
                              'Second': Transmissive_boundary(domain)} )

        domain.update_boundary()

        assert domain.quantities['stage'].boundary_values[0] == 5. #Dirichlet
        assert domain.quantities['stage'].boundary_values[1] == 5. #Dirichlet
        assert domain.quantities['stage'].boundary_values[2] ==\
               domain.get_conserved_quantities(2, edge=0)[0] #Transmissive (4.5)
        assert domain.quantities['stage'].boundary_values[3] ==\
               domain.get_conserved_quantities(2, edge=1)[0] #Transmissive (4.5)
        assert domain.quantities['stage'].boundary_values[4] ==\
               domain.get_conserved_quantities(3, edge=1)[0] #Transmissive (-1.5)
        assert domain.quantities['stage'].boundary_values[5] ==\
               domain.get_conserved_quantities(3, edge=2)[0] #Transmissive (-1.5)

        #Check enumeration
        for k, ((vol_id, edge_id), _) in enumerate(domain.boundary_objects):
            assert domain.neighbours[vol_id, edge_id] == -k-1




    def test_distribute_first_order(self):
        """Domain implements a default first order gradient limiter
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third'}


        domain = Domain(points, vertices, boundary,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'])
        domain.check_integrity()


        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        assert num.allclose( domain.quantities['stage'].centroid_values,
                             [2,5,3,0] )

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2],
                                          [3,3,3], [4, 4, 4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40, 40, 40]])


        domain.distribute_to_vertices_and_edges()

        #First order extrapolation
        assert num.allclose( domain.quantities['stage'].vertex_values,
                             [[ 2.,  2.,  2.],
                              [ 5.,  5.,  5.],
                              [ 3.,  3.,  3.],
                              [ 0.,  0.,  0.]])




    def test_update_conserved_quantities(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third'}


        domain = Domain(points, vertices, boundary,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'])
        domain.check_integrity()


        domain.set_quantity('stage', [1,2,3,4], location='centroids')
        domain.set_quantity('xmomentum', [1,2,3,4], location='centroids')
        domain.set_quantity('ymomentum', [1,2,3,4], location='centroids')


        #Assign some values to update vectors
        #Set explicit_update

        for name in domain.conserved_quantities:
            domain.quantities[name].explicit_update = num.array([4.,3.,2.,1.])
            domain.quantities[name].semi_implicit_update = num.array([1.,1.,1.,1.])


        #Update with given timestep (assuming no other forcing terms)
        domain.timestep = 0.1
        domain.update_conserved_quantities()

        sem = num.array([1.,1.,1.,1.])/num.array([1, 2, 3, 4])
        denom = num.ones(4, num.Float)-domain.timestep*sem

#        x = array([1, 2, 3, 4]) + array( [.4,.3,.2,.1] )
#        x /= denom

        x = num.array([1., 2., 3., 4.])
        x /= denom
        x += domain.timestep*num.array( [4,3,2,1] )

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].centroid_values, x)


    def test_set_region(self):
        """Set quantities for sub region
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = { (0, 0): 'Third',
                     (0, 2): 'First',
                     (2, 0): 'Second',
                     (2, 1): 'Second',
                     (3, 1): 'Second',
                     (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary,
                        conserved_quantities =\
                        ['stage', 'xmomentum', 'ymomentum'])
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5],
                                      [0,0,9], [-6, 3, 3]])

        assert num.allclose( domain.quantities['stage'].centroid_values,
                             [2,5,3,0] )

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2],
                                          [3,3,3], [4, 4, 4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40, 40, 40]])


        domain.distribute_to_vertices_and_edges()

        #First order extrapolation
        assert num.allclose( domain.quantities['stage'].vertex_values,
                             [[ 2.,  2.,  2.],
                              [ 5.,  5.,  5.],
                              [ 3.,  3.,  3.],
                              [ 0.,  0.,  0.]])

        domain.build_tagged_elements_dictionary({'mound':[0,1]})
        domain.set_region([add_to_verts])

        self.failUnless(domain.test == "Mound",
                        'set region failed')


    def test_track_speeds(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.timestepping_statistics(track_speeds=True)



    def test_region_tags(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        domain.set_region([set_bottom_friction, set_top_friction])
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),\
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]])

        domain.set_region([set_all_friction])
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])


    def test_region_tags2(self):
        """
        get values based on triangle lists.
        """
        from mesh_factory import rectangular
        from shallow_water import Domain

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        domain.set_region('top', 'friction', 1.0)
        domain.set_region('bottom', 'friction', 0.09)
        
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]])
        
        domain.set_region([set_bottom_friction, set_top_friction])
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]])

        domain.set_region([set_all_friction])
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Domain,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
