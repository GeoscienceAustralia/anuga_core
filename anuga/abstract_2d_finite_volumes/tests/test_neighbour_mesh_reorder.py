#!/usr/bin/env python



#FIXME: Seperate the tests for mesh and general_mesh

#FIXME (Ole): Maxe this test independent of anything that inherits from General_mesh (namely shallow_water)
import unittest
from math import sqrt
import copy


from anuga.abstract_2d_finite_volumes.neighbour_mesh import *
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, rectangular
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_periodic
from anuga.config import epsilon

from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geometry.polygon import is_inside_polygon
from anuga.utilities.numerical_tools import ensure_numeric

import numpy as num

def compare_meshes(mesh_1, mesh_2):

    # test that reorder is consistent
    num.testing.assert_allclose(mesh_1.nodes, mesh_2.nodes)
    num.testing.assert_allclose(mesh_1.triangles, mesh_2.triangles)

    assert mesh_1.boundary == mesh_2.boundary

    assert mesh_1.tagged_elements.keys() == mesh_2.tagged_elements.keys()
    for k in mesh_1.tagged_elements:
        num.testing.assert_array_equal(mesh_1.tagged_elements[k], mesh_2.tagged_elements[k])


    num.testing.assert_allclose(mesh_1.number_of_triangles_per_node, mesh_2.number_of_triangles_per_node)
    num.testing.assert_allclose(mesh_1.node_index, mesh_2.node_index)

    def get_triangle_set(mesh):
        triangle_set = set()
        for i in range(mesh.number_of_nodes):
            first = mesh.node_index[i]
            count = mesh.number_of_triangles_per_node[i]
            for index in mesh.vertex_value_indices[first:first+count]:
                triangle_set.add((i, index // 3, index % 3))
        return triangle_set


    mesh_1_triangle_set = get_triangle_set(mesh_1)
    mesh_2_triangle_set = get_triangle_set(mesh_2)

    assert mesh_1_triangle_set == mesh_2_triangle_set

    num.testing.assert_allclose(mesh_1.neighbours, mesh_2.neighbours)
    num.testing.assert_allclose(mesh_1.neighbour_edges, mesh_2.neighbour_edges)
    num.testing.assert_allclose(mesh_1.number_of_boundaries, mesh_2.number_of_boundaries)

    num.testing.assert_allclose(mesh_1.surrogate_neighbours, mesh_2.surrogate_neighbours)

    num.testing.assert_allclose(mesh_1.vertex_coordinates, mesh_2.vertex_coordinates)
    num.testing.assert_allclose(mesh_1.edge_midpoint_coordinates, mesh_2.edge_midpoint_coordinates)
    num.testing.assert_allclose(mesh_1.normals, mesh_2.normals)
    num.testing.assert_allclose(mesh_1.areas, mesh_2.areas)
    num.testing.assert_allclose(mesh_1.edgelengths, mesh_2.edgelengths)
    num.testing.assert_allclose(mesh_1.centroid_coordinates, mesh_2.centroid_coordinates)
    num.testing.assert_allclose(mesh_1.radii, mesh_2.radii)



class Test_Mesh_Reorder(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_reorder(self):
        """test_reorder(self):

        Check that reordering works

        """

        # Build test mesh

        from numpy import array

        nodes = array([[0., 0.],
        [0., 1.],
        [1., 0.],
        [1., 1.],
        [2., 0.],
        [2., 1.]])

        triangles = array([[2, 3, 0],
        [1, 0, 3],
        [4, 5, 2],
        [3, 2, 5]])

        boundary = {(0, 1): 'bottom',
        (1, 1): 'top',
        (1, 2): 'left',
        (2, 1): 'bottom',
        (2, 2): 'right',
        (3, 1): 'top'}

        tagged_elements = {"south" : [0,2]}

        mesh = Mesh(nodes, triangles,
                    boundary = boundary,
                    tagged_elements=tagged_elements)

        new_order = [2, 3, 0, 1]   # new_index = new_order[old_index]

        new_nodes = array([[0., 0.],
        [0., 1.],
        [1., 0.],
        [1., 1.],
        [2., 0.],
        [2., 1.]])

        new_triangles = array([[4, 5, 2],
        [3, 2, 5],
        [2, 3, 0],
        [1, 0, 3]])

        new_boundary = {(2, 1): 'bottom',
        (3, 1): 'top',
        (3, 2): 'left',
        (0, 1): 'bottom',
        (0, 2): 'right',
        (1, 1): 'top'}

        tagged_elements = {"south" : [2,0]}


        # Create a new_mesh with correctly reordered input
        new_mesh = Mesh(new_nodes, new_triangles,
                        boundary = new_boundary,
                        tagged_elements=tagged_elements)



        reorder_mesh = mesh.reorder(new_order, in_place = True)

        compare_meshes(new_mesh, reorder_mesh)



    def test_reorder_larger(self):
        """test_reorder(self):

        Check that reordering works

        """

        # Build test mesh

        from numpy import array

        nodes, triangles, boundary = rectangular(2,3, 2,3)

        tagged_elements = {"south" : [1,0,7,6],
                           "north" : [5,4,11,10]}

        mesh = Mesh(nodes, triangles, boundary, tagged_elements=tagged_elements)

        N = len(mesh)
        seed = 42

        import random

        rng = random.Random(seed)          # independent RNG
        new_order = list(range(N))
        rng.shuffle(new_order)  # new_index = new_order[old_index]

        new_order = num.array(new_order)

        inv_order = num.empty_like(new_order)
        inv_order[new_order] = num.arange(new_order.size)

        new_nodes = nodes.copy()
        new_triangles = triangles.copy()[new_order]
        new_boundary = {(int(inv_order[i]), j): v for (i, j), v in boundary.items()}

        new_tagged_elements = {k : num.array([inv_order[i] for i in v]) for k, v in tagged_elements.items()}
        #new_boundary = dict(sorted(new_boundary.items()))

        # Create a new_mesh with reordered input
        new_mesh = Mesh(new_nodes, new_triangles,
                        boundary = new_boundary,
                        tagged_elements=new_tagged_elements)

        # from original mesh create a reordered mesh
        reorder_mesh = mesh.reorder(new_order, in_place = True)

        compare_meshes(new_mesh, reorder_mesh)

    def test_basic_mesh_reorder_neighbours_consistent(self):
        """Basic_mesh.reorder() must produce correct neighbours even when
        _neighbours has not been accessed before the reorder call.

        Regression test for a bug where the lazy-built _neighbours were
        constructed from an unreordered _triangle_neighbours cache after
        reorder(), causing ghost-layer BFS to follow wrong adjacency and
        produce more ghost triangles than distribute() for the same mesh.
        """
        import numpy as np
        from anuga.abstract_2d_finite_volumes.basic_mesh import (
            rectangular_cross_basic_mesh)

        # Build a small mesh.  _neighbours starts as None (not yet accessed).
        bm = rectangular_cross_basic_mesh(5, 5, len1=5.0, len2=5.0)
        assert bm._neighbours is None, \
            "Precondition: _neighbours should be None before first access"

        # Capture ground-truth neighbours BEFORE reordering.
        nbrs_before = bm.neighbours.copy()   # triggers _build_neighbours

        # Reset so we can test the lazy path.
        bm2 = rectangular_cross_basic_mesh(5, 5, len1=5.0, len2=5.0)
        assert bm2._neighbours is None

        # Apply a non-trivial permutation.
        N = bm2.number_of_triangles
        rng = np.random.default_rng(42)
        new_order = rng.permutation(N)

        reordered = bm2.reorder(new_order, in_place=False)

        # After reorder the neighbours must be self-consistent:
        # for each triangle i and each edge j, if reordered.neighbours[i,j] == k
        # then triangle k must be adjacent to triangle i.
        nbrs = reordered.neighbours
        tris = reordered.triangles
        for i in range(N):
            for j in range(3):
                k = nbrs[i, j]
                if k < 0:
                    continue  # boundary edge, skip
                # Triangles i and k must share exactly 2 nodes.
                shared = set(tris[i]) & set(tris[k])
                assert len(shared) == 2, (
                    f"Triangle {i} and neighbour {k} share {len(shared)} "
                    f"nodes (expected 2) after reorder — neighbours are stale")

    # def test_reorder_larger_16_16(self):
    #     """Test larger mesh which failed in sequential_dist example

    #     """

    #     # Build test mesh

    #     from anuga import rectangular_cross_domain

    #     N = 29
    #     M = 29
    #     verbose = True

    #     domain = rectangular_cross_domain(N, M, 1.0, 1.0, verbose=verbose)
    #     mesh = domain.mesh



#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Mesh)
    runner = unittest.TextTestRunner()#verbosity=2)
    runner.run(suite)
