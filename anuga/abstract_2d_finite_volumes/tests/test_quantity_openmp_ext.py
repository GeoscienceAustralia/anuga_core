"""Direct unit tests for the quantity_openmp_ext kernels.

Gradients, extrapolation and the limiters were only exercised through whole
domain evolution before (audit issue #334). These tests drive each kernel on
a fixed four-triangle mesh with values chosen so the expected result can be
worked out by hand from the kernel's definition.

Mesh (points and elements as in test_quantity.py):

    d(0,4)
    |  \\
    b(0,2)---e(2,2)
    |  \\  1  |  \\
    | 0  \\   | 2  \\
    a(0,0)---c(2,0)---f(4,0)

    elements: 0 = bac, 1 = bce, 2 = ecf, 3 = dbe

Triangle 1 has three neighbours (0, 2 and 3); every other triangle has one
neighbour (triangle 1) and two boundary edges. Edge i of a triangle is the
edge opposite its vertex i, and neighbours[k, i] is the triangle across
edge i, or a negative number for a boundary edge.
"""
import numpy as num
import pytest

from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.abstract_2d_finite_volumes import quantity_openmp_ext as ext


@pytest.fixture
def mesh():
    points = [[0.0, 0.0], [0.0, 2.0], [2.0, 0.0], [0.0, 4.0], [2.0, 2.0], [4.0, 0.0]]
    elements = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = Generic_Domain(points, elements)
    domain.check_integrity()
    domain.beta_w = 1.0
    return domain


def _linear(x, y):
    return 2.0 * x + 3.0 * y


def _quantity(mesh, centroid_values=None):
    q = Quantity(mesh)
    if centroid_values is not None:
        q.centroid_values[:] = centroid_values
    return q


# ---------------------------------------------------------------------------
# Gradients
# ---------------------------------------------------------------------------

def test_compute_local_gradients_recovers_a_plane_from_vertex_values(mesh):
    q = _quantity(mesh)
    xy = mesh.vertex_coordinates            # (3N, 2)
    q.vertex_values[:] = _linear(xy[:, 0], xy[:, 1]).reshape(-1, 3)

    ext.compute_local_gradients(q)

    assert num.allclose(q.x_gradient, 2.0)
    assert num.allclose(q.y_gradient, 3.0)


def test_compute_gradients_is_exact_on_an_interior_triangle(mesh):
    c = mesh.centroid_coordinates
    q = _quantity(mesh, _linear(c[:, 0], c[:, 1]))

    ext.compute_gradients(q)

    # Triangle 1 has three true neighbours, so the three-point fit through
    # their centroids recovers the plane exactly.
    assert mesh.number_of_boundaries[1] == 0
    assert abs(q.x_gradient[1] - 2.0) < 1e-12
    assert abs(q.y_gradient[1] - 3.0) < 1e-12


def test_compute_gradients_projects_onto_the_neighbour_direction_with_one_neighbour(mesh):
    c = mesh.centroid_coordinates
    q = _quantity(mesh, _linear(c[:, 0], c[:, 1]))

    ext.compute_gradients(q)

    # With a single neighbour the kernel uses _gradient2: the gradient along
    # the centroid-to-centroid direction, and zero across it.
    for k in (0, 2, 3):
        assert mesh.number_of_boundaries[k] == 2
        n = int(mesh.surrogate_neighbours[k][mesh.surrogate_neighbours[k] != k][0])
        dx = c[n] - c[k]
        dq = q.centroid_values[n] - q.centroid_values[k]
        expected = dx * dq / (dx @ dx)
        assert num.allclose([q.x_gradient[k], q.y_gradient[k]], expected)


# ---------------------------------------------------------------------------
# Extrapolation
# ---------------------------------------------------------------------------

def test_extrapolate_from_gradient_is_the_plane_through_the_centroid(mesh):
    q = _quantity(mesh, [1.0, 2.0, 3.0, 4.0])
    q.x_gradient[:] = [1.0, 0.0, -1.0, 0.5]
    q.y_gradient[:] = [0.0, 1.0, 2.0, -0.5]

    ext.extrapolate_from_gradient(q)

    c = mesh.centroid_coordinates
    xy = mesh.vertex_coordinates.reshape(-1, 3, 2)
    for k in range(4):
        expected = (q.centroid_values[k]
                    + q.x_gradient[k] * (xy[k, :, 0] - c[k, 0])
                    + q.y_gradient[k] * (xy[k, :, 1] - c[k, 1]))
        assert num.allclose(q.vertex_values[k], expected)
        # Edge i is the midpoint of the two vertices other than i.
        v = q.vertex_values[k]
        assert num.allclose(q.edge_values[k], [0.5 * (v[1] + v[2]),
                                               0.5 * (v[2] + v[0]),
                                               0.5 * (v[0] + v[1])])
        # A plane through the centroid averages back to the centroid value.
        assert abs(q.vertex_values[k].mean() - q.centroid_values[k]) < 1e-12


# ---------------------------------------------------------------------------
# Limiters. All of them scale the deviations of the extrapolated values from
# the centroid by one factor phi per triangle, phi = min over the three
# values of min(beta * r, 1), where r is how far the deviation may go before
# it leaves the admissible range.
# ---------------------------------------------------------------------------

def _set_triangle_1(q, vertex=None, edge=None):
    """Triangle 1 sits between values 0.5, 1.5 and 1.0, so its admissible
    range from all neighbours is [0.5, 1.5] around its own value 1.0."""
    q.centroid_values[:] = [0.5, 1.0, 1.5, 1.0]
    if vertex is not None:
        q.vertex_values[1] = vertex
    if edge is not None:
        q.edge_values[1] = edge
    q.x_gradient[:] = 1.0
    q.y_gradient[:] = 1.0


def test_limit_vertices_by_all_neighbours(mesh):
    q = _quantity(mesh)
    # Deviations +1.0, -0.25, 0: the first overshoots qmax by 2x (r = 0.5),
    # the second stays inside (r = 2), so phi = 0.5.
    _set_triangle_1(q, vertex=[2.0, 0.75, 1.0])

    ext.limit_vertices_by_all_neighbours(q)

    assert num.allclose(q.vertex_values[1], [1.5, 0.875, 1.0])
    assert q.x_gradient[1] == 0.5 and q.y_gradient[1] == 0.5
    v = q.vertex_values[1]
    assert num.allclose(q.edge_values[1], [0.5 * (v[1] + v[2]),
                                           0.5 * (v[2] + v[0]),
                                           0.5 * (v[0] + v[1])])


def test_limit_vertices_by_all_neighbours_leaves_admissible_values_alone(mesh):
    q = _quantity(mesh)
    _set_triangle_1(q, vertex=[1.4, 0.6, 1.0])

    ext.limit_vertices_by_all_neighbours(q)

    assert num.allclose(q.vertex_values[1], [1.4, 0.6, 1.0])
    assert q.x_gradient[1] == 1.0


def test_limit_vertices_by_all_neighbours_honours_beta(mesh):
    q = _quantity(mesh)
    mesh.beta_w = 0.5
    _set_triangle_1(q, vertex=[1.4, 0.6, 1.0])   # r = 1.25, 1.25 and (dq = 0) 1

    ext.limit_vertices_by_all_neighbours(q)

    # phi = min over vertices of min(beta * r, 1). A vertex with no deviation
    # keeps r = 1, so with beta < 1 the limiter never lets phi exceed beta:
    # here 0.5, not the 0.625 the two deviating vertices alone would allow.
    assert num.allclose(q.vertex_values[1], 1.0 + 0.5 * num.array([0.4, -0.4, 0.0]))
    assert q.x_gradient[1] == 0.5


def test_limit_edges_by_all_neighbours(mesh):
    q = _quantity(mesh)
    _set_triangle_1(q, edge=[2.0, 0.75, 1.0])    # same deviations, at the edges

    ext.limit_edges_by_all_neighbours(q)

    assert num.allclose(q.edge_values[1], [1.5, 0.875, 1.0])
    assert q.x_gradient[1] == 0.5
    # Vertices are reconstructed from the edges: v_i = e_j + e_k - e_i.
    e = q.edge_values[1]
    assert num.allclose(q.vertex_values[1], [e[1] + e[2] - e[0],
                                             e[2] + e[0] - e[1],
                                             e[0] + e[1] - e[2]])


def test_limit_edges_by_neighbour_bounds_each_edge_by_the_triangle_across_it(mesh):
    q = _quantity(mesh)
    q.centroid_values[:] = [0.5, 1.0, 1.5, 1.0]
    # Edge i of triangle 1 faces neighbours[1, i]. Push every edge up by 0.4:
    # the edge facing the 1.5 neighbour is fine (r = 1.25), the edge facing
    # the 0.5 neighbour has qmax = qc, so r = 0 and phi collapses to 0.
    q.edge_values[1] = 1.4

    ext.limit_edges_by_neighbour(q)

    assert num.allclose(q.edge_values[1], 1.0)
    assert num.allclose(q.vertex_values[1], 1.0)


def test_limit_edges_by_neighbour_treats_a_boundary_edge_as_its_own_value(mesh):
    q = _quantity(mesh)
    q.centroid_values[:] = [1.0, 1.0, 1.0, 1.0]
    # Triangle 0 has two boundary edges. On a boundary edge the neighbour
    # value is the centroid itself, so any deviation there gives r = 0.
    boundary_edges = num.where(mesh.neighbours[0] < 0)[0]
    assert len(boundary_edges) == 2
    q.edge_values[0] = 1.0
    q.edge_values[0, boundary_edges[0]] = 1.3

    ext.limit_edges_by_neighbour(q)

    assert num.allclose(q.edge_values[0], 1.0)


def test_limit_gradient_by_neighbour_ignores_boundary_edges(mesh):
    q = _quantity(mesh)
    q.centroid_values[:] = [1.0, 1.0, 1.0, 1.0]
    q.x_gradient[:] = 1.0
    q.y_gradient[:] = 1.0
    boundary_edges = num.where(mesh.neighbours[0] < 0)[0]
    interior_edge = num.where(mesh.neighbours[0] >= 0)[0][0]
    q.edge_values[0] = 1.0
    q.edge_values[0, boundary_edges[0]] = 1.3     # unconstrained: no neighbour
    q.edge_values[0, interior_edge] = 1.0          # dq = 0: no constraint either

    ext.limit_gradient_by_neighbour(q)

    # Unlike limit_edges_by_neighbour, a boundary edge imposes nothing, so
    # phi stays 1 and the deviation survives.
    assert q.edge_values[0, boundary_edges[0]] == pytest.approx(1.3)
    assert q.x_gradient[0] == 1.0


def test_bound_vertices_below_by_constant(mesh):
    q = _quantity(mesh, [1.0, 1.0, 1.0, 1.0])
    q.vertex_values[:] = 1.0          # start with no deviation anywhere
    q.x_gradient[:] = 1.0
    q.y_gradient[:] = 1.0
    # Deviations -0.8, +0.2, 0 with a floor of 0.6: the first may only drop
    # to 0.6, so r = 0.5 and everything is scaled by 0.5.
    q.vertex_values[0] = [0.2, 1.2, 1.0]

    ext.bound_vertices_below_by_constant(q, 0.6)

    assert num.allclose(q.vertex_values[0], [0.6, 1.1, 1.0])
    assert q.x_gradient[0] == 0.5
    # Only downward deviations are constrained: a triangle above the floor
    # is untouched.
    assert num.allclose(q.vertex_values[1], q.centroid_values[1])


def test_bound_vertices_below_by_quantity(mesh):
    q = _quantity(mesh, [1.0, 1.0, 1.0, 1.0])
    q.vertex_values[:] = 1.0
    q.x_gradient[:] = 1.0
    q.y_gradient[:] = 1.0
    q.vertex_values[0] = [0.2, 1.2, 1.0]
    floor = _quantity(mesh)
    floor.vertex_values[:] = 0.0
    floor.vertex_values[0] = [0.6, 0.0, 0.0]      # per-vertex floor

    ext.bound_vertices_below_by_quantity(q, floor)

    assert num.allclose(q.vertex_values[0], [0.6, 1.1, 1.0])
    assert q.x_gradient[0] == 0.5
