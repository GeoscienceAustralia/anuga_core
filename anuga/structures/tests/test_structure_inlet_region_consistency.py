"""The parallel and sequential structure operators must select the same
inlet triangles (anuga-community/anuga_core#367).

The sequential Structure_operator builds each inlet from the exchange line
extended back by the apron; the parallel one computed the same polygon and
then handed its inlet the bare line. The factory returns the parallel class
whenever mpi4py imports, serial domain or not, so the exchange region of the
same script depended on whether mpi4py was installed.
"""
import numpy as np
import pytest

import anuga
from anuga import rectangular_cross_domain, Reflective_boundary
from anuga.structures.boyd_box_operator import Boyd_box_operator as Sequential_boyd

try:
    import mpi4py  # noqa: F401
    from anuga.parallel.parallel_operator_factory import Boyd_box_operator as Factory_boyd
    from anuga.parallel.parallel_boyd_box_operator import Parallel_Boyd_box_operator
    HAVE_MPI4PY = True
except ImportError:
    HAVE_MPI4PY = False


def make_domain():
    d = rectangular_cross_domain(40, 20, len1=200.0, len2=50.0)
    d.set_flow_algorithm('DE0')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    return d


def inlet_sets(op):
    return [sorted(int(k) for k in inlet.triangle_indices) for inlet in op.inlets]


CULVERT = dict(end_points=[[60.0, 25.0], [140.0, 25.0]], width=10.0, height=1.0,
               apron=5.0, losses=1.5, manning=0.013, verbose=False)
SKEW = dict(exchange_lines=[[[60.0, 20.0], [60.0, 30.0]], [[140.0, 30.0], [140.0, 20.0]]],
            width=10.0, height=1.0, apron=5.0, losses=1.5, manning=0.013, verbose=False)


def test_sequential_inlet_is_the_line_plus_the_apron():
    d = make_domain()
    op = Sequential_boyd(d, **CULVERT)
    # 5 m x 2.5 m rectangles of four 3.125 m^2 triangles; the apron polygon
    # is 10 m x 5 m (16 triangles inside) and its four edges lie on mesh
    # lines, adding the triangles along them on the far side
    for inlet in op.inlets:
        assert inlet.region.type == 'polygon'
        assert np.isclose(inlet.area, 26 * 3.125), inlet.area


@pytest.mark.skipif(not HAVE_MPI4PY, reason='needs mpi4py for the parallel operator')
def test_factory_operator_selects_the_same_inlet_triangles_as_the_sequential_one():
    for kwargs in (CULVERT, SKEW):
        seq = Sequential_boyd(make_domain(), **kwargs)
        fac = Factory_boyd(make_domain(), **kwargs)
        assert isinstance(fac, Parallel_Boyd_box_operator)
        assert inlet_sets(fac) == inlet_sets(seq), kwargs
        for a, b in zip(fac.inlets, seq.inlets):
            assert np.isclose(a.area, b.area)


@pytest.mark.skipif(not HAVE_MPI4PY, reason='needs mpi4py for the parallel operator')
def test_no_apron_means_the_line_alone_in_both():
    kwargs = dict(CULVERT, apron=None)
    seq = Sequential_boyd(make_domain(), **kwargs)
    fac = Factory_boyd(make_domain(), **kwargs)
    # apron=None is replaced by the width, so both still build polygons;
    # the point is only that they agree
    assert inlet_sets(fac) == inlet_sets(seq)
