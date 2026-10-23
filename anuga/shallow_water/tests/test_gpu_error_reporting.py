"""Errors in the unified (mode-2) C layer reach Python with their message.

gpu_set_error() records the reason on the domain struct and the Cython
wrappers raise RuntimeError with it (audit issue #341). This drives a real
error through that path: an inlet operator with a triangle index outside
the mesh is rejected by gpu_inlet_operator_init before anything is mapped.

It also guards against gpu_set_error echoing through itself, which is what
the first version of it did and which turned every error into a stack
overflow.
"""
import numpy as num
import pytest

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary


@pytest.fixture
def unified_domain():
    points, vertices, boundary = rectangular_cross(4, 4, len1=4.0, len2=4.0)
    d = Domain(points, vertices, boundary)
    d.set_compute_mode('legacy')
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    d.set_compute_mode('unified')
    return d


def test_bad_inlet_index_raises_with_the_c_layer_message(unified_domain):
    from anuga.shallow_water import sw_domain_gpu_ext as ext
    gpu_dom = unified_domain.gpu_interface.gpu_dom
    n = unified_domain.number_of_elements
    indices = num.array([0, n + 5], dtype=num.intc)          # second one is out of range
    areas = num.array([1.0, 1.0], dtype=num.float64)

    with pytest.raises(RuntimeError) as info:
        ext.init_inlet_operator(gpu_dom, indices, areas)

    assert 'out of range' in str(info.value)
    assert 'out of range' in ext.get_last_error(gpu_dom)


def test_a_good_init_leaves_no_error_and_frees_its_slot(unified_domain):
    from anuga.shallow_water import sw_domain_gpu_ext as ext
    gpu_dom = unified_domain.gpu_interface.gpu_dom
    indices = num.array([0, 1], dtype=num.intc)
    areas = num.array([1.0, 1.0], dtype=num.float64)

    op_id = ext.init_inlet_operator(gpu_dom, indices, areas)
    assert op_id >= 0
    ext.finalize_inlet_operator(gpu_dom, op_id)
