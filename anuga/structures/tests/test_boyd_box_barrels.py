"""Boyd box culverts with several barrels.

N identical barrels side by side must carry exactly N times one barrel: each
barrel has its own walls, so the wetted perimeter scales with N exactly as the
flow area does and the hydraulic radius is that of one barrel. The perimeter
used to be taken as N*width + 2*depth -- N boxes treated as one wide box with
two walls -- which overstated the hydraulic radius and the outlet-controlled
discharge (3 barrels carried up to 3.6x one barrel on a 200 m culvert).

Checked on the Python function and, through a domain, on the shared C kernel
that mode 1 uses for fully local culverts (and mode 2 on the GPU).
"""
import numpy as np
import pytest

import anuga
from anuga.structures.boyd_box_operator import boyd_box_function

CASES = {
    'inlet control, free outlet': dict(driving_energy=0.6, delta_total_energy=0.6,
                                       outlet_enquiry_depth=0.0, length=20.0),
    'outlet control, part full': dict(driving_energy=1.5, delta_total_energy=0.3,
                                      outlet_enquiry_depth=0.5, length=20.0),
    'outlet control, submerged': dict(driving_energy=2.0, delta_total_energy=0.2,
                                      outlet_enquiry_depth=1.2, length=20.0),
    'long culvert, part full': dict(driving_energy=1.5, delta_total_energy=0.3,
                                    outlet_enquiry_depth=0.5, length=200.0),
    'long culvert, submerged': dict(driving_energy=2.0, delta_total_energy=0.2,
                                    outlet_enquiry_depth=1.2, length=200.0),
}


def _box(barrels, case, blockage=0.0):
    return boyd_box_function(width=1.2, depth=0.9, blockage=blockage, barrels=barrels,
                             flow_width=1.2, sum_loss=1.5, manning=0.013,
                             **CASES[case])


@pytest.mark.parametrize('case', sorted(CASES))
@pytest.mark.parametrize('blockage', [0.0, 0.3])
def test_n_barrels_carry_n_times_one(case, blockage):
    q1, v1, d1, a1, p1 = _box(1, case, blockage)
    q3, v3, d3, a3, p3 = _box(3, case, blockage)
    assert q1 > 0.0
    assert q3 == pytest.approx(3.0 * q1, rel=1e-12)
    # same velocity and depth in every barrel as in the single one
    assert v3 == pytest.approx(v1, rel=1e-5)   # reported velocity carries a small regularisation
    assert d3 == pytest.approx(d1, rel=1e-12)


def _domain_with_culvert(barrels):
    d = anuga.rectangular_cross_domain(50, 10, len1=100.0, len2=20.0)
    d.set_flow_algorithm('DE0')
    d.set_store(False)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('friction', 0.01)
    d.set_quantity('stage', lambda x, y: np.where(x < 50.0, 2.0, 1.2))
    d.set_boundary({t: anuga.Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.set_evolve_max_timestep(0.01)
    op = anuga.Boyd_box_operator(d, end_points=[[20.0, 10.0], [80.0, 10.0]],
                                 width=1.2, height=0.9, barrels=barrels,
                                 losses=1.5, manning=0.013,
                                 use_velocity_head=False, smoothing_timescale=0.0,
                                 verbose=False)
    return d, op


def test_c_kernel_matches_python_for_several_barrels():
    """A long submerged culvert (outlet, i.e. friction, control: the regime the
    old perimeter got wrong) with 3 barrels. The shared C kernel, which mode 1
    uses for fully local culverts and mode 2 on the GPU, must give exactly the
    Python operator's discharge; and the hydraulic discharge must scale 3:1."""
    import anuga.structures.structure_operator as so
    from anuga.structures.structure_operator import _can_use_c_culvert

    def one_update(n, use_c):
        saved = so._c_culvert_funcs
        if not use_c:
            so._c_culvert_funcs = (None, None)
        try:
            d, op = _domain_with_culvert(n)
            assert _can_use_c_culvert(op) == use_c
            d.timestep = 0.01
            d.yieldstep = 0.01
            op()                          # one update on the initial state
            return op.discharge
        finally:
            so._c_culvert_funcs = saved

    q_c, q_py = one_update(3, True), one_update(3, False)
    assert q_c > 0.0
    assert q_c == pytest.approx(q_py, rel=1e-12)

    q = {}
    for n in (1, 3):
        d, op = _domain_with_culvert(n)
        q[n] = op.discharge_routine()[0]
    assert q[3] == pytest.approx(3.0 * q[1], rel=1e-12)
