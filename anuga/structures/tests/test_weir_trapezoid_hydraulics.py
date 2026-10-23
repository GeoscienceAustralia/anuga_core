"""Hydraulics of the weir/orifice trapezoid culvert, against physical checks.

Fixed together, in the Python function and the shared C kernel (mode 1 for fully
local culverts, mode 2 on the GPU):

* several barrels: only the bottom width was multiplied by the number of
  barrels, so N barrels acted as one wide trapezoid. Areas, perimeters and the
  critical depth are now per barrel, times N;
* running full, the wetted perimeter left out the roof (outlet submerged and
  outlet flowing full), and the Python's part-full perimeter wrongly included it;
* the C submerged-inlet (orifice) area already held blockage and barrels and was
  multiplied by both again;
* the barrel-friction limit was skipped once the head difference reached the
  driving energy, so Q jumped (as for the Boyd box).
"""
import pytest

from anuga.structures.boyd_box_operator import boyd_box_function
from anuga.structures.weir_orifice_trapezoid_operator import weir_orifice_trapezoid_function

G = 9.8

# (driving energy, outlet water level above the inlet invert, outlet depth)
REGIMES = {
    'weir inlet, free outlet': (0.6, -1.0, 0.0),
    'orifice inlet, free outlet': (3.0, -1.0, 0.0),
    'outlet control, part full': (1.5, 1.0, 0.4),
    'outlet control, submerged': (2.0, 1.8, 1.2),
}


def trap(regime, barrels=1.0, blockage=0.0, z=0.0, width=1.0, depth=1.0,
         length=100.0, manning=0.02):
    E, outlet_level, tw = REGIMES[regime]
    return weir_orifice_trapezoid_function(
        width=width, depth=depth, blockage=blockage, barrels=barrels, z1=z, z2=z,
        flow_width=width, length=length, driving_energy=E,
        delta_total_energy=E - outlet_level, outlet_enquiry_depth=tw,
        sum_loss=1.5, manning=manning, g=G)


@pytest.mark.parametrize('regime', sorted(REGIMES))
@pytest.mark.parametrize('z', [0.0, 1.0])
@pytest.mark.parametrize('blockage', [0.0, 0.3])
def test_n_barrels_carry_n_times_one(regime, z, blockage):
    q1 = trap(regime, 1.0, blockage, z)[0]
    q3 = trap(regime, 3.0, blockage, z)[0]
    assert q1 > 0.0
    assert q3 == pytest.approx(3.0 * q1, rel=1e-9)


@pytest.mark.parametrize('tw', [1.2, 0.95])
def test_a_rectangular_barrel_running_full_matches_the_boyd_box(tw):
    """Same section, same energy equation: a long submerged (or full) culvert
    is friction-controlled, and its wetted perimeter must include the roof."""
    E, length, n = 2.0, 300.0, 0.024
    dE = 0.4
    q_trap = weir_orifice_trapezoid_function(
        width=1.2, depth=0.9, blockage=0.0, barrels=1.0, z1=0.0, z2=0.0,
        flow_width=1.2, length=length, driving_energy=E, delta_total_energy=dE,
        outlet_enquiry_depth=tw, sum_loss=1.5, manning=n, g=G)[0]
    q_box = boyd_box_function(width=1.2, depth=0.9, blockage=0.0, barrels=1.0,
                              flow_width=1.2, length=length, driving_energy=E,
                              delta_total_energy=dE, outlet_enquiry_depth=tw,
                              sum_loss=1.5, manning=n)[0]
    assert q_trap == pytest.approx(q_box, rel=1e-9)


@pytest.mark.parametrize('length,manning', [(50.0, 0.013), (200.0, 0.024), (500.0, 0.024)])
def test_no_jump_where_the_head_difference_reaches_the_driving_energy(length, manning):
    E = 1.5

    def q(dE):
        return weir_orifice_trapezoid_function(
            width=1.0, depth=1.0, blockage=0.0, barrels=1.0, z1=0.0, z2=0.0,
            flow_width=1.0, length=length, driving_energy=E, delta_total_energy=dE,
            outlet_enquiry_depth=0.0, sum_loss=1.5, manning=manning, g=G)[0]

    below, above = q(E - 0.01), q(E + 0.01)
    assert below > 0.0
    assert above == pytest.approx(below, rel=0.03)


def _c_kernel_discharge(regime, barrels, blockage, z, length=100.0, manning=0.02):
    """Raw discharge from the shared C kernel (huge inlets, tiny timestep, so no
    volume limit applies), velocity head off, in the same state as `trap`."""
    ext = pytest.importorskip('anuga.shallow_water.sw_domain_gpu_ext')
    E, outlet_level, tw = REGIMES[regime]
    z_out = outlet_level - tw
    dt, big = 1.0e-6, 1.0e12
    res = ext.culvert_apply_one_host_py(
        2, G, 1.0, 1.0, 0.0, z, z,                   # weir trapezoid; width, height, diameter, z1, z2
        length, manning, 1.5, blockage, barrels,
        0, 1, 1, 1, 1.0e9, 0.0,
        1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0, 0,
        0.0, 0.0, dt,
        E, 0.0, 0.0, 0.0, E, E, 0.0, 0.0, big,
        outlet_level, 0.0, 0.0, z_out, tw, tw, 0.0, 0.0, big)
    return res[10]


@pytest.mark.parametrize('regime', sorted(REGIMES))
@pytest.mark.parametrize('barrels,blockage,z', [(1.0, 0.0, 0.0), (3.0, 0.0, 1.0),
                                                (2.0, 0.5, 0.0), (3.0, 0.3, 0.5)])
def test_c_kernel_matches_python(regime, barrels, blockage, z):
    q_py = trap(regime, barrels, blockage, z)[0]
    q_c = _c_kernel_discharge(regime, barrels, blockage, z)
    assert q_c == pytest.approx(q_py, rel=1e-6)
