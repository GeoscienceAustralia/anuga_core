"""The barrel-friction limit of Boyd box culverts.

Outlet control computes the flow the barrel's friction allows, Q_outlet_tailwater,
from the energy equation over the culvert length. Boyd pipe applies Q = min(Q,
Q_outlet_tailwater) in every regime; box used to apply it only while
the head difference delta_total_energy stayed below the driving energy. Once it
reached the driving energy (a free or drawn-down outlet) the friction limit was
skipped and Q jumped (1.9x on a 200 m, n = 0.024 culvert, 2.9x at 500 m),
over-predicting long or rough culverts. (The weir/orifice trapezoid has the
same gap; it is fixed together with its other outlet-control problems.)

Checked on the Python functions and on the shared C kernel that mode 1 uses for
fully local culverts and mode 2 on the GPU.
"""
import pytest

from anuga.structures.boyd_box_operator import boyd_box_function

G = 9.8
E = 1.5                       # driving energy [m]
LONG = dict(length=200.0, manning=0.024, sum_loss=1.5)


def box(delta_total_energy, length=200.0, manning=0.024):
    return boyd_box_function(width=1.0, depth=1.0, blockage=0.0, barrels=1.0,
                             flow_width=1.0, length=length, driving_energy=E,
                             delta_total_energy=delta_total_energy,
                             outlet_enquiry_depth=0.0, sum_loss=1.5,
                             manning=manning)[0]


@pytest.mark.parametrize('length,manning', [(50.0, 0.013), (200.0, 0.024), (500.0, 0.024)])
def test_no_jump_where_the_head_difference_reaches_the_driving_energy(length, manning):
    below = box(E - 0.01, length, manning)
    above = box(E + 0.01, length, manning)
    assert below > 0.0
    assert above == pytest.approx(below, rel=0.03)


def _c_kernel_discharge(ctype, delta_total_energy):
    """Raw discharge from the shared C kernel: huge inlets and a tiny timestep,
    so no volume limit touches it. Water level E above the inlet invert, the
    outlet dry with its invert E - delta_total_energy."""
    ext = pytest.importorskip('anuga.shallow_water.sw_domain_gpu_ext')
    dt, big = 1.0e-6, 1.0e12
    z_out = E - delta_total_energy
    res = ext.culvert_apply_one_host_py(
        ctype, G, 1.0, 1.0, 0.0, 0.0, 0.0,           # width, height, diameter, z1, z2
        LONG['length'], LONG['manning'], LONG['sum_loss'], 0.0, 1.0,
        0, 1, 1, 1, 1.0e9, 0.0,                      # velocity head off, jet, ..., max_velocity, smoothing
        1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0, 0,         # outward vectors, inverts
        0.0, 0.0, dt,
        E, 0.0, 0.0, 0.0, E, E, 0.0, 0.0, big,       # inlet 0: stage, mom, z, averages, area
        z_out, 0.0, 0.0, z_out, 0.0, 0.0, 0.0, 0.0, big)
    return res[10]


@pytest.mark.parametrize('dE', [E + 0.01, E + 0.1, E + 1.0])
def test_c_kernel_matches_python_above_the_driving_energy(dE):
    assert _c_kernel_discharge(0, dE) == pytest.approx(box(dE), rel=1e-9)
