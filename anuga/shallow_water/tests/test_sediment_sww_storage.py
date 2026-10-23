"""Elevation must be stored per timestep once the bed can move.

The bed-exchange and bedload kernels write `bed_centroid_values` /
`bed_edge_values` in place, and those ARE the `elevation` Quantity's arrays.
So an evolving bed is time-varying -- but `quantities_to_be_stored` defaults to
`'elevation': 1`, meaning write it once.  An sww from such a run then records
the initial bed and omits every change, which looks complete and is wrong.
"""

import pytest

import anuga


def _domain(n=6):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    return d


def test_plain_domain_still_stores_elevation_statically():
    # No sediment: nothing moves the bed, so the default must be untouched.
    d = _domain()
    assert d.quantities_to_be_stored['elevation'] == 1


def test_adding_a_sediment_class_upgrades_elevation_to_dynamic():
    d = _domain()
    assert d.quantities_to_be_stored['elevation'] == 1
    d.add_sediment_fraction(name='sand', diameter=200e-6)
    assert d.sediment_bed_evolution is True
    assert d.quantities_to_be_stored['elevation'] == 2, \
        'an evolving bed must be stored per timestep, not once'


def test_fixed_bed_leaves_elevation_static():
    # Phase 3 / RDy26 v1.0 configuration: the bed does not move, so storing it
    # once is correct and we must not inflate the file.
    d = _domain()
    d.set_sediment_parameters(bed_evolution=False)
    d.add_sediment_fraction(name='sand', diameter=200e-6)
    assert d.quantities_to_be_stored['elevation'] == 1


def test_enabling_bed_evolution_afterwards_upgrades_it():
    d = _domain()
    d.set_sediment_parameters(bed_evolution=False)
    d.add_sediment_fraction(name='sand', diameter=200e-6)
    assert d.quantities_to_be_stored['elevation'] == 1
    d.set_sediment_parameters(bed_evolution=True)
    assert d.quantities_to_be_stored['elevation'] == 2


@pytest.mark.parametrize('flag', [3, 4])
def test_deliberate_flags_are_not_rewritten(flag):
    # 3 = centroid-only dynamic, 4 = overwritten each yieldstep. Both are
    # explicit choices; upgrading them to 2 would override the user.
    d = _domain()
    d.quantities_to_be_stored['elevation'] = flag
    d.add_sediment_fraction(name='sand', diameter=200e-6)
    assert d.quantities_to_be_stored['elevation'] == flag


def test_an_explicit_2_is_left_alone():
    d = _domain()
    d.quantities_to_be_stored['elevation'] = 2
    d.add_sediment_fraction(name='sand', diameter=200e-6)
    assert d.quantities_to_be_stored['elevation'] == 2


def test_the_evolving_bed_actually_reaches_the_sww(tmp_path):
    """End to end: the stored elevation must change over time, not repeat.

    This is the test that would have caught the original defect -- the flag
    itself is an implementation detail, whereas this asserts the property a
    modeller depends on.
    """
    netCDF4 = pytest.importorskip('netCDF4')

    d = _domain(8)
    d.set_quantity('elevation', lambda x, y: -0.05 * x)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    d.add_sediment_fraction(name='sand', diameter=200e-6, initial_concentration=0.02)
    d.set_bedload()   # default: wong_parker_eq24

    d.set_name('sed_bed_storage')
    d.set_datadir(str(tmp_path))
    for _ in d.evolve(yieldstep=0.5, finaltime=1.5):
        pass

    with netCDF4.Dataset(tmp_path / 'sed_bed_storage.sww') as ds:
        elev = ds.variables['elevation'][:]

    # Flag 1 stores a single (vertices,) row; flag 2 stores (time, vertices).
    assert elev.ndim == 2, \
        'elevation was stored statically, so bed evolution is invisible'
    assert elev.shape[0] > 1
