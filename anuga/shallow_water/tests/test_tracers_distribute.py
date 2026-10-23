"""Tracers ride the quantity partitioning to the sub-domains (issue #278).

The parallel end-to-end check needs MPI and lives in
anuga/parallel/tests/test_parallel_tracers.py. These test the partition
machinery directly, in one process, so the pieces are covered wherever the
suite runs.

The design under test: tracers are not Quantity objects (#276), so rather than
building a second partitioning path they travel as reserved entries in the
dicts that already cross the wire -- the concentrations in `quantities`, the
names and shared beta in `kwargs`.
"""

import numpy as num
import pytest

import anuga

pytest.importorskip('anuga.parallel.distribute_mesh',
                    reason='the parallel module is not importable here')

from anuga.parallel.distribute_mesh import (      # noqa: E402
    TRACER_QUANTITY_PREFIX,
    partition_mesh,
    pop_tracer_quantities,
    restore_tracers,
    tracer_partition_metadata,
    tracer_partition_quantities,
)


def _domain(n=6):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    return d


def _wedge(d):
    x = d.centroid_coordinates[:, 0]
    y = d.centroid_coordinates[:, 1]
    return 0.2 + 0.5 * x + 0.3 * y


def _with_tracers(n=6):
    d = _domain(n)
    d.add_tracer('salinity', beta=1.5)
    d.add_tracer('dye', beta=1.5)
    d.set_tracer('salinity', _wedge(d))
    d.set_tracer('dye', 0.25)
    return d


# --- what gets handed to the partitioner -----------------------------------

def test_a_domain_without_tracers_carries_nothing_extra():
    d = _domain()
    assert tracer_partition_quantities(d) == {}
    assert tracer_partition_metadata(d) is None


def test_the_metadata_records_the_names_in_order_and_the_shared_beta():
    d = _with_tracers()
    meta = tracer_partition_metadata(d)
    assert meta['names'] == ['salinity', 'dye'], 'slot order must be preserved'
    assert meta['beta'] == 1.5


def test_the_values_are_keyed_by_slot_not_by_name():
    """A name is free text; the slot is what has to survive."""
    d = _with_tracers()
    q = tracer_partition_quantities(d)
    assert sorted(q) == [TRACER_QUANTITY_PREFIX + '0',
                         TRACER_QUANTITY_PREFIX + '1']
    assert num.allclose(q[TRACER_QUANTITY_PREFIX + '0'], _wedge(d))
    assert num.allclose(q[TRACER_QUANTITY_PREFIX + '1'], 0.25)


def test_partition_mesh_includes_the_tracers():
    """They must ride the same reordering as stage and friction."""
    d = _with_tracers()
    before = _wedge(d).copy()

    _, _, quantities, _, epart = partition_mesh(d, 3)

    key = TRACER_QUANTITY_PREFIX + '0'
    assert key in quantities, 'tracers were left out of the partition'
    assert 'stage' in quantities, 'the ordinary quantities are still there'
    # Same permutation the quantities got.
    assert num.allclose(quantities[key][:, 0], before[epart])


# --- and what comes back out on the other side -----------------------------

def test_popping_returns_slot_order_and_removes_the_keys():
    d = _with_tracers()
    quantities = {'stage': 'S'}
    quantities.update(tracer_partition_quantities(d))

    values = pop_tracer_quantities(quantities)

    assert list(quantities) == ['stage'], 'reserved keys were left behind'
    assert len(values) == 2
    assert num.allclose(values[0], _wedge(d))
    assert num.allclose(values[1], 0.25)


def test_popping_a_dict_with_no_tracers_returns_none():
    assert pop_tracer_quantities({'stage': 'S'}) is None


def test_non_contiguous_slots_are_refused():
    quantities = {TRACER_QUANTITY_PREFIX + '0': num.zeros(3),
                  TRACER_QUANTITY_PREFIX + '2': num.zeros(3)}
    with pytest.raises(ValueError, match='slots'):
        pop_tracer_quantities(quantities)


def test_restore_rebuilds_names_order_beta_and_values():
    src = _with_tracers()
    meta = tracer_partition_metadata(src)
    values = pop_tracer_quantities(dict(tracer_partition_quantities(src)))

    target = _domain()                      # same mesh, no tracers yet
    restore_tracers(target, meta, values)

    assert target.get_tracer_names() == ['salinity', 'dye']
    assert target.beta_tracer == 1.5
    assert num.allclose(target.get_tracer('salinity'), _wedge(src))
    assert num.allclose(target.get_tracer('dye'), 0.25)


def test_restore_derives_the_conserved_variable_from_the_local_depth():
    """m = h*c must be rebuilt from the sub-domain's own h, not carried."""
    src = _with_tracers()
    meta = tracer_partition_metadata(src)
    values = pop_tracer_quantities(dict(tracer_partition_quantities(src)))

    target = _domain()
    target.set_quantity('stage', 2.5)       # deliberately a different depth
    restore_tracers(target, meta, values)

    h = (target.quantities['stage'].centroid_values
         - target.quantities['elevation'].centroid_values)
    assert num.allclose(target.tracer_conserved_values[0], h * _wedge(src))
    assert not num.allclose(target.tracer_conserved_values[0],
                            src.tracer_conserved_values[0]), \
        'm was copied across rather than recomputed from the local depth'


def test_restore_is_a_no_op_for_a_partition_with_no_tracers():
    """Partition files written before tracers existed must still load."""
    target = _domain()
    restore_tracers(target, None, None)
    assert target.number_of_tracers == 0


def test_half_a_partition_is_refused():
    d = _with_tracers()
    target = _domain()
    with pytest.raises(ValueError, match='cannot be reconstructed'):
        restore_tracers(target, tracer_partition_metadata(d), None)


def test_a_name_count_mismatch_is_refused():
    d = _with_tracers()
    meta = tracer_partition_metadata(d)
    values = pop_tracer_quantities(dict(tracer_partition_quantities(d)))
    with pytest.raises(ValueError, match='name'):
        restore_tracers(_domain(), meta, values[:1])


# --- the serial bypass ------------------------------------------------------

def test_distribute_on_one_process_leaves_the_domain_alone():
    from anuga.parallel.parallel_api import distribute
    from anuga import numprocs

    d = _with_tracers()
    if numprocs == 1:
        assert distribute(d) is d
        assert d.get_tracer_names() == ['salinity', 'dye']


def test_the_baseline_capture_does_not_use_a_collective():
    """add_tracer runs on rank 0 alone, before distribute() -- see #278.

    If the baseline went through the reducing get_tracer_mass, rank 0 would
    enter a collective the other ranks never reach, and the run would hang
    rather than fail. Pin the local, collective-free path.
    """
    d = _domain()
    d.add_tracer('c', initial_value=0.4)
    s = d.get_tracer_index('c')
    expected = float(num.sum(d.tracer_conserved_values[s] * d.areas))
    assert num.isclose(d._tracer_initial_mass[s], expected, rtol=1e-12)
    assert num.isclose(d._local_tracer_mass(s), expected, rtol=1e-12)
