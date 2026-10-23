"""Runner for the negative-cells MPI deadlock regression test.

Reproduces the exact condition that used to deadlock update_conserved_quantities()
in parallel: the "possible loss of conservation" warning computed the total
water volume with the collective get_water_volume() (an MPI_Allreduce), but the
call sat inside `if num_negative_ids > 0` — a per-rank condition. When only some
ranks clamp negative-depth cells (routine in wetting/drying), only those ranks
entered the allreduce while the others advanced to the ghost exchange, hanging
the run forever.

Here we force ONLY rank 0 to have negative-depth cells, leave the other ranks
clean, then call update_conserved_quantities(). Before the fix this hangs; after
it (the warning uses a local volume, no collective) every rank returns and the
process exits 0. Run under mpiexec by test_parallel_negative_cells_deadlock.py;
the wrapper's subprocess timeout turns a regression back into a hang -> failure.
"""

import anuga
from anuga.parallel import distribute, myid, numprocs, finalize, barrier

# Abort the whole job (not hang) if any rank raises, so a genuine assertion
# failure surfaces instead of deadlocking the surviving ranks.
import sys
from anuga.utilities.parallel_abstraction import global_except_hook
sys.excepthook = global_except_hook


def main():
    # Small flat domain: 1 m of water over a 0 m bed.
    if myid == 0:
        domain = anuga.rectangular_cross_domain(10, 10)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_boundary({'left': anuga.Reflective_boundary(domain),
                             'right': anuga.Reflective_boundary(domain),
                             'top': anuga.Reflective_boundary(domain),
                             'bottom': anuga.Reflective_boundary(domain)})
    else:
        domain = None

    domain = distribute(domain)

    # Legacy (mode 1) CPU kernel regardless of ANUGA_DEFAULT_COMPUTE_MODE, and no
    # flux change this step so we isolate the negative-cell clamping path.
    domain.set_compute_mode('legacy')
    domain.timestep = 0.0

    # Force negative depth on a few of RANK 0's local cells only. The other
    # rank(s) stay wet, so num_negative_ids > 0 on rank 0 and == 0 elsewhere:
    # the divergent condition that used to deadlock the collective.
    if myid == 0:
        stage_c = domain.quantities['stage'].centroid_values
        bed_c = domain.quantities['elevation'].centroid_values
        stage_c[:3] = bed_c[:3] - 0.5

    barrier()
    # THE call under test. Before the fix, rank 0 blocks here in the volume
    # allreduce while the others move on -> hang.
    domain.update_conserved_quantities()
    barrier()

    if myid == 0:
        print('NEGATIVE_CELLS_DEADLOCK_OK')

    finalize()


if __name__ == '__main__':
    main()
