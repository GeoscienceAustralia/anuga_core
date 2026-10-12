"""Generate a standalone ANUGA run script from a TOML scenario.

`anuga_toml_run --emit-script run.py cfg.toml` writes a plain Python script that
reproduces the runner's phases (mesh -> initial conditions -> forcing/structures
-> boundaries -> evolve) by parsing the same TOML through ``PrepareData`` and
calling the same ``setup_*`` helpers. It runs identically to ``anuga_toml_run``
because it calls the same code, and it is an ordinary editable script: users can
delete phases, splice in their own operators/boundaries, or replace a
``setup_*`` call with hand-written ANUGA code.

The script is a base for adaptation, not a faithful flattening of every TOML
field into explicit API calls.
"""

import os

_TEMPLATE = '''\
#!/usr/bin/env python
"""ANUGA run script generated from {config!r} by ``anuga_toml_run --emit-script``.

Equivalent to running::

    anuga_toml_run {config}

It parses the same TOML via ``PrepareData`` and drives the standard phases
(mesh -> initial conditions -> forcing/structures -> boundaries -> evolve).
Everything below is ordinary Python -- edit it freely: delete phases you do not
need, splice in your own operators/boundaries, or replace a ``setup_*`` call
with hand-written ANUGA code. Run it directly or under MPI::

    python {script}
    mpiexec -np 4 python {script}     # legacy (CPU) parallel

For ``compute_mode = "unified"`` with GPU offload, match the MPI rank count to
the number of GPUs.
"""

import os

import anuga
from anuga import myid, numprocs, finalize, barrier
from anuga.scenario.prepare_data import PrepareData
from anuga.scenario import (
    setup_mesh,
    setup_initial_conditions,
    setup_riverwalls,
    setup_rainfall,
    setup_inlets,
    setup_bridges,
    setup_culverts,
    setup_weirs,
    setup_pumping_stations,
    setup_erosion,
    setup_sediment,
    setup_boundary_conditions,
)

# Resolve TOML-relative paths from this script's directory (assumes the script
# lives alongside the config). Remove/adjust if you keep them elsewhere.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------------------
# Parse the TOML and pre-process inputs (mesh geometry, output directory,
# spatial data). ``project`` carries every parsed field.
# ---------------------------------------------------------------------------
project = PrepareData({config!r}, output_log='Simulation_logfile.log')

# ---------------------------------------------------------------------------
# Mesh (build + partition/distribute under MPI)
# ---------------------------------------------------------------------------
domain = setup_mesh.setup_mesh(project)

# Compute path: 'legacy' (CPU OpenMP) or 'unified' (shared CPU/GPU kernels).
domain.set_compute_mode(getattr(project, 'compute_mode', 'legacy'))
domain.set_omp_num_threads(project.omp_num_threads)

# ---------------------------------------------------------------------------
# Initial conditions (quantities set at centroids) + riverwalls
# ---------------------------------------------------------------------------
setup_initial_conditions.setup_initial_conditions(domain, project)
setup_riverwalls.setup_riverwalls(domain, project)

# ---------------------------------------------------------------------------
# Forcing & hydraulic structures. Each helper is a no-op when the TOML has no
# data for it. Add your own operators here, e.g.::
#     from anuga import Rate_operator
#     Rate_operator(domain, rate=0.001, polygon=my_polygon)
# ---------------------------------------------------------------------------
setup_rainfall.setup_rainfall(domain, project)
setup_inlets.setup_inlets(domain, project)
setup_bridges.setup_bridges(domain, project)
setup_culverts.setup_culverts(domain, project)
setup_weirs.setup_weirs(domain, project)
setup_pumping_stations.setup_pumping_stations(domain, project)
setup_erosion.setup_erosion(domain, project)

# Sediment transport (suspended load, optional bedload and bed evolution),
# from the [sediment] table. A no-op without one.
setup_sediment.setup_sediment(domain, project)

# ---------------------------------------------------------------------------
# Boundary conditions (tags -> boundary objects, from the TOML)
# ---------------------------------------------------------------------------
setup_boundary_conditions.setup_boundary_conditions(domain, project)

# ---------------------------------------------------------------------------
# Evolve. yieldstep/outputstep/finaltime come from the TOML; replace with
# literals if you prefer. Customise the per-step reporting freely.
# ---------------------------------------------------------------------------
barrier()
for t in domain.evolve(yieldstep=project.yieldstep,
                       outputstep=project.outputstep,
                       finaltime=project.finaltime):
    if myid == 0:
        print(domain.timestepping_statistics())

# ---------------------------------------------------------------------------
# Finalise: merge the per-process SWW files (parallel) and shut down MPI.
# ---------------------------------------------------------------------------
barrier()
if numprocs > 1:
    domain.sww_merge(delete_old=True)
finalize()
'''


def build_run_script(config_basename, script_path=None):
    """Return the text of a standalone ANUGA run script for *config_basename*.

    *script_path* (optional) is only used to fill in the ``python <script>``
    example in the generated docstring; it does not have to exist.
    """
    script = os.path.basename(script_path) if script_path else 'run.py'
    return _TEMPLATE.format(config=config_basename, script=script)


def write_run_script(config_basename, script_path):
    """Write the generated run script to *script_path* and return that path."""
    text = build_run_script(config_basename, script_path=script_path)
    with open(script_path, 'w', encoding='utf-8') as fh:
        fh.write(text)
    return script_path
