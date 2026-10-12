#!/usr/bin/env python3
"""
anuga_toml_run — run an ANUGA scenario from a TOML configuration file.

Usage (serial):
    anuga_toml_run  path/to/scenario.toml

Usage (parallel, N processes):
    mpirun -np N anuga_toml_run  path/to/scenario.toml

The working directory is changed to the directory that contains the TOML
file before the simulation starts, so all relative paths inside the TOML
are resolved relative to that directory.

A 'user_functions.py' module may be placed alongside the TOML file to
provide custom callbacks (print_velocity_statistics, print_operator_inputs).
If the file is absent those hooks are silently skipped.
"""

import argparse
import os
import sys
import time

parser = argparse.ArgumentParser(
    description='Run an ANUGA scenario defined by a TOML configuration file.')
parser.add_argument(
    'config',
    metavar='CONFIG.toml',
    help='Path to the TOML scenario configuration file.')
parser.add_argument(
    '-v', '--verbose', action='store_true',
    help='Echo the detailed library output (mesh generation, riverwall setup, '
         'raster writing, ...) to the terminal. By default that goes only to '
         'the run log file and the terminal shows just the phase banners and '
         'summary. The full log is always written either way.')
parser.add_argument(
    '-n', '--dry-run', action='store_true',
    help='Do not run the simulation. Preview the scenario in the format given '
         'by --format (default: open an HTML summary in the browser).')
parser.add_argument(
    '--format', choices=('text', 'html', 'browser'), default='browser',
    help='Dry-run output. "browser" (default): write the HTML summary and open '
         'it. "html": write the HTML summary only. "text": print the '
         'highlighted, folded config to the terminal (no browser needed).')
parser.add_argument(
    '--summary-output', metavar='FILE.html', default=None,
    help='For --format html/browser: where to write the summary '
         '(default: <config>_summary.html next to the TOML).')
# --format text options
parser.add_argument(
    '--full', action='store_true',
    help='--format text: do not collapse repeated blocks; show every one.')
parser.add_argument(
    '--no-color', action='store_true',
    help='--format text: disable syntax highlighting.')
parser.add_argument(
    '--no-pager', action='store_true',
    help='--format text: write straight to stdout instead of a pager.')
parser.add_argument(
    '--threshold', type=int, default=6, metavar='N',
    help='--format text: collapse a run only when it has more than N identical '
         'blocks (default: 6).')
parser.add_argument(
    '--emit-script', metavar='FILE',
    help='Write a standalone, editable ANUGA Python script equivalent to this '
         'run (mesh -> initial conditions -> forcing -> boundaries -> evolve) '
         'and exit without running. Use it as a base for customisation.')
args = parser.parse_args()

config_path = os.path.abspath(args.config)
if not os.path.exists(config_path):
    sys.exit(f'ERROR: config file not found: {config_path}')

# Change to the scenario directory so all relative paths in the TOML resolve
# correctly.  Do this before any anuga imports that may write files.
scenario_dir = os.path.dirname(config_path)
os.chdir(scenario_dir)
config_basename = os.path.basename(config_path)

# ---------------------------------------------------------------------------
# Emit a standalone run script and exit (no mesh build / evolve). Done before
# the heavy ANUGA imports so it stays fast.
# ---------------------------------------------------------------------------
if args.emit_script:
    from anuga.scenario.emit_script import write_run_script
    out = os.path.abspath(args.emit_script)
    write_run_script(config_basename, out)
    print(f'Run script written to: {out}')
    sys.exit(0)

# ---------------------------------------------------------------------------
# Dry run: preview the scenario without building a mesh or evolving.  Done
# before the heavy ANUGA/mesh imports so it stays fast.
#   text    -> highlighted, folded config to the terminal
#   html    -> write the HTML summary file
#   browser -> write the HTML summary and open it
# ---------------------------------------------------------------------------
if args.dry_run:
    if args.format == 'text':
        from anuga.utilities.toml_view import render, page
        with open(config_path, encoding='utf-8') as _fh:
            _text = _fh.read()
        _color = (not args.no_color) and (sys.stdout.isatty() or not args.no_pager)
        _out = render(_text, collapse=not args.full, color=_color,
                      threshold=args.threshold)
        if args.no_pager:
            sys.stdout.write(_out)
        else:
            page(_out)
    else:
        from anuga.scenario.scenario_summary import write_scenario_summary
        out = write_scenario_summary(
            config_path, output_html=args.summary_output,
            base_dir=scenario_dir, open_browser=(args.format == 'browser'))
        print(f'Scenario summary written to: {out}')
    sys.exit(0)

# ---------------------------------------------------------------------------
# ANUGA imports (after chdir so parallel init finds the right cwd)
# ---------------------------------------------------------------------------

import anuga
from anuga.parallel import myid, numprocs, finalize, barrier
from anuga.operators.collect_max_quantities_operator import \
    Collect_max_quantities_operator

from anuga.scenario import (
    setup_boundary_conditions,
    setup_rainfall,
    setup_inlets,
    setup_bridges,
    setup_culverts,
    setup_weirs,
    setup_pumping_stations,
    setup_mesh,
    setup_initial_conditions,
    setup_riverwalls,
    setup_erosion,
    setup_sediment,
    raster_outputs,
)
from anuga.scenario.prepare_data import PrepareData

# ---------------------------------------------------------------------------
# Optional user_functions module (lives alongside the TOML file)
# ---------------------------------------------------------------------------

sys.path.insert(0, scenario_dir)
try:
    import user_functions
    _have_user_functions = True
except ImportError:
    _have_user_functions = False

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

t0 = time.time()

# Where banners/summary go. Default is the real terminal; PrepareData installs a
# saved terminal handle in _terminal below once it knows the verbosity, so the
# banners still show even when the run's detailed stdout is redirected to the
# log file (quiet mode).
_terminal = sys.__stdout__


def progress(msg):
    """Write to the terminal only (regardless of stdout redirection)."""
    _terminal.write(msg + '\n')
    _terminal.flush()


def emit(msg=''):
    """Runner output that belongs in both the log and the terminal (once).

    print() goes to the log (fd 1 in quiet mode; also the terminal via the tee
    in --verbose); in quiet mode we additionally echo to the saved terminal.
    """
    print(msg)
    if not args.verbose:
        progress(msg)


_SECTION_BAR = '=' * 64


def section(title, n=None, total=None):
    """Print a consistent phase banner (log + terminal, rank 0 only)."""
    if myid != 0:
        return
    tag = f'[{n}/{total}]  ' if n is not None else ''
    emit('')
    emit(_SECTION_BAR)
    emit(f'  {tag}{title}')
    emit(_SECTION_BAR)


def compute_water_balance(domain):
    """Return (v0, fs, bf, vol, imbalance) for the mass-balance identity
    ``V = V0 + BF + FS``, or None if unavailable (non-DE).

    Collective: the volume getters reduce across ranks, so this must be called
    on every rank (print on rank 0 only). v0=initial volume, fs=fractional-step
    inflow (rainfall/inlets/operators), bf=net boundary flux, vol=current volume.
    """
    stats = domain.report_water_volume_statistics(verbose=False, returnStats=True)
    if stats is None:
        return None
    vol, bf, fs = stats
    # Prefer the volume captured before evolve (seeded on all ranks); fall back
    # to volume_history[0] only if that seed is somehow absent.
    v0 = getattr(domain, '_toml_initial_volume', None)
    if v0 is None:
        v0 = domain.volume_history[0] if getattr(domain, 'volume_history', None) else 0.0
    return (v0, fs, bf, vol, vol - v0 - bf - fs)


def quantity_summary(name):
    """One-line description of an initial-condition quantity's data list.

    Each entry is [polygon, value]; value is a constant or a file path. Returns
    None if the quantity was not configured.
    """
    entries = getattr(project, f'{name}_data', None) or []
    if not entries:
        return None
    if len(entries) == 1:
        val = entries[0][1]
        return os.path.basename(val) if isinstance(val, str) else ('%g' % val)
    # Multiple entries: a spatial source (file) dominates; otherwise it's a set
    # of constant-value zones, so report the count and numeric range.
    files = [e[1] for e in entries if isinstance(e[1], str)]
    if files:
        return '%s (+%d more)' % (os.path.basename(files[0]), len(entries) - 1)
    nums = [e[1] for e in entries if isinstance(e[1], (int, float))]
    return '%d zones (%g-%g)' % (len(entries), min(nums), max(nums))


# ---------------------------------------------------------------------------
# Load configuration
# ---------------------------------------------------------------------------

project = PrepareData(config_basename, output_log='Simulation_logfile.log',
                      echo_terminal=args.verbose)
# In quiet mode PrepareData redirected fd 1 to the log file; use the terminal
# handle it saved so banners/summary still reach the console.
_terminal = getattr(project, 'terminal', sys.__stdout__)

# ---------------------------------------------------------------------------
# [1/6] Mesh
# ---------------------------------------------------------------------------

section('MESH', 1, 6)
domain = setup_mesh.setup_mesh(project)
# Phase timings recorded by setup_mesh (build vs partition/distribute).
mesh_build_time = getattr(domain, '_mesh_build_time', None)
mesh_distribute_time = getattr(domain, '_mesh_distribute_time', None)

# Propagate the scenario's coordinate reference system to the domain so it is
# written into the SWW file (zone / hemisphere / EPSG). project.proj4string is
# derived from projection_information (UTM zone int, "EPSG:<code>", or a proj4
# string), so pyproj gives a single EPSG code covering all three forms.
epsg = None
try:
    from pyproj import CRS
    epsg = CRS.from_proj4(project.proj4string).to_epsg()
    if epsg is not None:
        domain.set_epsg(epsg)
    else:
        emit('   (could not resolve an EPSG code from projection %r)'
             % project.proj4string)
except Exception as _e:
    emit('   (could not set domain CRS: %s)' % _e)

# Global mesh size. After distribute, len(domain)/get_extent() are this rank's
# partition, so reduce across ranks: total non-ghost triangles (sum) and the
# global extent (min/max). This is collective — run on every rank.
_ext = domain.get_extent()  # (xmin, xmax, ymin, ymax) for this rank
if numprocs > 1:
    try:
        from mpi4py import MPI
        _c = MPI.COMM_WORLD
        _ntri = _c.allreduce(int(domain.number_of_full_triangles), op=MPI.SUM)
        _xmin = _c.allreduce(_ext[0], op=MPI.MIN)
        _xmax = _c.allreduce(_ext[1], op=MPI.MAX)
        _ymin = _c.allreduce(_ext[2], op=MPI.MIN)
        _ymax = _c.allreduce(_ext[3], op=MPI.MAX)
    except Exception:
        _ntri = None
else:
    _ntri = len(domain)
    _xmin, _xmax, _ymin, _ymax = _ext

if myid == 0:
    if _ntri is not None:
        emit('   {:,} triangles   extent {:.0f} x {:.0f} m{}'.format(
            _ntri, _xmax - _xmin, _ymax - _ymin,
            '   EPSG:%d' % epsg if epsg else ''))
    elif epsg:
        emit('   EPSG:%d' % epsg)
    emit('   background resolution %g m^2   |   interior regions %d   |   '
         'riverwalls %d'
         % (getattr(project, 'default_res', float('nan')),
            len(getattr(project, 'interior_regions_data', []) or []),
            len(getattr(project, 'riverwall_csv_files', []) or [])))

# ---------------------------------------------------------------------------
# [2/6] Initial conditions (quantities + riverwalls, added after distribute)
# ---------------------------------------------------------------------------

section('INITIAL CONDITIONS', 2, 6)
setup_initial_conditions.setup_initial_conditions(domain, project)
setup_riverwalls.setup_riverwalls(domain, project)
if myid == 0:
    _parts = ['%s: %s' % (q, quantity_summary(q))
              for q in ('elevation', 'friction', 'stage')
              if quantity_summary(q) is not None]
    if _parts:
        emit('   ' + '   '.join(_parts))

# ---------------------------------------------------------------------------
# [3/6] Forcing & structures (rainfall, inlets, bridges, pumps, erosion)
# ---------------------------------------------------------------------------

section('FORCING & STRUCTURES', 3, 6)
setup_rainfall.setup_rainfall(domain, project)
setup_inlets.setup_inlets(domain, project)
setup_bridges.setup_bridges(domain, project)
setup_culverts.setup_culverts(domain, project)
setup_weirs.setup_weirs(domain, project)
setup_pumping_stations.setup_pumping_stations(domain, project)
# Erosion operators change elevation during the run; added after the other
# forcing terms and before boundary conditions.
setup_erosion.setup_erosion(domain, project)
# Sediment transport from the [sediment] table; a no-op without one.
setup_sediment.setup_sediment(domain, project)
if myid == 0:
    _counts = [
        ('rainfall inputs', len(getattr(project, 'rain_data', []) or [])),
        ('inlets', len(getattr(project, 'inlet_data', []) or [])),
        ('culverts', len(getattr(project, 'culvert_data', []) or [])),
        ('weirs', len(getattr(project, 'weir_data', []) or [])),
        ('bridges', len(getattr(project, 'bridge_data', []) or [])),
        ('pumping stations', len(getattr(project, 'pumping_station_data', []) or [])),
        ('erosion operators', len(getattr(project, 'erosion_data', []) or [])),
    ]
    _shown = ['%s %d' % (n, c) for n, c in _counts if c]
    emit('   ' + ('   '.join(_shown) if _shown else '(none)'))

# ---------------------------------------------------------------------------
# [4/6] Boundary conditions
# ---------------------------------------------------------------------------

section('BOUNDARY CONDITIONS', 4, 6)
setup_boundary_conditions.setup_boundary_conditions(domain, project)
if myid == 0:
    _bd = getattr(project, 'boundary_data', []) or []
    _pairs = '   '.join('%s=%s' % (r[0], r[1]) for r in _bd)
    emit('   %d boundaries:  %s' % (len(_bd), _pairs))

# ---------------------------------------------------------------------------
# Track maximum quantities
# ---------------------------------------------------------------------------

max_quantities = Collect_max_quantities_operator(
    domain,
    update_frequency=project.max_quantity_update_frequency,
    collection_start_time=project.max_quantity_collection_start_time,
    velocity_zero_height=1.0e-03)

# ---------------------------------------------------------------------------
# [5/6] Evolve
# ---------------------------------------------------------------------------

if hasattr(project, 'compute_mode'):
    domain.set_compute_mode(project.compute_mode)   # 'legacy' or 'unified'
else:
    domain.set_multiprocessor_mode(project.multiprocessor_mode)  # Excel back-compat
domain.set_omp_num_threads(project.omp_num_threads)

# Guard the classic MPI-GPU footgun: mode-2 (unified) with GPU offload hangs /
# garbles unless the MPI rank count matches the number of GPUs (ranks map
# round-robin device_id = rank % num_devices, so oversubscribed ranks share a
# device and the omp-target runtime deadlocks). Every rank evaluates the same
# condition, so aborting here is collective — no deadlock, and a clear message
# instead of a silent hang. See claude/KNOWN_ISSUES.md.
if (numprocs > 1
        and getattr(domain, 'get_compute_mode', lambda: 'legacy')() == 'unified'
        and anuga.gpu_offload_enabled()):
    try:
        from anuga.shallow_water.sw_domain_gpu_ext import get_num_gpu_devices
        _ndev = int(get_num_gpu_devices())
    except Exception:
        _ndev = 0
    if _ndev <= 0 or numprocs != _ndev:
        if myid == 0:
            _terminal.write(
                '\nERROR: compute_mode = "unified" with GPU offload needs the MPI '
                'rank count to\n       match the number of GPUs. This run has '
                '%d rank(s) and %d GPU(s), which\n       would hang. Re-run with '
                '-np %d, or set compute_mode = "legacy" for a\n       multi-rank '
                'CPU run.\n' % (numprocs, _ndev, _ndev if _ndev > 0 else 1))
            _terminal.flush()
        finalize()
        sys.exit(1)

section('EVOLVE', 5, 6)
if myid == 0:
    _mode = domain.get_compute_mode() if hasattr(domain, 'get_compute_mode') else 'legacy'
    if _mode == 'unified':
        _mode_label = ('unified, GPU offload' if anuga.gpu_offload_enabled()
                       else 'unified, CPU multicore')
    else:
        _mode_label = 'legacy (CPU OpenMP)'
    _ostep = project.outputstep if project.outputstep is not None else project.yieldstep
    emit('   algorithm %s   compute %s   OMP threads %s'
         % (project.flow_algorithm, _mode_label, domain.omp_num_threads))
    emit('   yieldstep %g s   outputstep %g s   finaltime %g s'
         % (project.yieldstep, _ostep, project.finaltime))

# V0 for the end-of-run water balance is captured at the first yieldstep (t=0)
# inside the loop, where the clamped 'height' is valid. Computing it here, before
# evolve, would take get_water_volume()'s pre-evolve path (unclamped
# stage-elevation) and count negative "depth" in dry topographic cells — a large
# bogus V0. Seed the slot so compute_water_balance has a defined initial, and so
# the first-yield capture below triggers regardless of
# report_mass_conservation_statistics.
domain._toml_initial_volume = None

barrier()
evolve_start = time.time()
for t in domain.evolve(yieldstep=project.yieldstep,
                       finaltime=project.finaltime,
                       outputstep=project.outputstep):
    # First yield is t=0 (initial state): record V0 on every rank
    # (get_water_volume is collective) so the balance is correct even when
    # per-yieldstep reporting is off.
    if domain._toml_initial_volume is None:
        domain._toml_initial_volume = domain.get_water_volume()
    if myid == 0:
        _stats = domain.timestepping_statistics()
        # Blank line before each yieldstep, indented to the section's 3 spaces.
        _stats = '\n' + '\n'.join('   ' + _ln for _ln in _stats.splitlines())
        print(_stats)                       # -> log (and terminal via tee if -v)
        if not args.verbose:
            progress(_stats)                # keep the quiet terminal informed

    if project.report_mass_conservation_statistics:
        _wb = compute_water_balance(domain)   # collective (all ranks)
        if myid == 0 and _wb is not None:
            _v0, _fs, _bf, _vol, _resid = _wb
            _den = max(abs(_vol), abs(_v0) + abs(_fs) + abs(_bf), 1.0)
            _bal = ('   └─ balance:  V=%.2f  FS=%.2f  BF=%.2f  '
                    'imbalance=%.2e (%.1e rel)'
                    % (_vol, _fs, _bf, _resid, _resid / _den))
            print(_bal)
            if not args.verbose:
                progress(_bal)

    if project.report_peak_velocity_statistics and _have_user_functions:
        user_functions.print_velocity_statistics(domain, max_quantities)

    if project.report_smallest_edge_timestep_statistics:
        domain.report_cells_with_small_local_timestep()

    if project.report_operator_statistics and _have_user_functions:
        user_functions.print_operator_inputs(domain)

barrier()
evolve_time = time.time() - evolve_start

# Final water balance (DE only). Collective: compute on ALL ranks, print on
# rank 0 below. Identity: final volume = initial + boundary flux + fractional-
# step (rainfall/inlet/operator) inflow.
water_balance = None
try:
    water_balance = compute_water_balance(domain)
except Exception as _e:
    if myid == 0:
        emit(f'Water balance unavailable: {_e}')

# ---------------------------------------------------------------------------
# [6/6] Run summary (rank 0), to both the log and the terminal.
# ---------------------------------------------------------------------------
if myid == 0:
    section('RUN SUMMARY', 6, 6)
    emit('   Phase timings (s)')
    if mesh_build_time is not None:
        emit('     mesh construction   %10.2f' % mesh_build_time)
    if mesh_distribute_time is not None and numprocs > 1:
        emit('     distribute          %10.2f' % mesh_distribute_time)
    emit('     evolve              %10.2f' % evolve_time)
    emit('     total (wall)        %10.2f' % (time.time() - t0))

    if water_balance is not None:
        v0, fs, bf, vol, resid = water_balance
        denom = max(abs(vol), abs(v0) + abs(fs) + abs(bf), 1.0)
        rel = resid / denom
        verdict = 'OK conserved' if abs(rel) < 1e-6 else 'CHECK imbalance'
        emit('')
        emit('   Water balance (m^3)')
        emit('     initial volume          %14.2f' % v0)
        emit('     rainfall + inlets (FS)  %14.2f' % fs)
        emit('     net boundary flux (BF)  %14.2f' % bf)
        emit('     final volume            %14.2f' % vol)
        emit('     imbalance               %14.2f   (%.1e rel)  %s'
             % (resid, rel, verdict))

# ---------------------------------------------------------------------------
# Outputs (max-quantity CSVs, parallel SWW merge, GeoTIFF rasters)
# ---------------------------------------------------------------------------

section('OUTPUTS')
if myid == 0:
    emit('   writing to %s' % project.output_dir)

max_quantity_file_start = domain.get_datadir() + '/Max_quantities_'
max_quantities.export_max_quantities_to_csv(max_quantity_file_start)

os.chdir(project.output_dir)
if myid == 0 and numprocs > 1:
    print('Number of processors %g ' % numprocs)
    print('That took %.2f seconds' % (time.time() - t0))
    print('Communication time %.2f seconds' % domain.communication_time)
    print('Reduction Communication time %.2f seconds'
          % domain.communication_reduce_time)
    print('Broadcast time %.2f seconds'
          % domain.communication_broadcast_time)

    anuga.utilities.sww_merge.sww_merge_parallel(
        project.scenario,
        np=numprocs, verbose=True, delete_old=True)

if myid == 0:
    try:
        raster_outputs.make_me_some_tifs(
            sww_file='./' + project.scenario + '.sww',
            bounding_polygon=project.bounding_polygon,
            proj4string=project.proj4string,
            cell_size=project.output_tif_cellsize)
    except Exception as e:
        print('GeoTif creation failed: ' + str(e))
        print('You can try manually using raster_outputs.py or '
              'anuga.utilities.plot_utils.Make_Geotif')

barrier()
finalize()
