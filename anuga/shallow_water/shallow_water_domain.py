"""
Finite-volume computations of the shallow water wave equation.

Title: ANGUA shallow_water_domain - 2D triangular domains for finite-volume
       computations of the shallow water wave equation.


Author: Ole Nielsen, Ole.Nielsen@ga.gov.au
        Stephen Roberts, Stephen.Roberts@anu.edu.au
        Duncan Gray, Duncan.Gray@ga.gov.au
        Gareth Davies, gareth.davies.ga.code@gmail.com

CreationDate: 2004

Description:
    This module contains a specialisation of class Generic_Domain from
    module generic_domain.py consisting of methods specific to the
    Shallow Water Wave Equation

    U_t + E_x + G_y = S

    where

    U = [w, uh, vh]
    E = [uh, u^2h + gh^2/2, uvh]
    G = [vh, uvh, v^2h + gh^2/2]
    S represents source terms forcing the system
    (e.g. gravity, friction, wind stress, ...)

    and _t, _x, _y denote the derivative with respect to t, x and y
    respectively.


    The quantities are

    symbol    variable name    explanation

    x         x                horizontal distance from origin [m]
    y         y                vertical distance from origin [m]
    z         elevation        elevation of bed on which flow is modelled [m]
    h         height           water height above z [m]
    w         stage            absolute water level, w = z+h [m]
    u                          speed in the x direction [m/s]
    v                          speed in the y direction [m/s]
    uh        xmomentum        momentum in the x direction [m^2/s]
    vh        ymomentum        momentum in the y direction [m^2/s]


    eta                        mannings friction coefficient [to appear]
    nu                         wind stress coefficient [to appear]

    The conserved quantities are w, uh, vh

Reference:
    Catastrophic Collapse of Water Supply Reservoirs in Urban Areas,
    Christopher Zoppou and Stephen Roberts,
    Journal of Hydraulic Engineering, vol. 127, No. 7 July 1999

    Hydrodynamic modelling of coastal inundation.
    Nielsen, O., S. Roberts, D. Gray, A. McPherson and A. Hitchman
    In Zerger, A. and Argent, R.M. (eds) MODSIM 2005 International Congress on
    Modelling and Simulation. Modelling and Simulation Society of Australia and
    New Zealand, December 2005, pp. 518-523. ISBN: 0-9758400-2-9.
    http://www.mssanz.org.au/modsim05/papers/nielsen.pdf

    See also: https://anuga.anu.edu.au and https://en.wikipedia.org/wiki/ANUGA_Hydro


Constraints: See license in the user guide
"""

from __future__ import annotations




# Decorator added for profiling
#------------------------------


def profileit(name):
    def inner(func):
        def wrapper(*args, **kwargs):
            import cProfile
            prof = cProfile.Profile()
            retval = prof.runcall(func, *args, **kwargs)

            # Note use of name from outer scope
            print(str(args[1])+'_'+name)
            prof.dump_stats(str(args[1])+'_'+name)
            return retval
        return wrapper
    return inner
#-----------------------------

import warnings

import numpy as num
import sys
import os
import time
from typing import TYPE_CHECKING
from collections.abc import Callable, Iterator
from numpy.typing import ArrayLike

if TYPE_CHECKING:
    from datetime import datetime as DateTime
    from zoneinfo import ZoneInfo as ZoneInfoType
    # Only for set_collect_max_quantities()'s return annotation. The real import is
    # done inside that method, because operators import this module — a module-level
    # import here would be circular.
    from anuga.operators.collect_max_quantities_operator import (
        Collect_max_quantities_operator,
    )


try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo

try:
    import dill as pickle
except ImportError:
    import pickle

from anuga.abstract_2d_finite_volumes.generic_domain \
                    import Generic_Domain

from anuga.config import MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU
from anuga.config import LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2

from anuga.shallow_water.forcing import Cross_section
from anuga.utilities.numerical_tools import mean
from anuga.file.sww import SWW_file

import anuga.utilities.log as log

from anuga.utilities.parallel_abstraction import size, rank, get_processor_name
from anuga.utilities.parallel_abstraction import finalize, send, receive
from anuga.utilities.parallel_abstraction import pypar_available, barrier


# The centroid arrays the GPU interface round-trips between host and device —
# see gpu_domain_sync_to_device()/_from_device() in gpu/gpu_domain_core.c.  A host
# write to any other quantity (elevation, friction, a user tracer) is not
# device-resident state and so needs no sync.  Keep this in step with the C.
GPU_SYNCED_QUANTITIES = frozenset(('stage', 'xmomentum', 'ymomentum', 'height'))

# Absolute lower bound (m^3) on the volume added by clamping negative-depth cells
# below which the "possible loss of conservation" warning is never raised. This
# rejects pure floating-point noise (femto/pico-litre deficits in a nearly-dry
# domain) that can otherwise be a large *fraction* of an essentially-zero total
# volume. Any physically meaningful conservation loss is many orders larger.
_negative_volume_noise_floor = 1.0e-9


#-----------------------------------------------------
# Code for profiling cuda version
#-----------------------------------------------------
def nvtxRangePush(*arg):
    pass
def nvtxRangePop(*arg):
    pass

try:
    from cupy.cuda.nvtx import RangePush as nvtxRangePush
    from cupy.cuda.nvtx import RangePop  as nvtxRangePop
except ImportError:
    pass

try:
    from nvtx import range_push as nvtxRangePush
    from nvtx import range_pop  as nvtxRangePop
except ImportError:
    pass


#-----------------------------------------------------
# Process-global GPU offload control
#
# Whether mode 2 ('unified') offloads to a GPU is a *process-level* OpenMP
# setting (the target-offload runtime ICV), not a per-domain property: a single
# process cannot run one domain on the GPU and another on the CPU with the same
# unified kernels. These module functions own that process-wide decision; the
# per-domain choice (legacy vs unified) lives on Domain.set_compute_mode().
#-----------------------------------------------------

def _gpu_ext_or_none():
    try:
        from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext
        return gpu_ext
    except Exception:
        return None


def gpu_offload_supported() -> bool:
    """True if this build/run can offload mode 2 to a GPU device.

    Reflects build + hardware: False for a CPU-only build (``gpu_offload=false`` /
    ``CPU_ONLY_MODE`` — the standard pip/conda install), when no device is
    present, or when offload was disabled at launch via
    ``OMP_TARGET_OFFLOAD=disabled`` (``gpu_available()`` covers all three). This
    is a static capability and is *not* affected by :func:`set_gpu_offload`.
    """
    ge = _gpu_ext_or_none()
    return bool(ge.gpu_available()) if ge is not None else False


def gpu_offload_enabled() -> bool:
    """Return the resolved process-global offload state for mode 2 ('unified').

    True only when the build supports offload *and* it has not been switched off
    via :func:`set_gpu_offload`. Always False on a CPU-only build.
    """
    if not gpu_offload_supported():
        return False
    ge = _gpu_ext_or_none()
    return bool(ge.get_offload_enabled()) if ge is not None else False


def gpu_startup_banner(numprocs: int,
                       num_devices: int,
                       device_id: int,
                       offload_active: bool,
                       omp_num_threads: str = '1') -> list:
    """Build the mode-2 ('unified') startup banner as a list of lines.

    Pure function of its arguments so the rank/device mismatch warning can be tested
    without an actual multi-GPU machine.

    `numprocs` is the MPI rank count and `num_devices` the number of GPUs the runtime
    can actually see.  These are DIFFERENT NUMBERS, and the banner used to print the
    rank count labelled as "GPU(s)" — so a 4-rank run on a 1-GPU box cheerfully
    reported "4 GPU(s)".  That mattered because ranks are assigned to devices
    round-robin (``device_id = rank % num_devices``), so an oversubscribed run silently
    puts several ranks on one device.  The banner was the natural place to notice that,
    and instead it concealed it.  See issue #194.

    Why the warning is worded the way it is.  Measured on one RTX 5070, a 160k-triangle
    mode-2 evolve, ranks all sharing the single device:

        ranks   no MPS     with MPS
          1      3.80 s     3.51 s
          2      7.13 s     3.48 s
          4     11.21 s     3.72 s

    So the *reliable* cost of oversubscribing is speed: without MPS the ranks time-slice
    the device and it is ~3x slower.  NVIDIA MPS lets their kernels run concurrently and
    restores parity — but never beats one rank per GPU, because a single rank with the
    whole mesh already saturates the device; splitting it creates no new parallelism.
    So MPS is a way to stop losing, not a way to go faster.

    Hangs and garbled results have been reported in this configuration on real hardware,
    but did NOT reproduce in the runs above (every rank count gave a bit-identical
    checksum).  The warning therefore leads with the measured slowdown and reports the
    hangs as a possibility, rather than promising a failure that may not arrive — an
    oversubscribed run that quietly works but is 3x slow is the likelier outcome, and is
    exactly the one a user would otherwise never notice.
    """

    bar = '+==============================================================================+'
    lines = [bar]

    if not offload_active:
        lines.append("| ANUGA compute mode: 'unified' CPU multicore (gpu_ext kernels, no offload)   |")
        lines.append(f'| OMP_NUM_THREADS={omp_num_threads}')
    elif device_id < 0:
        lines.append('| WARNING: No GPU devices found, running on CPU via OpenMP target offloading  |')
    elif num_devices <= 0:
        # Device count unavailable (query failed). Say so rather than inventing one —
        # printing numprocs here is exactly the bug this function exists to remove.
        lines.append(f'| GPU interface initialized: {numprocs} MPI rank(s), device count unknown, '
                     f'OpenMP target offloading')
    else:
        lines.append(f'| GPU interface initialized: {numprocs} MPI rank(s) on {num_devices} GPU(s), '
                     f'OpenMP target offloading')

        if numprocs > num_devices:
            lines.append(f'| WARNING: {numprocs} MPI ranks but only {num_devices} GPU(s). Ranks map round-robin')
            lines.append(f'|          (rank % {num_devices}), so several ranks share one device. Mode-2 MPI')
            lines.append('|          expects ONE RANK PER GPU. Measured: up to ~3x SLOWER than one')
            lines.append('|          rank per GPU (the ranks time-slice the device); hangs and wrong')
            lines.append('|          results have also been reported in this configuration.')
            lines.append(f'|          Re-run with -np {num_devices}. NVIDIA MPS removes the slowdown')
            lines.append('|          but never beats one rank per GPU, so it is not a way to go faster.')
        elif numprocs < num_devices:
            idle = num_devices - numprocs
            lines.append(f'| NOTE: {idle} GPU(s) idle — {numprocs} rank(s) for {num_devices} device(s). '
                         f'Use -np {num_devices} to use them all.')

    lines.append(bar)
    return lines


# Boundary classes the unified (mode-2) path evaluates on the device, by class
# name. Everything else is evaluated on the host each step and synced back,
# which works but costs a host<->device round trip per step. This is the ONE
# list; it used to be copied into nine methods with three different spellings
# of the Flather class, so Flather boundaries silently took the host path in
# most of them.
GPU_BOUNDARY_TYPES = frozenset({
    'Reflective_boundary', 'Dirichlet_boundary', 'Transmissive_boundary',
    'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary',
    'Time_boundary', 'File_boundary', 'Field_boundary',
    'Absorbing_wave_boundary', 'Characteristic_wave_boundary',
    'Flather_external_stage_zero_velocity_boundary',
})


def _strict_compute_mode_from_env() -> bool:
    """ANUGA_STRICT_COMPUTE_MODE=1 makes every domain strict by default."""
    return os.environ.get('ANUGA_STRICT_COMPUTE_MODE', '').strip().lower() in ('1', 'true', 'yes', 'on')


def set_gpu_offload(enable: bool = True, verbose: bool = True, strict: bool = False) -> bool:
    """Enable or disable GPU offload for mode-2 ('unified') domains, process-wide.

    This is a *process-level* switch, not per-domain: it affects every domain
    that runs in 'unified' mode, because OpenMP target offload is a process-wide
    runtime setting (one process cannot run the same unified kernels on a GPU for
    one domain and on the CPU for another).

    The setting is honoured at GPU-domain init, where arrays are mapped to the
    chosen device. **Call it before building the first 'unified' domain** — once
    a domain's data is mapped to a device, switching this domain's offload would
    leave data and execution on different devices.

    Implemented via the OpenMP default-device ICV (``omp_set_default_device`` /
    a process-global flag in ``sw_domain_gpu_ext``) — robust and re-enableable,
    unlike mutating ``OMP_TARGET_OFFLOAD`` after the runtime has initialised.

    Parameters
    ----------
    enable : bool
        True to offload 'unified' domains to a GPU (GPU build + device required);
        False to force 'unified' to run CPU-multicore.
    verbose : bool
        Print a one-line confirmation.

    strict : bool
        If True, requesting ``enable=True`` on a build or machine that cannot
        offload raises ``RuntimeError`` instead of warning and falling back to
        the CPU. Use it in batch jobs, where a "GPU job" that quietly runs on
        the CPU is the expensive failure. ``ANUGA_STRICT_COMPUTE_MODE=1`` in
        the environment has the same effect.

    Returns
    -------
    bool
        The resolved offload state (:func:`gpu_offload_enabled`). Requesting
        ``enable=True`` on a build without offload support warns and returns
        False — never hard-fails unless ``strict``.
    """
    import warnings

    ge = _gpu_ext_or_none()
    was_supported = gpu_offload_supported()
    strict = strict or _strict_compute_mode_from_env()

    if enable:
        # Clear any prior disable BEFORE checking capability, so re-enabling is
        # not blocked by our own OMP_TARGET_OFFLOAD=disabled.
        os.environ.pop('OMP_TARGET_OFFLOAD', None)
        if not gpu_offload_supported():
            msg = ("set_gpu_offload(True): this ANUGA build has no GPU offload support "
                   "(built with gpu_offload=false) or no device is present; 'unified' "
                   "domains will run on CPU multicore. Rebuild with -Dgpu_offload=true "
                   "and a GPU-capable compiler to enable offload.")
            if strict:
                raise RuntimeError(msg + " (strict mode: refusing to fall back to the CPU)")
            warnings.warn(msg, stacklevel=2)
            enable = False
    elif was_supported:
        # Disabling offload on a GPU build. This forces the unified kernels onto
        # the host, but a GPU build (nvc -mp=gpu,multicore) runs the `omp target`
        # regions on the host through a slow fallback that does NOT scale with
        # threads — so this is for correctness A/B (GPU vs CPU give the same
        # results), NOT for performance. For fast CPU multicore, use a
        # gpu_offload=false (gcc) build instead.
        try:
            from anuga import myid
        except Exception:
            myid = 0
        if myid == 0:
            warnings.warn(
                "set_gpu_offload(False) on a GPU build: 'unified' will run on the "
                "host, but the nvc GPU build's host fallback is slow and does not "
                "scale with threads (use it for correctness checks, not timing). "
                "For fast CPU multicore, build with -Dgpu_offload=false.",
                stacklevel=2)

    if not enable:
        # OMP_TARGET_OFFLOAD=disabled is what actually keeps the `omp target`
        # regions (solver AND operators) on the host. The default-device flag
        # alone does NOT stop the solver offloading, so it must be set here.
        # The OpenMP runtime reads this at its first target region, so call
        # set_gpu_offload() before the first evolve()/domain build.
        os.environ['OMP_TARGET_OFFLOAD'] = 'disabled'

    # Keep the C-side flag consistent so gpu_domain_init picks the host device
    # and the inlet/culvert operators route there too (gpu_compute_device).
    if ge is not None:
        ge.set_offload_enabled(bool(enable))

    state = bool(enable)
    if verbose:
        print(f"GPU offload {'enabled' if state else 'disabled'} "
              f"(process-wide; 'unified' domains run on {'GPU' if state else 'CPU multicore'})")
    return state


# Process-wide OpenMP thread count for ANUGA kernels. This is the single source
# of truth read back by ``Domain.omp_num_threads`` (a property), so that setting
# it once via ``anuga.set_omp_num_threads(n)`` is reflected by every domain in
# the session — including ones already constructed (important in notebooks).
# Initialised from OMP_NUM_THREADS so introspection is sane before the first
# call / domain construction.
try:
    _omp_num_threads = int(os.environ.get('OMP_NUM_THREADS', 1))
except (ValueError, TypeError):
    _omp_num_threads = 1


def get_omp_num_threads() -> int:
    """Return the current process-wide OpenMP thread count for ANUGA kernels."""
    return _omp_num_threads


def set_omp_num_threads(omp_num_threads: int | None = None, verbose: bool = True) -> int:
    """Set the OpenMP thread count for ANUGA kernels (process-wide).

    ``OMP_NUM_THREADS`` / ``omp_set_num_threads`` controls the whole process, so
    this is a module-level setting, not per-domain — it affects every domain's
    OpenMP regions (both the legacy ``sw_domain_openmp_ext`` solver and the
    unified ``gpu_ext`` kernels), and is reflected by the ``omp_num_threads``
    property of every existing domain. ``Domain.set_omp_num_threads`` delegates
    here.

    Parameters
    ----------
    omp_num_threads : int or None
        Thread count. If None, use ``OMP_NUM_THREADS`` from the environment,
        defaulting to 1 when unset.
    verbose : bool
        Print a one-line confirmation.

    Returns
    -------
    int
        The thread count applied.
    """
    if omp_num_threads is None:
        omp_num_threads = os.environ.get('OMP_NUM_THREADS', None)
        if verbose:
            print(f'Using OMP_NUM_THREADS from environment: {omp_num_threads}')

    if omp_num_threads is None:
        omp_num_threads = 1  # Default to 1 if not set

    try:
        omp_num_threads = int(omp_num_threads)
    except (ValueError, TypeError):
        raise ValueError('OMP_NUM_THREADS must be an integer')

    # omp_set_num_threads is process-global: one call covers every OpenMP region
    # in the process, so routing through the legacy extension also sets the
    # thread count for the unified gpu_ext kernels.
    from .sw_domain_openmp_ext import set_omp_num_threads as set_omp_num_threads_ext
    set_omp_num_threads_ext(omp_num_threads)
    # Keep the env var consistent so banners / introspection / any subprocess
    # report the same count (the runtime ICV is already set above).
    os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)
    # Record the process-wide count so every domain's omp_num_threads property
    # reflects this call (including domains constructed before it).
    global _omp_num_threads
    _omp_num_threads = omp_num_threads

    if verbose:
        print(f'Setting omp_num_threads to {omp_num_threads}')
    return omp_num_threads


# [D-5] default thickness of the near-bed layer as a fraction of the depth.
# 0.2 is what van Rijn's trench and Wang & Ribberink's settling flume are
# reproduced with; 0 falls back to the reference height a.
TWO_LAYER_FRACTION = 0.2

class Domain(Generic_Domain):
    """Object which encapulates the shallow water model


    This class is a specialization of class Generic_Domain from
    module generic_domain.py consisting of methods specific to the
    Shallow Water Wave Equation

    Shallow Water Wave Equation

    .. math::

        U_t + E_x + G_y = S

    where

    .. math::

        U = [w, uh, vh]^T

    .. math::

        E = [uh, u^2h + gh^2/2, uvh]

    .. math::

        G = [vh, uvh, v^2h + gh^2/2]



    S represents source terms forcing the system
    (e.g. gravity, friction, wind stress, ...)

    and _t, _x, _y denote the derivative with respect to t, x and y
    respectively.


    The quantities are

    .. list-table::
        :widths: 25 25 50
        :header-rows: 1

        * - symbol
          - variable name
          - explanation
        * - x
          - x
          - horizontal distance from origin [m]
        * - y
          - y
          - vertical distance from origin [m]
        * - z
          - elevation
          - elevation of bed on which flow is modelled [m]
        * - h
          - height
          - water height above z [m]
        * - w
          - stage
          - absolute water level, w = z+h [m]
        * - u
          -
          - speed in the x direction [m/s]
        * - v
          -
          - speed in the y direction [m/s]
        * - uh
          - xmomentum
          - momentum in the x direction [m^2/s]
        * - vh
          - ymomentum
          - momentum in the y direction [m^2/s]
        * -
          -
          -
        * - eta
          -
          - mannings friction coefficient [to appear]
        * - nu
          -
          - wind stress coefficient [to appear]

    The conserved quantities are w, uh, vh

    """

    def __init__(self,
                 coordinates: ArrayLike | str | None = None,
                 vertices: ArrayLike | None = None,
                 boundary: dict | None = None,
                 tagged_elements: dict | None = None,
                 geo_reference=None,
                 use_inscribed_circle: bool = False,
                 mesh_filename: str | None = None,
                 use_cache: bool = False,
                 verbose: bool = False,
                 conserved_quantities: list[str] | None = None,
                 evolved_quantities: list[str] | None = None,
                 other_quantities: list[str] | None = None,
                 full_send_dict: dict | None = None,
                 ghost_recv_dict: dict | None = None,
                 starttime: float = 0,
                 processor: int = 0,
                 numproc: int = 1,
                 number_of_full_nodes: int | None = None,
                 number_of_full_triangles: int | None = None,
                 ghost_layer_width: int = 2,
                 **kwargs) -> None:

        """Instantiate a shallow water domain.

        :param coordinates: vertex locations for the mesh
        :param vertices: vertex indices defining the triangles of the mesh
        :param boundary: boundaries of the mesh
        """

        # Define quantities for the shallow_water domain
        if conserved_quantities is None:
            conserved_quantities = ['stage', 'xmomentum', 'ymomentum']

        if evolved_quantities is None:
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum']

        if other_quantities is None:
            other_quantities = ['elevation', 'friction', 'height',
                                'xvelocity', 'yvelocity', 'x', 'y']

        # Selective array allocation per quantity type.
        # Quantities not listed default to 'evolved' (all arrays) for
        # backward compatibility with user-defined quantities.
        self._quantity_type_map = {
            'stage':       'evolved',
            'xmomentum':   'evolved',
            'ymomentum':   'evolved',
            # elevation: centroid + edge only; gradients lazy (erosion operators
            # trigger allocation via compute_local_gradients when needed)
            'elevation':   'edge_diagnostic',
            'friction':    'centroid_only',
            'height':      'edge_diagnostic',
            'xvelocity':   'edge_diagnostic',
            'yvelocity':   'edge_diagnostic',
            'x':           'coordinate',
            'y':           'coordinate',
        }


        Generic_Domain.__init__(self,
                            coordinates,
                            vertices,
                            boundary,
                            conserved_quantities,
                            evolved_quantities,
                            other_quantities,
                            tagged_elements,
                            geo_reference,
                            use_inscribed_circle,
                            mesh_filename,
                            use_cache,
                            verbose,
                            full_send_dict,
                            ghost_recv_dict,
                            starttime,
                            processor,
                            numproc,
                            number_of_full_nodes=number_of_full_nodes,
                            number_of_full_triangles=number_of_full_triangles,
                            ghost_layer_width=ghost_layer_width)

        #-------------------------------
        # Operator Data Structures
        #-------------------------------
        self.fractional_step_operators = []
        self.kv_operator = None
        self.max_quantities_operator = None
        self.dplotter = None
        self.gpu_culvert_manager = None  # Initialized when GPU mode + Boyd operators

        #-------------------------------
        # Set flow defaults
        #-------------------------------
        self.set_flow_algorithm()

        #-------------------------------
        # Set default multiprocessor mode
        # 1. Openmp
        # 2. Cuda
        #-------------------------------
        self.gpu_interface = None
        self.use_c_rk_loop = True  # Use C RK loop (faster) vs Python-orchestrated GPU loop
        # Default compute mode: 'legacy' (mode 1). Set ANUGA_DEFAULT_COMPUTE_MODE=unified
        # to default new domains to mode 2 (the migration target); the device interface
        # is then built lazily at first evolve() so construction needs no boundaries.
        # SERIAL ONLY: under MPI (numprocs > 1) the env opt-in is ignored and domains
        # stay 'legacy'. Mode-2 parallel is validated for purpose-built setups (it must
        # be selected explicitly), but is not yet robust as a blanket default for every
        # parallel test's evolve pattern — defaulting all parallel domains to mode 2
        # deadlocks the MPI path. Parallel unified is a later migration step.
        try:
            from anuga import numprocs
        except Exception:
            numprocs = 1
        if (numprocs == 1
                and os.environ.get('ANUGA_DEFAULT_COMPUTE_MODE', 'legacy').lower() == 'unified'):
            self.set_compute_mode('unified')
        else:
            self.set_compute_mode('legacy')

        #-------------------------------
        # C extension domain structure
        # Will be setup by setup_domain_openmp_ext
        #-------------------------------
        self._Domain_C_struct = None

        #-------------------------------
        # Generic passive tracers (see add_tracer).
        #
        # Ns = 0 is the default and costs nothing: the kernels read
        # number_of_tracers and skip the tracer work entirely, and the six
        # arrays below are passed to the C struct as NULL.
        #-------------------------------
        self.number_of_tracers = 0
        self.beta_tracer = 1.0
        self._tracer_names = []
        # (tracer index, boundary tag) -> callable, for time-varying
        # inflow concentrations. See set_tracer_boundary().
        self._tracer_boundary_functions = {}
        self.tracer_centroid_values = None      # c            (ns, N)
        self.tracer_edge_values = None          # c at edges   (ns, 3N)
        self.tracer_boundary_values = None      # c at bdry    (ns, boundary_length)
        self.tracer_explicit_update = None      # dm/dt        (ns, N)
        self.tracer_conserved_values = None     # m = h*c      (ns, N)
        self.tracer_backup_values = None        # m backup     (ns, N)
        self.tracer_external_source = None      # S_ms [G-3]   (ns, N) [m/s]
                                                # allocated by add_tracer, with
                                                # the other tracer blocks

        #-------------------------------
        # Phase 3: suspended sediment. A sediment class IS a tracer -- class s
        # occupies tracer slot s -- with per-class settling parameters. Zero
        # classes costs nothing: the source kernel returns on one test.
        #-------------------------------
        self.n_sediment_classes = 0
        self.sediment_c_max = 0.30              # [L-2]; FG21 0.30, aS16 0.20
        self._sediment_names = []
        # rho_s and the settling kwargs are kept per fraction so that
        # R and v_s can be recomputed if the domain-wide rho_w changes
        # after registration -- otherwise they would silently go stale.
        self._sediment_rho_s = []
        self._sediment_settling_kwargs = []
        self.sediment_gamma0 = 0.0024           # [E-1] empirical, FG21
        # Erosion route, spec 4.1.1 -- see set_bed_material().
        self.sediment_erosion_mode = 0          # 0 non-cohesive, 1 cohesive
        self.sediment_tau_crit = 0.088          # [Pa], aS16's value
        self.sediment_K_e = 0.2e-6 / 0.088**0.5  # [E-5]
        self.sediment_rho_w = 1000.0
        self.sediment_K_partheniades = 1.0e-4   # [E-4] kg m-2 s-1 (RDy26)
        # [E-6] de Leeuw et al. (2020) entrainment, erosion mode 3; the
        # constants of their Eq 26a until set_bed_material says otherwise.
        (self.sediment_dl_A, self.sediment_dl_alpha, self.sediment_dl_beta,
         self.sediment_dl_threshold) = self.DE_LEEUW_FITS['de_leeuw_2020']
        self.sediment_dl_ks = 0.0               # <= 0: ks_factor * diameter
        self.sediment_dl_ks_factor = 3.0        # k_s = 3 D84
        self.sediment_deposition_mode = 0       # 0 = [D-1], 1 = [D-2]
        self.sediment_tau_d = 0.0               # [D-2] critical depo stress [Pa]
        # Bed shear closure, spec 3.1/3.4 (divergence D1). 0 = [T-1] quadratic
        # drag (default); 1 = [T-7] depth-slope, for reproducing anugaSed.
        self.sediment_shear_closure = 0
        self.sediment_max_slope = 0.0           # cap on S for [T-7]/[T-7e]; 0 = none
        self.sediment_slope_frozen = 0          # depth-slope S fixed at setup
        self.sediment_settling_velocity = None  # v_s        (ncl,)
        self.sediment_d_star = None             # d*(Z)      (ncl,)
        self.sediment_diameter = None           # d_g   [m]  (ncl,)
        self.sediment_R = None                  # R          (ncl,)
        self.sediment_tau_c_star = None         # tau_c*     (ncl,)
        self.sediment_reference_height = None   # a     [m]  (ncl,)
        self.sediment_d_star_mode = 0           # 0 constant, 1 Rouse [S-4]
        # [D-5] the two-layer suspension with the log-law velocity split is
        # the DEFAULT closure: against van Rijn's trench and pick-up flume
        # and Wang & Ribberink's settling flume the instantaneous exchange
        # deposits at the equilibrium stratification d*, which is 2 to 4
        # times the measured adaptation rate for w_s/u* of 0.2 to 1.
        # adaptation='none' restores it.
        self.sediment_adaptation_mode = 4       # 0 none, 1 Armanini, 2 constant, 3 carried, 4 two-layer
        self.sediment_adaptation_alpha = 1.0    # [D-3] alpha for mode 2
        self.sediment_nearbed_base = -1        # [D-4] first near-bed tracer, or -1
        self.sediment_layer_fraction = TWO_LAYER_FRACTION   # [D-5] near-bed layer h1/h; 0: a/h
        self.sediment_exchange_factor = 1.0    # [D-5] factor on the partition's relaxation rate
        self.sediment_velocity_profile = 1     # [D-5v] advect the layers at their log-law speeds
        self.tracer_speed_factor = None        # [D-5v] (n_tracers, n) speed factors, or None
        # van Rijn-style floor a >= sediment_a_h_floor * h, applied when
        # sediment_d_star_mode = 1. Standard practice, on by default. Set to 0
        # to reach anugaSed's regime (they use no floor); the d* fit covers
        # a/h down to 1e-3. This is the largest single divergence from
        # anugaSed -- roughly 8x less deposition at h = 1 m. See spec 12, D4b.
        self.sediment_a_h_floor = 0.01
        self.sediment_rouse_beta_mode = 0      # [S-2b] 1: Z / (1 + 2 (w_s/u*)^2)
        self.sediment_rouse_scale = 1.0        # [S-2b] factor on Z
        # [L-4] maximum packing fraction bounding the near-bed concentration
        # c_b = d* c. Without it the equilibrium Rouse d* makes the deposition
        # rate diverge as shear vanishes. Same constant that bounds E* in
        # [E-1]. Inactive when d* = 1, since c <= c_max = 0.3 < 0.65.
        self.sediment_c_pack = 0.65
        # [G-4] bed porosity lambda: the sediment VOLUME leaving suspension is
        # (1-lambda) dz, the remainder being pore space filled from the water
        # column. LM15 Example 2 uses 0.28.
        self.sediment_porosity = 0.30
        self.sediment_morphological_factor = 1.0   # M, bed change per step x M
        # Coupling stage, spec 2.4. True = evolving bed via [G-4] (Phase 4);
        # False = FIXED bed (Phase 3), which is RDy26 v1.0's configuration and
        # what the analytic constant-depth deposition solutions assume. Both are
        # published configurations, not a debug switch.
        self.sediment_bed_evolution = True
        # Bedload, spec 6. Off by default; see set_bedload().
        self.sediment_bedload_mode = 0
        self.sediment_bedload_K = 3.97
        self.sediment_bedload_m = 1.5
        self.sediment_bedload_tau_c_star = 0.0495
        self.sediment_bedload_h_min = 0.0          # [K-7] bedload depth ramp, 0 = off
        self.sediment_qbx = None
        self.sediment_qba = None
        self.sediment_qby = None
        # [L-5] non-erodible base, spec 4.5. Off by default: with no base the
        # bed is bottomless, which is what every published test case in the
        # spec assumes. See set_erodible_base().
        self.sediment_z_base = None
        self.sediment_has_z_base = 0
        self.sediment_shear_factor = None      # [T-16] per-centroid factor on the sediment's bed shear
        # The two user intents behind sediment_z_base, kept apart so they
        # compose: a base is a DEPTH limit, a region is a WHERE limit, and
        # setting one must not silently discard the other. Both are folded
        # into the single field the kernel reads by _rebuild_sediment_base().
        self._sediment_user_base = None      # (n,) from set_erodible_base
        self._sediment_erodible_mask = None  # (n,) bool from set_erodible_region
        # Bed composition (active layer), off by default: the bed is then one
        # shared material and every class may be entrained from it. See
        # set_bed_composition().
        self.sediment_bed_composition = 0
        self.sediment_active_layer = 0.0
        self.sediment_bed_active = None      # (ncl, n) solid volume per area [m]
        self.sediment_bed_substrate = None   # (ncl, n)
        # Spec 7, angle-of-repose relaxation. Off by default: it is a numerical
        # heuristic, not physics, and it suppresses knickpoint retreat that may
        # be real. See set_angle_of_repose().
        self.sediment_repose_tan = 0.0
        self.sediment_repose_relax = 1.0
        self.sediment_repose_max_sweeps = 50
        self.sediment_repose_dz = None
        self.sediment_bed_exhausted = None
        # Which boundary edges pass bedload (set_bedload open_boundaries):
        # int64 per boundary edge, 1 = open. The tags are kept so the array
        # can be rebuilt when it is first allocated or the mesh is reordered.
        self.sediment_bedload_open = None
        self.sediment_bedload_supply = None    # per boundary edge, m2/s, <0 = none
        self._sediment_bedload_supply = {}     # tag -> m2/s
        self._sediment_bedload_open_tags = ()
        # Scratch for the source kernel, (ncl, n). Allocated with the classes.
        self.sediment_source_limited = None
        # Scratch for the source kernel, (n): the per-cell slope of the
        # depth-slope closures, computed in a pass of its own.
        self.sediment_slope_work = None
        # Friction closure for the sediment kernel (spec 3.3). 'constant' is
        # the right default for ordinary flood work; see set_sediment_friction.
        self.sediment_friction_mode = 0
        self.sediment_manning_ll = 0.065
        self.sediment_wilson_bed = 0
        self.sediment_wilson_D = 1.0e-3
        self.tracer_boundary_flux = None        # d(mass)/dt across the
                                                # domain boundary, per cell,
                                                # accumulated by the kernel (ns, N)
        # Per-substep totals of the above, (max_time_substeps * ns), mirroring
        # boundary_flux_sum for water; and the running time integral, which
        # tracer_flux_integral_operator advances with the timestepping method's
        # own weights. _tracer_initial_mass is the baseline a conservation
        # check is measured against.
        self.tracer_boundary_flux_sum = None
        self._tracer_flux_integral = None
        self._tracer_initial_mass = None

        #-------------------------------
        # If environment variable OMP_NUM_THREADS is not set,
        # then set to default (1 thread). If a value is given to
        # the method, then it will override the default.
        #------------------------------
        self.set_omp_num_threads(verbose=False)

        #-------------------------------
        # datetime and timezone
        #-------------------------------
        self.set_timezone()

        #-------------------------------
        # Forcing Terms
        #
        # Gravity is now incorporated in
        # compute_fluxes routine
        #-------------------------------
        from .friction import manning_friction_semi_implicit
        self.forcing_terms.append(manning_friction_semi_implicit)

        # Vegetation drag (spec 8), off until set_vegetation_drag()
        self.vegetation_mode = 0
        self.vegetation_Cd = 1.68
        self.vegetation_bed_chezy = 65.0
        self.sediment_vegetation_shear = 1       # bed share; see set_vegetation_drag


        #-------------------------------
        # Stored output
        #-------------------------------
        self.set_store(True)
        self.set_store_centroids(True)
        self.set_store_vertices_uniquely(False)
        self.quantities_to_be_stored = {'elevation': 1,
                                        'friction':1,
                                        'stage': 2,
                                        'xmomentum': 2,
                                        'ymomentum': 2}

        #-------------------------------
        # Set up check pointing every n
        # yieldsteps
        #-------------------------------
        self.checkpoint = False
        self.yieldstep_counter = 0
        self.checkpoint_step = 10

        #-------------------------------
        # Useful auxiliary quantity
        # Set centroid and edge values directly from mesh coordinate arrays
        # without allocating vertex_values (kept lazy to save memory).
        # Edge formula matches _interpolate in quantity_openmp.c:
        #   edge[0] = 0.5*(v[1]+v[2]),  edge[1] = 0.5*(v[2]+v[0]),
        #   edge[2] = 0.5*(v[0]+v[1])
        #-------------------------------
        n = self.number_of_elements
        vx = self.vertex_coordinates[:, 0].reshape(n, 3)
        vy = self.vertex_coordinates[:, 1].reshape(n, 3)

        qx = self.quantities['x']
        qx.centroid_values[:]  = (vx[:, 0] + vx[:, 1] + vx[:, 2]) / 3.0
        qx.edge_values[:, 0]   = 0.5 * (vx[:, 1] + vx[:, 2])
        qx.edge_values[:, 1]   = 0.5 * (vx[:, 2] + vx[:, 0])
        qx.edge_values[:, 2]   = 0.5 * (vx[:, 0] + vx[:, 1])
        qx.set_boundary_values_from_edges()

        qy = self.quantities['y']
        qy.centroid_values[:]  = (vy[:, 0] + vy[:, 1] + vy[:, 2]) / 3.0
        qy.edge_values[:, 0]   = 0.5 * (vy[:, 1] + vy[:, 2])
        qy.edge_values[:, 1]   = 0.5 * (vy[:, 2] + vy[:, 0])
        qy.edge_values[:, 2]   = 0.5 * (vy[:, 0] + vy[:, 1])
        qy.set_boundary_values_from_edges()

        # For riverwalls, we need to know the 'edge_flux_type' for each edge
        # Edge-flux-type of 0 == Normal edge, with shallow water flux
        #                   1 == riverwall
        #                   2 == ?
        #                   etc
        # Lazy: allocated by _ensure_work_arrays() (no river walls) or by
        # create_riverwalls() (before evolve).  Saves ~108 MB at N=2.25M.
        self.edge_flux_type = None          # (3N,) int — 0=normal, 1=riverwall
        self.edge_river_wall_counter = None # (3N,) int — per-edge riverwall index
        self.number_of_riverwall_edges = 0

        # Riverwalls -- initialise with dummy values
        # Presently only works with DE algorithms, will fail otherwise
        import anuga.structures.riverwall
        self.riverwallData=anuga.structures.riverwall.RiverWall(self)
        self.create_riverwalls = self.riverwallData.create_riverwalls

        ## Keep track of the fluxes through the boundaries
        ## Only works for DE algorithms at present
        max_time_substeps=3 # Maximum number of substeps supported by any timestepping method
        # boundary_flux_sum holds boundary fluxes on each sub-step [unused substeps = 0.]
        self.boundary_flux_sum=num.array([0.]*max_time_substeps)
        from anuga.operators.boundary_flux_integral_operator import boundary_flux_integral_operator
        self.boundary_flux_integral=boundary_flux_integral_operator(self)
        # Make an integer counting how many times we call compute_fluxes_central -- so we know which substep we are on
        #self.call=1

        # List to store the volumes we computed before
        self.volume_history=[]

        # Work arrays actually read by the C extension — allocated lazily on the
        # first evolve step (via setup_Domain_C_struct → _ensure_work_arrays).
        self.x_centroid_work     = None  # (N,)   scratch for velocity extrapolation
        self.y_centroid_work     = None  # (N,)   scratch for velocity extrapolation

        # The following arrays were historically allocated but are confirmed dead
        # (C computation never reads them).  They remain None (→ NULL in C struct)
        # for the entire simulation lifetime, saving ~600+ MB at N=2.25M.
        self.edge_flux_work      = None  # (9N,)  — unused, NULL in C struct
        self.neigh_work          = None  # (9N,)  — unused, NULL in C struct
        self.pressuregrad_work   = None  # (3N,)  — unused, NULL in C struct


        #-----------------------------------
        # parameters for structures
        #-----------------------------------
        self.use_new_velocity_head = False


    #------------------------------------------------
    # Domain_C_struct is a cdef class with a custom __cinit__,
    # so Cython will not auto-generate a default pickling protocol for it;
    # when pickle reaches the Domain object and tries to pickle _Domain_C_struct,
    # you get TypeError: no default __reduce__ due to non-trivial __cinit__.
    # So we implement __getstate__ and __setstate__ to exclude it from pickling,
    # and recreate it lazily when needed.
    #------------------------------------------------
    def __getstate__(self):
        state = self.__dict__.copy()
        # Do not pickle the C wrapper; it can be recreated
        state.pop('_Domain_C_struct', None)
        # The mode-2 ('unified') GPU interface wraps a cdef GPUDomain with a
        # non-trivial __cinit__ (not picklable) and holds device handles. Drop
        # it; _ensure_gpu_interface() rebuilds it lazily after unpickling.
        state.pop('gpu_interface', None)
        state.pop('_gpu_boundary_info_initialized', None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Recreate C wrapper lazily when needed
        self._Domain_C_struct = None
        # Force the mode-2 device interface to be rebuilt on demand.
        self.gpu_interface = None

    #------------------------------------------------
    # Generic passive tracers
    #------------------------------------------------
    # Every tracer block, and how it is INDEXED -- which is what fixes both its
    # shape and the way reorder() has to permute it:
    #
    #   'cell'      one value per triangle       -> (ns, N)
    #   'edge'      three values per triangle    -> (ns, 3N)
    #   'boundary'  one value per boundary edge  -> (ns, boundary_length)
    #
    # Single source of truth on purpose: add_tracer sizes the arrays from this
    # and _reorder_tracer_arrays permutes from it, so a new block cannot be
    # added without declaring how it is indexed. Tracers are deliberately not
    # Quantity objects (#276), which means anything walking domain.quantities
    # misses them -- #277 was exactly that omission, in reorder().
    _TRACER_ARRAY_KINDS = {
        'tracer_centroid_values':  'cell',
        'tracer_edge_values':      'edge',
        'tracer_boundary_values':  'boundary',
        'tracer_explicit_update':  'cell',
        'tracer_conserved_values': 'cell',
        'tracer_backup_values':    'cell',
        'tracer_boundary_flux':    'cell',
        # The external source S_ms [G-3]. Allocated with the rest rather than
        # on first use: a lazily allocated array has to change the C struct and
        # the device mapping when it appears, and set_tracer_source did that by
        # discarding the GPU interface -- mid-run, from inside a fractional
        # step. Sediment_transport_operator then saw gpu_interface is None, took the CPU
        # path while the state was on the device, and the source contributed
        # exactly nothing on a GPU build (#288). Being here also means reorder()
        # permutes it, which it previously did not.
        'tracer_external_source':  'cell',
    }
    _TRACER_ARRAYS = tuple(_TRACER_ARRAY_KINDS)

    def add_tracer(self, name, beta=None, initial_value=0.0):
        """Register a passive tracer named `name` and return its index.

        A tracer is a depth-averaged concentration `c` advected with the water
        flux.  The conserved variable is `m = h*c`; `c` is derived from it each
        substep, exactly as height is derived from stage.

        Parameters
        ----------
        name : str
            Identifier for the tracer.  Must be unique on this domain.
        beta : float, optional
            Edge-reconstruction limiter coefficient.  0 selects first order,
            > 0 a limited second-order reconstruction.  **This is a single
            value shared by every tracer** (the C struct carries one
            `beta_tracer` scalar), so passing a value that disagrees with an
            already-registered tracer is an error rather than a silent
            last-writer-wins.  Defaults to leaving the current value alone.
        initial_value : float or array-like, optional
            Initial concentration `c`.  A scalar fills the domain; an array
            must have one value per cell.  `m = h*c` is seeded consistently.

        Returns
        -------
        int
            The tracer's index, i.e. its row in the `(ns, ...)` arrays.

        Notes
        -----
        Registering a tracer **reallocates** every tracer array, so any
        reference held to one of them beforehand becomes stale.  Add every
        tracer before seeding values, or re-fetch via `get_tracer`.
        """
        if getattr(self, 'sediment_nearbed_base', -1) >= 0:
            raise ValueError(
                'no tracer can be added once the layered-adaptation near-bed '
                'tracers are registered (at the first evolve): they must stay '
                'the last tracers. Add every tracer before evolving.')
        if not isinstance(name, str) or not name:
            raise ValueError('tracer name must be a non-empty string')
        if name in self._tracer_names:
            raise ValueError(
                'a tracer named %r is already registered (index %d)'
                % (name, self._tracer_names.index(name)))
        # A tracer is written to the sww as <name>_c, which is exactly the
        # variable a quantity's centroid values go to. A tracer named 'stage'
        # would therefore OVERWRITE the stage in the output -- no error, and a
        # file whose stage is silently something else. Reserve the names.
        if name in self.quantities:
            raise ValueError(
                'a tracer cannot be named %r: that is a quantity on this '
                'domain, and both are written to the sww as %s_c, so the '
                'tracer would overwrite it' % (name, name))
        if name.startswith('max_'):
            raise ValueError(
                "a tracer cannot be named %r: names beginning 'max_' are "
                'reserved for the running maxima that '
                'Collect_max_quantities_operator writes' % name)

        if beta is not None:
            beta = float(beta)
            if beta < 0.0:
                raise ValueError('beta must be >= 0, got %g' % beta)
            # One scalar serves all tracers, so a second, different value would
            # silently change the reconstruction of the tracers already added.
            if self.number_of_tracers > 0 and beta != self.beta_tracer:
                raise ValueError(
                    'beta_tracer is shared by all tracers on a domain: cannot '
                    'add %r with beta=%g while existing tracers use beta=%g'
                    % (name, beta, self.beta_tracer))
            self.beta_tracer = beta

        N = self.number_of_elements
        ns = self.number_of_tracers
        widths = {'cell': N, 'edge': 3 * N, 'boundary': self.boundary_length}
        shapes = {attr: (ns + 1, widths[kind])
                  for attr, kind in self._TRACER_ARRAY_KINDS.items()}

        # Grow each array by one row, preserving the tracers already there.
        # The kernels index these as centroid[s*N + k] etc., so they must stay
        # C-contiguous float64 -- num.zeros gives both.
        # The per-substep totals and the running integral are sized by tracer
        # count, not by cells, so they are rebuilt here rather than in the loop
        # below. The integral and the baseline mass are PRESERVED across a
        # later add_tracer: a tracer registered mid-setup must not silently
        # reset another tracer's accounting.
        max_substeps = len(self.boundary_flux_sum)
        new_sum = num.zeros((ns + 1) * max_substeps, dtype=num.float64)
        new_int = num.zeros(ns + 1, dtype=num.float64)
        new_m0 = num.zeros(ns + 1, dtype=num.float64)
        if ns > 0 and self._tracer_flux_integral is not None:
            new_int[:ns] = self._tracer_flux_integral
            new_m0[:ns] = self._tracer_initial_mass
        self.tracer_boundary_flux_sum = new_sum
        self._tracer_flux_integral = new_int
        self._tracer_initial_mass = new_m0

        for attr in self._TRACER_ARRAYS:
            new = num.zeros(shapes[attr], dtype=num.float64)
            if ns > 0:
                old = getattr(self, attr)
                new[:ns] = old
            setattr(self, attr, new)

        index = ns
        self._tracer_names.append(name)
        self.number_of_tracers = ns + 1

        # THE TRAP: the C struct is built once and cached, and evolve() never
        # passes update_domain_c_struct=True. Without this line a tracer
        # registered after the struct exists is invisible to the kernels --
        # no error, the tracer simply never moves. It also now holds pointers
        # into the arrays we just replaced. Invalidate so the next call
        # rebuilds it against the new arrays.
        self._Domain_C_struct = None

        # Same argument on the device side: the GPU interface has the OLD
        # tracer arrays mapped (or none at all, if it was built at Ns=0), and
        # the arrays it points at have just been freed. Tear it down so it is
        # rebuilt and re-mapped against the new ones. _ensure_gpu_interface()
        # recreates it on demand.
        self.gpu_interface = None
        # NB: this flag is tested with hasattr, not for truthiness (see
        # update_boundary), so it must be DELETED, not set False -- setting it
        # False would skip the re-initialisation that defines _gpu_all_on_gpu.
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

        if initial_value is not None:
            self.set_tracer(name, initial_value)

        # The operator that turns the kernel's per-substep boundary totals into
        # a time integral. Registered here rather than left to the user: without
        # it get_tracer_boundary_flux_integral() would silently read zero, which
        # is exactly the failure a conservation check exists to catch.
        from anuga.operators.tracer_flux_integral_operator import (
            tracer_flux_integral_operator)
        if not any(isinstance(op, tracer_flux_integral_operator)
                   for op in self.fractional_step_operators):
            tracer_flux_integral_operator(self)

        # Baseline for check_tracer_conservation, taken AFTER the initial value
        # is seeded so it is the mass actually present at t = 0.
        self._tracer_initial_mass[index] = self._local_tracer_mass(index)

        return index

    def settling_velocity(self, diameter, rho_s=2650.0, rho_w=1000.0,
                          nu=1.0e-6, C1=18.0, C2=0.4):
        """Settling velocity `v_s` from Ferguson & Church (2004), spec `[S-1]`.

            v_s = R g d^2 / ( C1 nu + sqrt(0.75 C2 R g d^3) )

        Smooth across the Stokes/turbulent transition and branch-free, which is
        why the spec prefers it over the Dietrich (1982) polynomial fit.

        `C1 = 18, C2 = 0.4` are the smooth-sphere constants; use `1.0, 1.1` for
        natural irregular grains. `R = rho_s/rho_w - 1` is submerged specific
        gravity.

        Verified against the spec: d = 4.5e-5 m quartz gives 1.75e-3 m/s,
        the value P13 report for Ferguson & Church.
        """
        import math
        from anuga.config import g as _g
        R = rho_s / rho_w - 1.0
        d = float(diameter)
        if d <= 0.0:
            raise ValueError('grain diameter must be > 0, got %g' % d)
        return (R * _g * d * d) / (C1 * nu + math.sqrt(0.75 * C2 * R * _g * d**3))

    def _allocate_sediment_fraction(self, name, diameter, d_star=1.0, beta=None,
                                    initial_concentration=0.0, rho_s=2650.0,
                                    tau_c_star=0.04,
                                    reference_height=None,
                                    **settling_kwargs):
        """Register one suspended sediment fraction and return its index.

        Private: the public entry point is
        :class:`~anuga.operators.sediment_operator.Sediment_transport_operator`,
        which calls this. Kept here because the bookkeeping it does -- growing
        the per-fraction parameter arrays and invalidating the C struct and the
        device mapping -- belongs with the rest of the domain's array
        management.

A sediment fraction is a tracer -- so it is transported by the machinery of
        Phases 1-2 -- plus the settling parameters the source term needs. The
        tracer is registered first, so fraction `s` always occupies tracer
        slot `s`; `add_tracer` and `add_sediment_fraction` must not be interleaved on
        the same domain if you rely on that.

        Parameters
        ----------
        name : str
            Class identifier, e.g. 'sand'. Also the tracer name.
        diameter : float
            Grain diameter `d_g` in metres. Settling velocity is computed once
            here via `settling_velocity` (`[S-1]`) rather than per cell.
        d_star : float, optional
            Ratio of near-bed to depth-averaged concentration in `[D-1]`.
            Default 1.0, the well-mixed limit. The Rouse profile of spec 4.3
            replaces this constant later.
        rho_s, rho_w : float, optional
            Sediment and water densities. Only their ratio matters: the
            submerged specific gravity `R = rho_s/rho_w - 1` is what enters the
            Shields stress, where water density cancels.
        tau_c_star : float, optional
            Critical Shields stress for entrainment `[E-1]`. Default 0.04,
            FG21's choice for suspension. Setting it to 0 disables entrainment
            for this class, leaving deposition only.
        reference_height : float, optional
            `a` in `[S-4]`, the near-bed reference height at which `c_b` is
            evaluated, in metres. Only used when `sediment_d_star_mode = 1`.
            Defaults to `2*diameter`. **This is a first-order choice, not a
            detail**: `d*` varies by up to 13x across plausible `a/h` at high
            Rouse number. aSM16 requires `a` but never states it, so the default
            here is `2*diameter`, which the D4b audit of `anugaSed` independently
            corroborates as their convention too. The kernel additionally applies
            the floor `a >= domain.sediment_a_h_floor * h` (default 0.01); set
            that to 0 to match `anugaSed`, which applies no floor.
        initial_concentration : float or array-like, optional
            Initial `c_s`; seeds `m = h*c` consistently.

        Notes
        -----
        Phase 3 is the FIXED-BED stage (spec 2.4): the bed does not evolve, and
        there is no bed->flow or sediment->momentum feedback. Entrainment draws
        from an inexhaustible bed and deposited mass leaves the system; the
        Exner bookkeeping of `[G-4]` is Phase 4.

        Entrainment uses the NON-COHESIVE Shields route `[E-1]`/`[E-2]`, which
        is a statement about the bed material (sand/gravel), not a numerical
        preference -- see spec 4.1.1. The cohesive Hanson & Simon route
        `[E-3]` is for silt and clay and is not implemented here.
        """
        if self.sediment_nearbed_base >= 0:
            raise ValueError(
                "with a layered adaptation every fraction must be added before "
                'the first evolve: its tracers are registered then')
        if self.number_of_tracers != self.n_sediment_classes:
            raise ValueError(
                'add_sediment_fraction requires fraction s to occupy tracer '
                'slot s, but this domain already has %d tracers and %d sediment '
                'fractions. Do not mix add_tracer() and add_sediment_fraction().'
                % (self.number_of_tracers, self.n_sediment_classes))

        if tau_c_star < 0.0:
            raise ValueError('tau_c_star must be >= 0, got %g' % tau_c_star)
        # rho_w is a property of the fluid, so it is domain-wide: see
        # set_sediment_parameters. Only rho_s varies per fraction.
        rho_w = self.sediment_rho_w
        v_s = self.settling_velocity(diameter, rho_s=rho_s, rho_w=rho_w,
                                     **settling_kwargs)
        if self.sediment_bed_composition:
            raise RuntimeError(
                'add every sediment fraction before set_bed_composition(): the '
                'bed composition is held per class, so it cannot gain a class '
                'afterwards')
        index = self.add_tracer(name, beta=beta,
                                initial_value=initial_concentration)

        ncl = self.n_sediment_classes + 1
        for attr, value in (('sediment_settling_velocity', v_s),
                            ('sediment_d_star', float(d_star)),
                            ('sediment_diameter', float(diameter)),
                            ('sediment_R', rho_s / rho_w - 1.0),
                            ('sediment_tau_c_star', float(tau_c_star)),
                            ('sediment_reference_height',
                             float(reference_height) if reference_height
                             is not None else 2.0 * float(diameter))):
            new = num.zeros(ncl, dtype=num.float64)
            if self.n_sediment_classes > 0:
                new[:self.n_sediment_classes] = getattr(self, attr)
            new[index] = value
            setattr(self, attr, new)

        if self.sediment_qbx is None:
            self.sediment_qbx = num.zeros(self.number_of_elements,
                                          dtype=num.float64)
            self.sediment_qby = num.zeros(self.number_of_elements,
                                          dtype=num.float64)
            self.sediment_qba = num.zeros(self.number_of_elements,
                                          dtype=num.float64)
            # [L-5] snapshot, int64 to match anuga_int.
            self.sediment_bed_exhausted = num.zeros(self.number_of_elements,
                                                    dtype=num.int64)
            # Spec 7 Jacobi scratch.
            self.sediment_repose_dz = num.zeros(self.number_of_elements,
                                                dtype=num.float64)
        self._build_bedload_open()

        # Scratch for the source kernel: every class's bed exchange is held
        # here until [L-5] has limited them together. Sized (ncl, n) and
        # REALLOCATED as classes are added -- the kernel dereferences it
        # whenever n_sediment_classes > 0, so it must never be short.
        self.sediment_source_limited = num.zeros(
            (ncl, self.number_of_elements), dtype=num.float64)
        if self.sediment_slope_work is None:
            self.sediment_slope_work = num.zeros(self.number_of_elements,
                                                 dtype=num.float64)
        self._freeze_bed_slope()

        self._sediment_names.append(name)
        self._sediment_rho_s.append(float(rho_s))
        self._sediment_settling_kwargs.append(dict(settling_kwargs))
        self.n_sediment_classes = ncl

        # add_tracer already invalidated both caches, but it did so BEFORE the
        # arrays above existed. Invalidate again so the rebuilt struct sees them.
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

        # The bed can now move, so elevation must be stored per timestep.
        self._sync_elevation_storage()

        return index

    def _sync_elevation_storage(self):
        """Store elevation per timestep once the bed can move.

        The bedload and bed-exchange kernels write the bed in place --
        `bed_centroid_values` and `bed_edge_values` ARE the `elevation`
        Quantity's arrays (see sw_domain_openmp_ext.pyx) -- so an evolving bed
        is a time-varying quantity. But `quantities_to_be_stored` defaults to
        `'elevation': 1`, and flag 1 means write it ONCE (sww.py sorts flag 1
        into `static_quantities`). The sww then records the initial bed and
        silently omits every change to it.

        That is worse than an obviously missing variable: the output looks
        complete, and the bed in it is plausible but wrong. So upgrade to
        flag 2 as soon as a domain has both sediment classes and an evolving
        bed.

        Only the default is upgraded. Flags 3 and 4 are deliberate choices
        (centroid-only, and overwrite-each-yieldstep), so they are left alone
        rather than silently rewritten. The flag is never downgraded either:
        turning bed evolution back off leaves elevation stored dynamically,
        which costs space but cannot mislead.
        """
        if self.n_sediment_classes <= 0 or not self.sediment_bed_evolution:
            return
        if self.quantities_to_be_stored.get('elevation') == 1:
            self.quantities_to_be_stored['elevation'] = 2

    def set_angle_of_repose(self, angle=None, relax=1.0, max_sweeps=50):
        """Relax bed slopes steeper than `angle` by moving material downslope.

        Spec 7, from FG21 §2.2.4. Where the centroid-to-centroid bed slope
        exceeds the critical angle, material is diffused downslope until it
        does not.

            domain.set_angle_of_repose(35.0)     # degrees; FG21 use 35
            domain.set_angle_of_repose(None)     # off again (the default)

        Parameters
        ----------
        angle : float
            Critical angle in DEGREES, in (0, 90). `None` or 0 disables it.
        relax : float
            Relaxation in (0, 1], default 1.0. The kernel already divides by
            the edge count for stability, so 1.0 is the fastest STABLE setting
            rather than an aggressive one -- measured 793 sweeps to converge an
            over-steep cone against 2400+ at 0.3. Lower it only if you see
            something pathological.
        max_sweeps : int
            Hard cap on sweeps per timestep, default 50. Reaching it is
            reported, because it means the bed may still be over-steep.

        Notes
        -----
        FG21 are explicit that this is **a numerical heuristic, not physics**:
        real bed slope failures are advective. It exists to stop the rest of
        the model breaking on over-steep slopes. It has a side effect worth
        knowing before you switch it on -- it limits the steepness of canyon
        walls and knickpoints, and so suppresses knickpoint retreat that may be
        real. That is why it is off by default.

        Mass is conserved: material removed from an over-steep cell is
        deposited on its neighbour, never discarded. This is the one place this
        module differs sharply from `sanddune_erosion_operator`, which lowers
        an over-steep cell and lets the material vanish.

        Respects `[L-5]`: a cell cannot slump away material it is not allowed
        to lose, so a locked cell or one at its base stays put and its
        neighbours relax around it.

        **On the sweep count.** This is an explicit diffusion solve, so
        convergence from a badly over-steep bed is slow: an over-steep cone
        needed 793 sweeps to reach the critical angle from cold. That is not
        what the per-timestep cap is sized for. In a running model the bed is
        already near-relaxed and each step needs a handful of sweeps; the cap
        is there for the pathological case, and hitting it is not fatal --
        progress carries over, so the bed keeps relaxing on subsequent steps.
        It is reported so that you know relaxation is lagging rather than
        finished.

        If you START from a bed steeper than the critical angle, expect the
        cap to be hit on the first steps. Either accept that it settles over
        the first few, or raise `max_sweeps` for that run.
        """
        if angle is None or angle == 0.0:
            self.sediment_repose_tan = 0.0
        else:
            if not 0.0 < angle < 90.0:
                raise ValueError(
                    'angle of repose must be in (0, 90) degrees, got %g'
                    % angle)
            self.sediment_repose_tan = float(num.tan(num.radians(angle)))

        if not 0.0 < relax <= 1.0:
            raise ValueError('relax must be in (0, 1], got %g' % relax)
        if max_sweeps < 1:
            raise ValueError('max_sweeps must be >= 1, got %d' % max_sweeps)
        self.sediment_repose_relax = float(relax)
        self.sediment_repose_max_sweeps = int(max_sweeps)

        if (self.sediment_repose_tan > 0.0
                and self.sediment_repose_dz is None):
            self.sediment_repose_dz = num.zeros(self.number_of_elements,
                                                dtype=num.float64)

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_shear_amplification(self, factor=None):
        """Multiply the bed shear the sediment sees by a per-cell factor `[T-16]`.

        The depth-averaged shear is the uniform-flow one and misses the
        local amplification at obstacles: at a bridge pier the horseshoe
        vortex raises the bed shear to two to four times the approach value
        within about a diameter, and that is what digs the scour hole. This
        method takes the factor as a field so that a case or a structure
        operator can put it where the physics is missing::

            domain.set_shear_amplification(lambda x, y: 1 + 3 * np.exp(-(r(x, y) - R) / R))
            domain.set_shear_amplification(None)     # off again (the default)

        Parameters
        ----------
        factor : None, float, array (n,) or callable(x, y)
            Per-centroid factor >= 0, applied to `f_c` in both the suspended
            exchange and the bedload kernel (so to `u_*`, the Shields stress,
            the Rouse number and every closure built on them). `None` or 1
            removes it.
        """
        if factor is None:
            self.sediment_shear_factor = None
        else:
            x = self.centroid_coordinates[:, 0]
            y = self.centroid_coordinates[:, 1]
            if callable(factor):
                f = num.asarray(factor(x, y), dtype=num.float64)
            else:
                f = num.asarray(factor, dtype=num.float64)
            f = num.broadcast_to(f, (self.number_of_elements,)).astype(num.float64).copy()
            if f.min() < 0.0:
                raise ValueError('the shear amplification must be >= 0')
            self.sediment_shear_factor = num.ascontiguousarray(f)
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_erodible_base(self, elevation=None, depth=None):
        """Set the non-erodible base -- bedrock -- below which no erosion acts.

        `[L-5]`. Without this the bed is bottomless: erosion lowers it for as
        long as the flow has the strength to, which is the right default for
        an alluvial channel and wrong wherever the erodible layer is finite --
        a scoured reach over an outcrop, a lined culvert, a dam apron, a soil
        layer of known depth over rock.

        The base is a per-CENTROID elevation, not a scalar, because bedrock is
        a surface. Give it either way round::

            domain.set_erodible_base(elevation=z_rock)  # absolute [m]
            domain.set_erodible_base(depth=0.5)         # 0.5 m below the
                                                        # elevation set so far
            domain.set_erodible_base()                  # remove the base

        Parameters
        ----------
        elevation : float or array (n,)
            Base elevation, in the same datum as the domain's elevation
            quantity. Scalar broadcasts.
        depth : float or array (n,)
            Erodible thickness below the CURRENT bed. The base is recorded as
            an elevation at the moment of the call, so later changes to the
            elevation quantity do not move it. Scalar broadcasts.

        Give exactly one. With neither, the base is removed and the bed is
        bottomless again.

        Notes
        -----
        The limiter acts on the SOURCE, never by clamping elevation. Erosion
        is scaled back to what the remaining thickness can supply, and the
        sediment that is not eroded never enters the water column, so the
        budget still closes exactly. Clamping z afterwards would leave
        suspended sediment that came from nowhere.

        Where several classes compete for the last of the material, they are
        scaled by a shared proportional factor rather than served in order:
        the bed carries no per-class stratigraphy, so no class has a better
        claim, and registration order must not change the answer. Deposition
        is never scaled -- it is what replenishes the bed.

        Bedload is limited too, and stays exactly conservative while it is:
        the limit applies to the transport vector and to whole edges, both of
        which the two cells sharing an edge see identically.
        """
        if elevation is not None and depth is not None:
            raise ValueError('give elevation or depth, not both')

        if elevation is None and depth is None:
            self._sediment_user_base = None
        else:
            n = self.number_of_elements
            z = self.quantities['elevation'].centroid_values
            if depth is not None:
                d = num.asarray(depth, dtype=num.float64)
                if num.any(d < 0.0):
                    raise ValueError('erodible depth must be >= 0')
                base = z - d
            else:
                base = num.asarray(elevation, dtype=num.float64)

            base = num.ascontiguousarray(
                num.broadcast_to(base, (n,)), dtype=num.float64).copy()

            # A base above the bed is not a configuration, it is a mistake:
            # the cell starts with negative erodible thickness and the
            # limiter would simply hold it there, silently.
            over = base - z
            if num.any(over > 0.0):
                worst = int(num.argmax(over))
                raise ValueError(
                    'erodible base is ABOVE the bed in %d of %d cells '
                    '(worst: cell %d, base %g > elevation %g). The base is '
                    'the floor of erosion, so it must lie at or below the '
                    'initial bed everywhere.'
                    % (int(num.sum(over > 0.0)), n, worst,
                       base[worst], z[worst]))

            self._sediment_user_base = base

        self._rebuild_sediment_base()

    def set_erodible_region(self, region=None, polygon=None, center=None,
                            radius=None, indices=None, erodible=True):
        """Restrict erosion to part of the domain (or lock part of it).

        The bed is erodible everywhere by default. Give a region to say
        otherwise::

            domain.set_erodible_region(polygon=breach)     # ONLY here erodes
            domain.set_erodible_region(polygon=apron,      # everywhere BUT here
                                       erodible=False)
            domain.set_erodible_region(my_region)          # a Region object
            domain.set_erodible_region()                   # remove the restriction

        Parameters
        ----------
        region : Region
            An already-built `anuga.abstract_2d_finite_volumes.region.Region`.
            This is the general form: `Region` also understands `line=`,
            `poly=` and `expand_polygon=`, which have no keyword here, so
            build one and pass it when you need them. It must belong to THIS
            domain -- its indices mean nothing on another mesh, and one built
            elsewhere is rejected rather than quietly mis-selecting.
        polygon : list of [x, y]
            Region boundary, as for the region-based operators.
        center, radius : [x, y], float
            A circular region instead.
        indices : array of int
            Triangle ids directly. Overrides the geometric arguments.
        erodible : bool
            `True` (default): the region named is the ONLY erodible part.
            `False`: the region named is the only LOCKED part.

        The keyword arguments are the same ones `Erosion_operator` and the
        other region-based operators take, and are resolved by the same
        `Region` class, so a polygon that selects a set of cells there selects
        the same set here.

        Notes
        -----
        A locked cell is held at the elevation it has WHEN THIS IS CALLED, by
        giving it zero erodible thickness -- the region restriction is
        `[L-5]` with the layer set to nothing, not a separate mechanism. So
        call it after the elevation is set.

        Locked means it cannot be SCOURED. Sediment may still settle onto it,
        which is what a concrete apron or a rock bar does in the field, and
        that new material is erodible again -- it is above the base. If you
        want a cell that neither erodes nor accretes, that is not this.

        Composes with `set_erodible_base()`: the base sets how DEEP erosion
        may go, the region sets WHERE it may happen, and setting one leaves
        the other in place. Where they disagree the stricter wins.
        """
        from anuga.abstract_2d_finite_volumes.region import Region

        if (region is None and polygon is None and center is None
                and radius is None and indices is None):
            self._sediment_erodible_mask = None
            self._rebuild_sediment_base()
            return

        if region is not None:
            if any(a is not None
                   for a in (polygon, center, radius, indices)):
                raise ValueError(
                    'give a Region or the arguments to build one, not both')
            if not isinstance(region, Region):
                # A list of points is the likely mistake, and it is a silent
                # one: Region would not be consulted and every cell would look
                # selected.
                raise TypeError(
                    'region must be a Region object; to pass a list of points '
                    'use set_erodible_region(polygon=...)')
            if getattr(region, 'domain', None) is not self:
                raise ValueError(
                    'that Region belongs to a different domain; its triangle '
                    'indices do not refer to this mesh')
        else:
            region = Region(self, indices=indices, polygon=polygon,
                            center=center, radius=radius)
        idx = region.indices

        n = self.number_of_elements
        if idx is None:
            # Region resolved to "everywhere".
            sel = num.ones(n, dtype=bool)
        else:
            idx = num.asarray(idx, dtype=num.int64)
            if idx.size == 0:
                # Almost always a coordinate mistake -- a polygon in the wrong
                # units or the wrong datum selects nothing, and the run then
                # quietly does no erosion at all (or, with erodible=False, is
                # unrestricted). Neither is what anyone meant.
                raise ValueError(
                    'the region selects no cells; check the polygon or centre '
                    'is in the same coordinates as the mesh')
            sel = num.zeros(n, dtype=bool)
            sel[idx] = True

        # erodible=True: the region is the erodible part. Otherwise it is the
        # locked part and everything else erodes.
        self._sediment_erodible_mask = sel if erodible else ~sel
        self._rebuild_sediment_base()

    # Depth used for "no limit" in the combined base field. The cap it implies
    # (thickness (1-lambda)/dt) is then so far above any physical erosion rate
    # that it never binds, which is what an absent base means.
    _SEDIMENT_UNLIMITED_DEPTH = 1.0e6

    def _rebuild_sediment_base(self):
        """Fold the base and the region into the one field the kernel reads."""
        base = self._sediment_user_base
        mask = self._sediment_erodible_mask

        if base is None and mask is None:
            self.sediment_z_base = None
            self.sediment_has_z_base = 0
        else:
            n = self.number_of_elements
            z = self.quantities['elevation'].centroid_values
            if base is None:
                combined = z - self._SEDIMENT_UNLIMITED_DEPTH
            else:
                combined = base.copy()
            if mask is not None:
                # Locked cells get zero thickness: the base IS the bed.
                combined[~mask] = z[~mask]
            self.sediment_z_base = num.ascontiguousarray(combined,
                                                         dtype=num.float64)
            self.sediment_has_z_base = 1
            if self.sediment_bed_exhausted is None:
                self.sediment_bed_exhausted = num.zeros(n, dtype=num.int64)

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def erodible_thickness(self):
        """Remaining erodible thickness per centroid [m], or None if no base.

        `elevation - sediment_z_base`, clipped at zero. Zero means the cell
        has reached bedrock and will not erode further.
        """
        if not self.sediment_has_z_base:
            return None
        z = self.quantities['elevation'].centroid_values
        return num.maximum(z - self.sediment_z_base, 0.0)

    def set_bed_composition(self, fractions=None, active_layer=0.01,
                            substrate=None):
        """Track the bed's composition per sediment class (an active layer).

        Without this the bed is one shared material: any class can be
        entrained from it wherever the shear allows, so mud keeps coming out
        of a sand-bedded channel for as long as the bed lasts. With it the bed
        holds an amount of each class in a thin active layer over a substrate
        (Hirano), and:

        * entrainment of a class is multiplied by its share of the active
          layer, so a class the bed does not hold is not entrained, and a
          bed that loses its fines armours;
        * a class can never take more out of a cell than the cell holds of it;
        * after each step the active layer is restored to `active_layer`:
          on aggradation the overflow goes down at the active composition, on
          degradation it is refilled from the substrate at the substrate
          composition.

        Each class is conserved exactly (the layers hold amounts, solid volume
        per bed area, not fractions).

        Call it AFTER every add_sediment_fraction, initialize_sediment_operator
        (for the porosity) and set_erodible_base (for the thickness, which it
        divides into the two layers)::

            domain.set_erodible_base(depth=1.0)
            domain.set_bed_composition({'mud': fmud, 'sand': 1 - fmud},
                                       active_layer=0.01)

        Parameters
        ----------
        fractions : dict
            Class name -> volume fraction of the active layer: a scalar, a
            per-centroid array, or a function f(x, y) of the centroid
            coordinates (as in set_quantity). Classes not given are 0. At
            every centroid the fractions must sum to 1 (within 1e-6); cells
            where the erodible thickness is 0 hold nothing whatever is given.
            ``None`` switches composition off again.
        active_layer : float
            Active-layer thickness [m] (bed, i.e. including pores). The
            Caltech Wax Lake model uses 0.01 m.
        substrate : dict, optional
            The substrate's composition, same form. Defaults to `fractions`.

        Notes
        -----
        Bedload and the angle-of-repose relaxation still move the bed as one
        material (phase 1); combining them with composition warns.
        """
        if fractions is None:
            self.sediment_bed_composition = 0
            self.sediment_active_layer = 0.0
            self.sediment_bed_active = None
            self.sediment_bed_substrate = None
            self._invalidate_sediment_caches()
            return

        ncl = self.n_sediment_classes
        if ncl == 0:
            raise RuntimeError('add the sediment fractions before set_bed_composition()')
        if not self.sediment_has_z_base:
            raise RuntimeError(
                'set_bed_composition() needs a finite erodible thickness to '
                'divide into the active layer and the substrate: call '
                'set_erodible_base() first')
        if not (active_layer > 0.0):
            raise ValueError('active_layer must be > 0')
        if substrate is None:
            substrate = fractions

        n = self.number_of_elements
        names = self.get_sediment_names()

        def fraction_table(spec, what):
            unknown = set(spec) - set(names)
            if unknown:
                raise ValueError('%s: unknown sediment class(es) %s; known: %s'
                                 % (what, sorted(unknown), names))
            F = num.zeros((ncl, n), dtype=num.float64)
            x = self.centroid_coordinates[:, 0]
            y = self.centroid_coordinates[:, 1]
            for s, name in enumerate(names):
                v = spec.get(name, 0.0)
                v = v(x, y) if callable(v) else v
                F[s] = num.broadcast_to(num.asarray(v, dtype=num.float64), (n,))
            if num.any(F < 0.0):
                raise ValueError('%s: fractions must be >= 0' % what)
            total = F.sum(axis=0)
            bad = num.abs(total - 1.0) > 1e-6
            if num.any(bad):
                k = int(num.argmax(bad))
                raise ValueError('%s: fractions sum to %g at cell %d (and %d '
                                 'others); they must sum to 1'
                                 % (what, total[k], k, int(bad.sum()) - 1))
            return F

        Fa = fraction_table(fractions, 'fractions')
        Fs = fraction_table(substrate, 'substrate')

        one_minus_lambda = 1.0 - self.sediment_porosity
        thick = self.erodible_thickness()
        t_act = num.minimum(thick, active_layer)
        t_sub = thick - t_act
        self.sediment_bed_active = num.ascontiguousarray(Fa * t_act * one_minus_lambda)
        self.sediment_bed_substrate = num.ascontiguousarray(Fs * t_sub * one_minus_lambda)
        self.sediment_active_layer = float(active_layer)
        self.sediment_bed_composition = 1

        if self.sediment_bedload_mode or self.sediment_repose_tan > 0.0:
            warnings.warn(
                'bed composition is tracked by the suspended exchange only: '
                'bedload and the angle-of-repose relaxation still move the bed '
                'as one material, so the layers and the elevation drift apart '
                'where they act', stacklevel=2)
        self._invalidate_sediment_caches()

    def get_bed_composition(self, layer='active'):
        """Volume fraction of each class in the active layer or the substrate.

        Returns {class name: (n,) array}; cells holding nothing in that layer
        give 0 for every class. The host copy is current at yieldsteps (on a
        GPU it is synced back with the bed).
        """
        if not self.sediment_bed_composition:
            return None
        A = self.sediment_bed_active if layer == 'active' else self.sediment_bed_substrate
        total = A.sum(axis=0)
        safe = num.where(total > 0.0, total, 1.0)
        return {name: num.where(total > 0.0, A[s] / safe, 0.0)
                for s, name in enumerate(self.get_sediment_names())}

    def _invalidate_sediment_caches(self):
        """New sediment arrays must reach the C struct and the device."""
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_deposition(self, law='d_star', tau_d=0.0, near_bed='constant',
                       reference_height_floor=0.01, adaptation='two_layer',
                       adaptation_alpha=None, layer_fraction=None,
                       exchange_factor=1.0, velocity_profile=None,
                       rouse_beta='none', rouse_scale=1.0):
        """Select the deposition law and its near-bed treatment (spec 4.4).

        Parameters
        ----------
        law : {'d_star', 'threshold'}
            `'d_star'` (default) -- `[D-1]`, `D = d*(Z) c v_s`.
            `'threshold'` -- `[D-2]`, `D = v_s c (1 - tau_b/tau_d)` for
            `tau_b < tau_d`, else zero; RDycore-sediment's form.
        tau_d : float
            Critical deposition stress in **pascals**, for `'threshold'`.
            `tau_d = 0` disables deposition entirely -- the hook RDycore's
            passive-transport benchmarks rely on.
        near_bed : {'constant', 'rouse'}
            How `d*` in `[D-1]` is obtained. `'constant'` uses the
            per-fraction value given to `add_sediment_fraction` (default 1.0, the well-mixed
            limit of P14/P13). `'rouse'` evaluates the fitted `[S-4]` profile
            per cell from the local Rouse number.
        reference_height_floor : float
            The van Rijn-style floor `a >= floor * h`, used by `'rouse'` and
            by the adaptation lag. Default 0.01. Set to 0 to reach anugaSed's
            regime, which applies no floor -- see spec 12, D4b, where this
            is the largest single divergence from them.
        adaptation : {'two_layer', 'carried', 'armanini', 'constant', 'none'}
            How the near-bed concentration follows a change of flow.
            `'two_layer'` (**the default**, `[D-5]`) splits the column into a
            near-bed layer and the rest, each fraction carrying its
            upper-layer mass in a tracer `<name>_upper`; deposition is `v_s`
            times the lower-layer concentration, and the exchange between
            the layers is set so the equilibrium is still the Rouse ratio
            `d*`. `'carried'` (`[D-4]`) carries only the ratio
            `r_b = c_b / c`, in a tracer `<name>_nearbed_ratio`, and relaxes
            it toward `d*` over a settling time when the flow slows.
            `'armanini'` (`[D-3]`) keeps no state and scales the exchange
            rate by Galappatti and Vreugdenhil's `alpha` in Armanini and Di
            Silvio's closed form. `'constant'` uses `adaptation_alpha`
            everywhere. `'none'` deposits at `d* c` the instant the flow
            changes, which was the default up to ANUGA 4.0.

            All of them leave every equilibrium concentration untouched and
            change only the transient. The layered ones register their
            tracers at the first `evolve`, so every fraction must be added
            before it. See :ref:`sediment_deposition` for which to choose,
            and the Notes below for why the default changed.
        adaptation_alpha : float, optional
            `alpha` for `adaptation='constant'`; must be > 0.
        layer_fraction : float, optional
            `'two_layer'` only: the near-bed layer thickness as a fraction
            of the depth, in [0, 0.5], default 0.2. Zero takes the reference
            height `a/h` with its floor instead. The exchange is re-derived
            so the equilibrium stays the Rouse ratio whatever the thickness.
        rouse_beta : {'none', 'van_rijn'}, optional
            Correction to the Rouse number of the `'rouse'` fit. `'van_rijn'`
            divides `Z` by `1 + 2 (w_s/u*)^2` (at most 2), van Rijn's
            (1984b) allowance for sediment mixing more strongly than
            momentum, which flattens the profile and lowers `d*`.
        rouse_scale : float, optional
            Factor on the (corrected) Rouse number, default 1. Van Rijn's
            trench profiles need about 0.65 on top of `'van_rijn'`.
        velocity_profile : bool, optional
            `'two_layer'` only (`[D-5v]`): advect each layer at its log-law
            mean velocity rather than the depth-averaged one, so the
            near-bed layer, which holds most of the sediment, lags the flow
            and the transport is the profile integral of `u c` rather than
            `u h c`. Default `None`, meaning on with `'two_layer'` and off
            otherwise. A case that sets an inflow concentration to carry a
            known load must divide it by the factor this applies.
        exchange_factor : float, optional
            `'two_layer'` only: a factor on the rate at which the partition
            relaxes toward its equilibrium, default 1 (the two-box
            estimate). Changes the transient only.

        Notes
        -----
        **Why the default is not `'none'`.** `[D-1]` deposits at `d* c`, and
        `d*` is the *equilibrium* stratification of the suspension, not the
        rate at which it reaches that equilibrium. Using one for the other
        deposits too fast, by a factor that grows with the settling velocity
        relative to the shear: 1.2 at `w_s/u* = 0.05`, 2.0 at 0.2, 4.1 at
        1.0. Three flume datasets agree. In Wang and Ribberink's (1986)
        settling flume, which has deposition and nothing else, the measured
        decay coefficient is 1.44 and 1.55 where their own theory gives 1.47
        and `'none'` gives 2.99; van Rijn's pick-up flume reaches its
        equilibrium load within about 15 depths where the flume took more
        than 40, and his trench fills about 25 % too fast (issue #389).
        `'two_layer'` with the velocity split gives 1.50 on the first, and
        the trench and pick-up flume improve with it too.

        Results are unchanged wherever `d* = 1`, the well-mixed limit that
        `near_bed='constant'` gives by default, because the partition
        reduces to `[D-1]` exactly there.

        See :ref:`sediment_deposition` for the closures side by side, and
        the sediment validation cases under `validation_tests/sediment/` for
        the measurements.
        """
        laws = {'d_star': 0, 'threshold': 1}
        if law not in laws:
            raise ValueError('unknown deposition law %r; expected one of %r'
                             % (law, sorted(laws)))
        modes = {'constant': 0, 'rouse': 1}
        if near_bed not in modes:
            raise ValueError('unknown near_bed %r; expected one of %r'
                             % (near_bed, sorted(modes)))
        if tau_d < 0.0:
            raise ValueError('tau_d must be >= 0 Pa, got %g' % tau_d)
        adapt = {'none': 0, 'armanini': 1, 'constant': 2, 'carried': 3,
                 'two_layer': 4}
        if adaptation not in adapt:
            raise ValueError('unknown adaptation %r; expected one of %r'
                             % (adaptation, sorted(adapt)))
        if adaptation == 'constant':
            if adaptation_alpha is None or adaptation_alpha <= 0.0:
                raise ValueError("adaptation='constant' needs adaptation_alpha > 0, got %r"
                                 % (adaptation_alpha,))
        elif adaptation_alpha is not None:
            raise ValueError("adaptation_alpha is only used with adaptation='constant'")
        self.sediment_deposition_mode = laws[law]
        self.sediment_tau_d = float(tau_d)
        self.sediment_d_star_mode = modes[near_bed]
        self.sediment_a_h_floor = float(reference_height_floor)
        if layer_fraction is None:
            layer_fraction = TWO_LAYER_FRACTION
        if not 0.0 <= layer_fraction <= 0.5:
            raise ValueError('layer_fraction must be in [0, 0.5], 0 meaning the '
                             'reference height a, got %r' % (layer_fraction,))
        if not exchange_factor > 0.0:
            raise ValueError('exchange_factor must be > 0, got %r' % (exchange_factor,))
        self.sediment_adaptation_mode = adapt[adaptation]
        self.sediment_adaptation_alpha = (float(adaptation_alpha)
                                          if adaptation == 'constant' else 1.0)
        self.sediment_layer_fraction = float(layer_fraction)
        self.sediment_exchange_factor = float(exchange_factor)
        if velocity_profile is None:
            velocity_profile = (adaptation == 'two_layer')
        if velocity_profile and adaptation != 'two_layer':
            raise ValueError("velocity_profile needs adaptation='two_layer'")
        self.sediment_velocity_profile = 1 if velocity_profile else 0
        self.tracer_speed_factor = None
        betas = {'none': 0, 'van_rijn': 1}
        if rouse_beta not in betas:
            raise ValueError("rouse_beta must be 'none' or 'van_rijn', got %r" % (rouse_beta,))
        if not rouse_scale > 0.0:
            raise ValueError('rouse_scale must be > 0, got %r' % (rouse_scale,))
        self.sediment_rouse_beta_mode = betas[rouse_beta]
        self.sediment_rouse_scale = float(rouse_scale)
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_shear_closure(self, closure='quadratic_drag', max_slope=None,
                          freeze_slope=False):
        """Select how bed shear stress is obtained (spec 3.1 / 3.4).

        Parameters
        ----------
        closure : {'quadratic_drag', 'depth_slope', 'energy_slope'}
            `'quadratic_drag'` (default) -- `[T-1]`, `tau_b = rho f_c |v|^2`,
            which makes no equilibrium assumption.
            `'depth_slope'` -- `[T-7]`, `tau_b = rho g h S` with `S` the bed
            slope magnitude, as `aSM16` Eqs 6-7 and hence anugaSed: the
            steady uniform (normal) flow approximation.
            `'energy_slope'` -- `[T-7e]`, the same form with `S` the
            **free-surface** slope, which under the shallow-water assumption
            is the energy grade line, so it drops `[T-7]`'s equilibrium
            assumption. It is what the older `Bed_shear_erosion_operator`
            used (`EN_slope`).
        max_slope : float, optional
            Cap on `S` for the two slope closures, applied per cell each
            step. `None` (default) or 0: no cap. A stated, spatially uniform
            bound in place of anugaSed's undocumented `S <- min(S, mean(S)/2)`.
        freeze_slope : bool, optional
            `'depth_slope'` only. Take `S` from the bed as it stands when
            this is called, or when the first grain size is registered,
            whichever is later, and keep it for the run instead of re-reading
            the evolving bed. See :meth:`bed_slope_magnitude`.

        Notes
        -----
        **Which to use.** Prefer `'quadratic_drag'`, or `'energy_slope'`
        wherever the bed slope is not a good proxy for the energy slope:
        backwater, a pool-riffle sequence, a flat bed that is drawing down, a
        dam break. `'depth_slope'` with `freeze_slope=True` is the
        reproducible way to run the anugaSed closure on a moving bed.

        **How `S` is formed.** The magnitude of the least-squares gradient of
        the bed, or of the stage for `'energy_slope'`, over the cell and its
        neighbours, one-sided at boundaries, so a plane gives its slope in
        every cell, walls and corners included. It is NOT read from the edge
        values, which the DE algorithms rebuild each step through the
        hydrodynamic limiter -- that gave zero slope along reflective walls.
        A slope discontinuity is resolved over one cell either side.

        **The feedback.** With the bed evolving, `'depth_slope'` feeds back
        on itself: erosion roughens the bed, a rougher bed has steeper local
        slopes, steeper slopes erode faster. That is the closure, not the
        discretisation; anugaSed contains it with a domain-global clamp, and
        `max_slope` and `freeze_slope` are the two explicit bounds offered in
        its place.

        **Why `[T-1]` is recommended.** Spec 3.4 keeps `[T-7]` only for
        reproducing published anugaSed results, for three reasons:
        normal-flow equilibrium is exactly what fails in the dam-breach and
        outburst floods this work targets; `S` should be the energy slope,
        not the bed slope, and substituting it into `[T-7]` recovers `[T-1]`
        identically; and the domain-global slope clamp anugaSed applies has
        no counterpart in their own manual.

        **This does not reproduce anugaSed exactly.** Their code additionally
        divides the elevation gradient by a domain-mean cell size and applies
        `S <- min(S, mean(S)/2)`. Neither is in `aSM16`; the first is
        dimensionally inconsistent and the second is the undocumented clamp of
        divergence D1a. `[T-7]` here follows the manual, not the code.
        """
        closures = {'quadratic_drag': 0, 'depth_slope': 1, 'energy_slope': 2}
        if closure not in closures:
            raise ValueError('unknown shear closure %r; expected one of %r'
                             % (closure, sorted(closures)))
        if max_slope is not None and max_slope < 0.0:
            raise ValueError('max_slope must be >= 0, got %g' % max_slope)
        if freeze_slope and closure != 'depth_slope':
            raise ValueError("freeze_slope applies to closure 'depth_slope' "
                             "only: the free surface of 'energy_slope' moves, "
                             "and 'quadratic_drag' uses no slope")
        self.sediment_shear_closure = closures[closure]
        self.sediment_max_slope = float(max_slope) if max_slope else 0.0
        self.sediment_slope_frozen = 1 if freeze_slope else 0
        if freeze_slope:
            self._freeze_bed_slope()
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def bed_slope_magnitude(self, values=None):
        """The per-cell slope magnitude the depth-slope closure uses, `[T-7]`.

        The least-squares gradient of the centroid `values` (default: the
        elevation quantity) over each cell and its neighbours, one-sided at
        boundaries, exactly as the sediment kernel computes it. Returns an
        array of length `number_of_elements`.
        """
        if values is None:
            values = self.quantities['elevation'].centroid_values
        f = num.asarray(values, dtype=float)
        cc = self.centroid_coordinates
        nb = self.neighbours
        n = len(f)
        sxx = num.zeros(n); sxy = num.zeros(n); syy = num.zeros(n)
        sxf = num.zeros(n); syf = num.zeros(n)
        for i in range(3):
            k = num.where(nb[:, i] >= 0)[0]
            j = nb[k, i]
            dx = cc[j, 0] - cc[k, 0]
            dy = cc[j, 1] - cc[k, 1]
            df = f[j] - f[k]
            sxx[k] += dx * dx; sxy[k] += dx * dy; syy[k] += dy * dy
            sxf[k] += dx * df; syf[k] += dy * df
        det = sxx * syy - sxy * sxy
        tr = sxx + syy
        gx = num.zeros(n); gy = num.zeros(n)
        plane = det > 1.0e-12 * tr * tr
        gx[plane] = (syy[plane] * sxf[plane] - sxy[plane] * syf[plane]) / det[plane]
        gy[plane] = (sxx[plane] * syf[plane] - sxy[plane] * sxf[plane]) / det[plane]
        line = ~plane & (tr > 0.0)
        gx[line] = sxf[line] / tr[line]
        gy[line] = syf[line] / tr[line]
        return num.sqrt(gx * gx + gy * gy)

    def _freeze_bed_slope(self):
        """Record the current bed slope for a frozen depth-slope closure.

        Called by set_shear_closure(freeze_slope=True) and again when the
        first grain size allocates the slope array; the later of the two
        wins, so the bed as it stands when the sediment is fully configured
        is what the run uses.
        """
        if not getattr(self, 'sediment_slope_frozen', 0):
            return
        if self.sediment_slope_work is None:
            self.sediment_slope_work = num.zeros(self.number_of_elements,
                                                 dtype=num.float64)
        self.sediment_slope_work[:] = self.bed_slope_magnitude()

    def initialize_sediment_operator(self, porosity=None, c_max=None,
                                     c_pack=None, bed_evolution=None,
                                     rho_w=None, morphological_factor=None,
                                     description=None, label=None,
                                     logging=False, verbose=False):
        """Switch sediment transport on, and return the operator.

        This is the entry point for sediment transport. It takes the
        DOMAIN-WIDE parameters -- the ones that describe the run as a whole
        rather than any one fraction -- and creates the single
        :class:`~anuga.operators.sediment_operator.Sediment_transport_operator`
        that carries the bed exchange.

        Grain sizes are added separately, with :meth:`add_sediment_fraction`:

        .. code-block:: python

            domain.initialize_sediment_operator(porosity=0.28, rho_w=1000.0)
            domain.add_sediment_fraction('sand', diameter=2.0e-4)
            domain.add_sediment_fraction('silt', diameter=2.0e-5)

        The split is the point. A parameter belongs to exactly one of the two
        calls, so there is never a question of which call wins: `porosity` and
        `rho_w` are properties of the run, `diameter` and `rho_s` are properties
        of a fraction.

        ONE OPERATOR PER DOMAIN. Calling this twice returns the same operator,
        applying any parameters given the second time -- the kernel makes a
        single pass over every registered fraction, so a second operator in
        the fractional-step list would apply the bed exchange twice per step.

        Calling it is optional in the simplest case: :meth:`add_sediment_fraction` will
        create the operator with default domain-wide parameters if none exists.
        The two may be called in either order.

        The closure choices -- which shear, erosion, deposition and bedload laws
        to use -- have their own setters, because each carries its own
        parameters and validation: :meth:`set_shear_closure`,
        :meth:`set_sediment_friction`, :meth:`set_bed_material`,
        :meth:`set_deposition` and :meth:`set_bedload`.

        Parameters
        ----------
        porosity : float, optional
            Bed porosity `lambda` in `[G-4]`. Default 0.30.
        c_max : float, optional
            `[L-2]`, the ceiling on depth-averaged concentration. Default 0.30.
        c_pack : float, optional
            `[L-4]`, maximum packing bounding near-bed concentration.
        bed_evolution : bool, optional
            `True` (default) evolves the bed; `False` is the fixed-bed stage.
        rho_w : float, optional
            Water density, kg/m3. Default 1000.
        description, label, logging, verbose
            Passed to the operator; see
            :class:`~anuga.operators.base_operator.Operator`.

        Returns
        -------
        Sediment_transport_operator
            The domain's sediment operator.

        See Also
        --------
        add_sediment_fraction : register one sediment fraction.
        set_sediment_parameters : change the domain-wide parameters later.
        sediment_summary : print the complete active configuration.
        """
        from anuga.operators.sediment_operator import (
            Sediment_transport_operator)

        if (porosity is not None or c_max is not None or c_pack is not None
                or bed_evolution is not None or rho_w is not None
                or morphological_factor is not None):
            self.set_sediment_parameters(porosity=porosity, c_max=c_max,
                                         c_pack=c_pack,
                                         bed_evolution=bed_evolution,
                                         rho_w=rho_w,
                                         morphological_factor=morphological_factor)

        return Sediment_transport_operator(
            self, description=description, label=label, logging=logging,
            verbose=verbose)

    def add_sediment_fraction(self, name, diameter, rho_s=2650.0,
                              tau_c_star=0.04, d_star=1.0, beta=None,
                              initial_concentration=0.0,
                              reference_height=None, **settling_kwargs):
        """Register one suspended sediment fraction and return its index.

        Everything here is a property of THIS fraction. The domain-wide
        parameters live on :meth:`initialize_sediment_operator`; in particular
        there is no `rho_w` here, because there is one fluid.

        .. code-block:: python

            domain.add_sediment_fraction('sand', diameter=2.0e-4)
            domain.add_sediment_fraction('silt', diameter=2.0e-5, tau_c_star=0.11)

        A sediment fraction is a tracer with settling parameters attached, so it
        inherits the transport, boundary and conservation machinery described
        under :ref:`tracers`, and takes the tracer slot of the same index. Do
        not interleave :meth:`add_tracer` and `add_sediment_fraction` on the same
        domain if you rely on that correspondence.

        If the domain has no sediment operator yet, one is created with the
        default domain-wide parameters.

        Parameters
        ----------
        name : str
            Identifier for this fraction, e.g. 'sand'. Also its tracer name.
        diameter : float
            Grain diameter in metres.
        rho_s : float, optional
            Sediment particle density, kg/m3. Default 2650 (quartz). Enters as
            the submerged specific gravity `R = rho_s/rho_w - 1`.
        tau_c_star : float, optional
            Critical Shields stress for entrainment `[E-1]`. Default 0.04.
            Setting it to 0 disables entrainment for this fraction, leaving
            deposition only.
        d_star : float, optional
            Ratio of near-bed to depth-averaged concentration in `[D-1]`.
            Default 1.0, the well-mixed limit. Ignored when the domain's
            `near_bed` mode is `'rouse'`, which computes it per cell.
        beta : float, optional
            Edge reconstruction limiter for this fraction's tracer.
        initial_concentration : float or array-like, optional
            Initial `c_s`; seeds `m = h*c` consistently.
        reference_height : float, optional
            `a` in `[S-4]`, in metres. Defaults to `2*diameter`.
        **settling_kwargs
            Passed to :meth:`settling_velocity` -- e.g. `shape='natural'`.

        Returns
        -------
        int
            The index of this fraction, which is also its tracer index.

        See Also
        --------
        initialize_sediment_operator : the domain-wide parameters.
        """
        # A domain-wide parameter here would otherwise fall into
        # settling_kwargs and surface as a TypeError from settling_velocity,
        # which says nothing about what the caller did wrong.
        domain_wide = sorted(set(settling_kwargs) &
                             {'porosity', 'c_max', 'c_pack', 'bed_evolution',
                              'rho_w'})
        if domain_wide:
            raise TypeError(
                'add_sediment_fraction() got %s, which %s of the run rather than '
                'of one fraction; pass %s to initialize_sediment_operator() '
                'or set_sediment_parameters()'
                % (', '.join(domain_wide),
                   'are properties' if len(domain_wide) > 1 else 'is a property',
                   'them' if len(domain_wide) > 1 else 'it'))

        if not any(isinstance(op, self._sediment_operator_class())
                   for op in getattr(self, 'fractional_step_operators', ())):
            self.initialize_sediment_operator()

        return self._allocate_sediment_fraction(
            name, diameter, rho_s=rho_s, tau_c_star=tau_c_star,
            d_star=d_star, beta=beta,
            initial_concentration=initial_concentration,
            reference_height=reference_height, **settling_kwargs)

    @staticmethod
    def _sediment_operator_class():
        from anuga.operators.sediment_operator import (
            Sediment_transport_operator)
        return Sediment_transport_operator

    def _set_sediment_rho_w(self, rho_w):
        """Set the domain-wide water density and refresh what derives from it.

        `rho_w` is a property of the fluid, so there is one per domain rather
        than one per fraction. Two quantities are computed from it at
        registration time -- the submerged specific gravity `R = rho_s/rho_w - 1`
        and the settling velocity `v_s` `[S-1]` -- so changing it afterwards has
        to recompute them, or already-registered fractions keep values from
        the old density and nothing says so.
        """
        if rho_w <= 0.0:
            raise ValueError('rho_w must be > 0, got %g' % rho_w)
        self.sediment_rho_w = float(rho_w)
        for i in range(self.n_sediment_classes):
            rho_s = self._sediment_rho_s[i]
            self.sediment_R[i] = rho_s / self.sediment_rho_w - 1.0
            self.sediment_settling_velocity[i] = self.settling_velocity(
                self.sediment_diameter[i], rho_s=rho_s,
                rho_w=self.sediment_rho_w,
                **self._sediment_settling_kwargs[i])

    def set_sediment_parameters(self, porosity=None, c_max=None, c_pack=None,
                                bed_evolution=None, rho_w=None,
                                morphological_factor=None):
        """Set the scalar sediment parameters, with validation.

        Everything here is a physical property of the run, not a numerical
        knob. All are optional; only what you pass is changed.

        Parameters
        ----------
        porosity : float
            Bed porosity `lambda` in `[G-4]`. The sediment VOLUME leaving
            suspension is `(1-lambda) dz`, the rest being pore space filled
            from the water column. Default 0.30; LM15 use 0.28.
        c_max : float
            `[L-2]`, the ceiling on depth-averaged volumetric concentration.
            Default 0.30 (FG21); aS16 use 0.20.
        c_pack : float
            `[L-4]`, maximum packing bounding the NEAR-BED concentration
            `c_b = d* c`. Default 0.65, the same constant that bounds `E*` in
            `[E-1]`. Only bites when `d* != 1`.
        bed_evolution : bool
            Spec 2.4's coupling stage. `True` (default) evolves the bed via
            `[G-4]`/`[G-5]`; `False` is the FIXED BED of Phase 3, which is
            RDycore v1.0's configuration and what the analytic
            constant-depth deposition solutions assume.
        rho_w : float
            Water density, used to form the dimensional bed shear stress.
        morphological_factor : float
            Morphological acceleration factor `M` (Delft3D's MORFAC), which
            multiplies every step's bed change so the bed reaches a
            morphological time `M` times longer than the hydrodynamic one
            simulated. Default 1 (off). See the Notes.

        Notes
        -----
        **The morphological factor.** `M` multiplies the bed change from the
        suspended exchange `[G-4]` and from bedload `[G-5]` alike, and the
        erodible base `[L-5]` is still respected. The water column is
        untouched, so the suspension keeps adapting on its own fast time
        scale while the bed runs ahead. That holds while the bed changes
        little over one hydrodynamic adjustment time: a tidal case needs an
        `M` for which `M` tidal cycles average out, while a flume in steady
        flow tolerates 10 or more. The bed and water-column sediment budgets
        then differ by exactly `M`, by construction.
        """
        if morphological_factor is not None:
            if not morphological_factor > 0.0:
                raise ValueError('morphological_factor must be > 0, got %g'
                                 % morphological_factor)
            self.sediment_morphological_factor = float(morphological_factor)
        if porosity is not None:
            if not 0.0 <= porosity < 1.0:
                raise ValueError('porosity must be in [0, 1), got %g' % porosity)
            self.sediment_porosity = float(porosity)
        if c_max is not None:
            if c_max <= 0.0:
                raise ValueError('c_max must be > 0, got %g' % c_max)
            self.sediment_c_max = float(c_max)
        if c_pack is not None:
            if c_pack <= 0.0:
                raise ValueError('c_pack must be > 0, got %g' % c_pack)
            self.sediment_c_pack = float(c_pack)
        if bed_evolution is not None:
            self.sediment_bed_evolution = bool(bed_evolution)
            # May have just been turned on after the classes were registered.
            self._sync_elevation_storage()
        if rho_w is not None:
            self._set_sediment_rho_w(rho_w)
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def sediment_summary(self):
        """Return the complete active sediment configuration as text.

        Every choice that affects the answer, with its units and the spec
        label it implements. Worth printing at the top of any run: the module
        has enough switches that "which erosion law was that?" is a real
        question six months later, and several of the choices are physics
        statements rather than tuning (spec 4.1.1).
        """
        if self.n_sediment_classes == 0:
            return 'sediment: no sediment fractions registered'

        ero = {0: "Shields / Smith-McLean, non-cohesive (sand, gravel), "
                  "gamma0=%.4g   [E-1]" % self.sediment_gamma0,
               1: "Hanson & Simon, cohesive (silt, clay)   [E-3]",
               2: "Partheniades (RDycore)   [E-4]",
               3: "de Leeuw et al. 2020, non-cohesive bed material load "
                  "(A=%.3g alpha=%.3g beta=%.3g, X0=%.3g, k_s=%s)   [E-6]"
                  % (self.sediment_dl_A, self.sediment_dl_alpha,
                     self.sediment_dl_beta, self.sediment_dl_threshold,
                     ('%.4g m' % self.sediment_dl_ks) if self.sediment_dl_ks > 0.0
                     else '%.3g d' % self.sediment_dl_ks_factor),
               }[self.sediment_erosion_mode]
        dep = {0: "D = d* c v_s   [D-1]", 1: "D = v_s c (1 - tau_b/tau_d)   [D-2]"
               }[self.sediment_deposition_mode]
        if self.sediment_deposition_mode == 0 and self.sediment_adaptation_mode in (3, 4):
            # the near-bed concentration is carried, not d* c: say so here
            # rather than leave [D-1] to be read as the whole story
            dep = ("D = c_b v_s with c_b carried   [D-1] under %s"
                   % ('[D-4]' if self.sediment_adaptation_mode == 3 else '[D-5]'))
        dstar = {0: "constant, per fraction", 1: "Rouse profile   [S-4]"
                 }[self.sediment_d_star_mode]
        shear = {0: "quadratic drag, tau_b = rho f_c |v|^2   [T-1]",
                 1: "depth-slope, tau_b = rho g h S (bed slope; aSM16)   [T-7]",
                 2: "energy-slope, tau_b = rho g h S (free surface)   [T-7e]"
                 }[self.sediment_shear_closure]
        if self.sediment_shear_factor is not None:
            f = self.sediment_shear_factor
            shear += ("; amplified by %.3g to %.3g per cell   [T-16]"
                      % (f.min(), f.max()))
        if self.sediment_shear_closure in (1, 2):
            bounds = []
            if self.sediment_max_slope > 0.0:
                bounds.append("S capped at %g" % self.sediment_max_slope)
            if self.sediment_slope_frozen:
                bounds.append("S frozen at the setup bed")
            if bounds:
                shear += "; " + ", ".join(bounds)
        fric = {0: "constant n, from the domain friction quantity",
                1: "larsen_lamb, n = %.5f   [T-13..15]" % self.sediment_manning_ll,
                2: "wilson, bed=%s, D=%.4g m   [T-8..10]"
                   % (['sand', 'gravel', 'boulder'][self.sediment_wilson_bed],
                      self.sediment_wilson_D)}[self.sediment_friction_mode]
        bl = ("off" if self.sediment_bedload_mode == 0 else
              ("Engelund-Hansen, TOTAL LOAD (suspended source disabled)   [K-5]"
               if self.sediment_bedload_mode == 2 else
               "power law, K=%.4g m=%.4g tau_c*=%.4g   [K-1]"
               % (self.sediment_bedload_K, self.sediment_bedload_m,
                  self.sediment_bedload_tau_c_star)))

        L = ['sediment configuration',
             # The name is free text ('sand', 'silt'), so carry the
             # diameter alongside it -- the name alone does not say what the
             # fraction actually is.
             '  fractions          : %d  --  %s'
             % (self.n_sediment_classes,
                ', '.join('%s (d=%.4g m)' % (nm, self.sediment_diameter[i])
                          for i, nm in
                          enumerate(self.get_sediment_names()))),
             '  erosion            : %s' % ero,
             '  deposition         : %s' % dep,
             '  near-bed d*        : %s' % dstar,
             '  shear closure      : %s' % shear,
             '  friction closure   : %s' % fric,
             '  bedload            : %s' % bl,
             *(['  bedload h_min [K-7]: %.4g m (none below, full above twice it)'
                % self.sediment_bedload_h_min]
               if self.sediment_bedload_mode and self.sediment_bedload_h_min > 0.0 else []),
             '  bed evolution      : %s  (%s)'
             % (self.sediment_bed_evolution,
                'Phase 4, evolving' if self.sediment_bed_evolution
                else 'Phase 3, FIXED bed'),
             '  porosity lambda    : %.4g' % self.sediment_porosity,
             '  morphological M    : %.4g%s' % (self.sediment_morphological_factor,
                                               '' if self.sediment_morphological_factor == 1.0
                                               else '  (bed change per step x M)'),
             '  c_max      [L-2]   : %.4g' % self.sediment_c_max,
             '  c_pack     [L-4]   : %.4g' % self.sediment_c_pack,
             '  rho_w              : %.4g kg/m3' % self.sediment_rho_w]
        if self.sediment_erosion_mode in (1, 2):
            L.append('  tau_crit           : %.4g Pa' % self.sediment_tau_crit)
        if self.sediment_erosion_mode == 1:
            L.append('  K_e        [E-5]   : %.4e m3/N/s' % self.sediment_K_e)
        if self.sediment_erosion_mode == 2:
            L.append('  K_p        [E-4]   : %.4e kg/m2/s'
                     % self.sediment_K_partheniades)
        if self.sediment_deposition_mode == 1:
            L.append('  tau_d      [D-2]   : %.4g Pa' % self.sediment_tau_d)
        if self.sediment_adaptation_mode == 1:
            L.append('  adaptation [D-3]   : Armanini & Di Silvio alpha(w_s/u*, a/h); '
                     'E - D = alpha v_s (c_eq - c)')
        elif self.sediment_adaptation_mode == 2:
            L.append('  adaptation [D-3]   : constant alpha = %.4g; E - D = alpha v_s (c_eq - c)'
                     % self.sediment_adaptation_alpha)
        elif self.sediment_adaptation_mode == 3:
            L.append('  adaptation [D-4]   : carried near-bed ratio r_b per fraction, settling lag '
                     'T = (z_c - a)/w_s when d* rises%s'
                     % ('' if self.sediment_nearbed_base >= 0 else ' (tracers registered at evolve)'))
        elif self.sediment_adaptation_mode == 4:
            L.append('  adaptation [D-5]   : two-layer suspension, near-bed layer %s, exchange '
                     'set to the Rouse equilibrium, rate x %.3g%s'
                     % ('%.3g h' % self.sediment_layer_fraction if self.sediment_layer_fraction > 0
                        else 'a thick', self.sediment_exchange_factor,
                        '' if self.sediment_nearbed_base >= 0 else ' (tracers registered at evolve)'))
            if self.sediment_velocity_profile:
                L.append('  velocity profile   : layers advected at their log-law mean speeds')
        if self.sediment_d_star_mode == 1:
            L.append('  a/h floor          : %.4g' % self.sediment_a_h_floor)
            if self.sediment_rouse_beta_mode or self.sediment_rouse_scale != 1.0:
                L.append('  Rouse number       : %sx %.3g'
                         % ('Z / (1 + 2 (w_s/u*)^2), ' if self.sediment_rouse_beta_mode else 'Z ',
                            self.sediment_rouse_scale))
        mask = self._sediment_erodible_mask
        if self._sediment_user_base is not None:
            t = self.erodible_thickness()
            # Report only where erosion is actually permitted; locked cells
            # carry zero thickness and would drag the minimum to 0 whatever
            # the layer is.
            tt = t if mask is None else t[mask]
            if len(tt):
                L.append('  erodible base [L-5]: thickness %.4g to %.4g m, '
                         '%d of %d erodible cells at bedrock'
                         % (tt.min(), tt.max(), int((tt <= 0.0).sum()), len(tt)))
            else:
                L.append('  erodible base [L-5]: set, but no cell is erodible')
        else:
            L.append('  erodible base [L-5]: none (unlimited depth)')
        if self.sediment_bed_composition:
            Fa = self.get_bed_composition('active')
            L.append('  bed composition    : active layer %.3g m; mean active '
                     'fractions %s' % (self.sediment_active_layer,
                     ', '.join('%s %.3g' % (k, float(num.mean(v))) for k, v in Fa.items())))
        if self.sediment_repose_tan > 0.0:
            L.append('  angle of repose    : %.1f degrees, relax %.2g, '
                     'max %d sweeps'
                     % (num.degrees(num.arctan(self.sediment_repose_tan)),
                        self.sediment_repose_relax,
                        self.sediment_repose_max_sweeps))
        else:
            L.append('  angle of repose    : off')
        if mask is None:
            L.append('  erodible region    : whole domain')
        else:
            L.append('  erodible region    : %d of %d cells erodible '
                     '(%d locked at their current bed)'
                     % (int(mask.sum()), len(mask), int((~mask).sum())))
        L.append('  ([E-1] and the like are cross-references to the term in '
                 'the physics;')
        L.append('   see the Sediment physics appendix -- the description '
                 'before each')
        L.append('   label is the whole story.)')
        L.append('  per fraction:')
        for i, nm in enumerate(self.get_sediment_names()):
            L.append('    %-10s d=%.4g m  v_s=%.4e m/s  R=%.4g  tau_c*=%.4g'
                     % (nm, self.sediment_diameter[i],
                        self.sediment_settling_velocity[i],
                        self.sediment_R[i], self.sediment_tau_c_star[i]))
        return '\n'.join(L)

    # [E-6] constant sets (A, alpha, beta, threshold).
    DE_LEEUW_FITS = {
        # de Leeuw et al. (2020), Eq 26a, fitted to sand and gravel with the
        # skin-friction shear velocity.
        'de_leeuw_2020': (4.74e-4, 1.5, 1.18, 0.015),
        # As implemented in the Delta-X Wax Lake Delta sediment model (Wang,
        # Salter & Lamb 2023, ORNL DAAC doi:10.3334/ORNLDAAC/2309), which
        # attributes it to Nghiem et al. (2022) and uses it for flocculated
        # mud (as bed material load) and sand alike.
        'nghiem_2022': (7.04e-4, 0.94475138, 1.81, 0.015),
    }

    def set_bed_material(self, material='noncohesive', tau_crit=0.088,
                         K_e=None, entrainment='smith_mclean',
                         de_leeuw_fit='de_leeuw_2020', A=None, alpha=None,
                         beta=None, threshold=None, skin_roughness=None,
                         skin_roughness_factor=3.0):
        """Select the erosion law by naming the BED MATERIAL (spec 4.1.1).

        `'noncohesive'` (default) -- sand, gravel, boulders, and mud carried
        as flocculated bed material load. The entrainment relation is chosen
        by `entrainment`:

        * `'smith_mclean'` (default) -- Shields entrainment via Smith & McLean
          / Parker, `[E-1]`/`[E-2]`, with a critical Shields stress per
          fraction (`tau_c_star` on `add_sediment_fraction`).
        * `'de_leeuw'` -- de Leeuw et al. (2020), `[E-6]`:
          ``E* = A X^beta / (1 + 3 A X^beta)`` with
          ``X = (u*_skin / v_s)^alpha Fr - threshold``, the near-bed
          concentration at 0.1 h (pair it with the Rouse near-bed profile
          and a 0.1 reference-height floor to match). `u*_skin` comes from
          Manning-Strickler on a skin roughness `k_s`; the per-fraction
          `tau_c_star` is not used. `de_leeuw_fit` picks the constants:
          `'de_leeuw_2020'` (their Eq 26a, sand and gravel) or
          `'nghiem_2022'` (as the Delta-X Wax Lake model uses them, for
          flocculated mud and sand); `A`, `alpha`, `beta` and `threshold`
          override them one by one. `skin_roughness` fixes `k_s` in metres
          for every fraction (de Leeuw use 3 D84 of the bed); by default it
          is `skin_roughness_factor` (3) times each fraction's diameter.

        `'partheniades'` -- `[E-4]`, `E = K_p (tau_b - tau_c)/tau_c`, the form
        RDycore-sediment uses. `K_e` here is the Partheniades coefficient as a
        **mass** flux in kg m-2 s-1 (RDy26 use 1e-4), a different quantity from
        `[E-5]`'s erodibility, and is divided internally by the class density.

        `'cohesive'` -- silt, clay, cohesive bank material. Excess DIMENSIONAL
        shear via Hanson & Simon `[E-3]`, `E = K_e (tau_b - tau_c)`, with the
        jet-test erodibility `[E-5]` `K_e = 0.2e-6 / sqrt(tau_c)` unless `K_e`
        is given explicitly. `tau_crit` is in **pascals**, not Shields units;
        aS16 use 0.088.

        Notes
        -----
        `[E-1]` and `[E-3]` are **not competing formulations of the same
        physics** -- they describe different sediment, calibrated from different
        experiments. Spec 4.1.1 puts it plainly: choosing between them is a
        statement about the bed, and getting it wrong is a physics error, not a
        tuning error. That is why this method is named for the material rather
        than for the equation.
        """
        materials = {'noncohesive': 0, 'cohesive': 1, 'partheniades': 2}
        if material not in materials:
            raise ValueError('unknown bed material %r; expected one of %r'
                             % (material, sorted(materials)))
        if tau_crit <= 0.0:
            raise ValueError('tau_crit must be > 0 Pa, got %g' % tau_crit)
        if entrainment not in ('smith_mclean', 'de_leeuw'):
            raise ValueError("entrainment must be 'smith_mclean' or 'de_leeuw', got %r"
                             % (entrainment,))
        if entrainment == 'de_leeuw' and material != 'noncohesive':
            raise ValueError("entrainment='de_leeuw' is a non-cohesive relation; "
                             "use it with material='noncohesive'")
        if entrainment == 'de_leeuw':
            if de_leeuw_fit not in self.DE_LEEUW_FITS:
                raise ValueError('unknown de_leeuw_fit %r; expected one of %r'
                                 % (de_leeuw_fit, sorted(self.DE_LEEUW_FITS)))
            fit = dict(zip(('A', 'alpha', 'beta', 'threshold'),
                           self.DE_LEEUW_FITS[de_leeuw_fit]))
            for key, value in (('A', A), ('alpha', alpha), ('beta', beta),
                               ('threshold', threshold)):
                if value is not None:
                    fit[key] = float(value)
            if fit['A'] <= 0.0 or fit['alpha'] <= 0.0 or fit['beta'] <= 0.0:
                raise ValueError('de Leeuw A, alpha and beta must be > 0, got %r' % fit)
            (self.sediment_dl_A, self.sediment_dl_alpha, self.sediment_dl_beta,
             self.sediment_dl_threshold) = (fit['A'], fit['alpha'], fit['beta'],
                                            fit['threshold'])
            if skin_roughness is not None and skin_roughness <= 0.0:
                raise ValueError('skin_roughness must be > 0 m, got %g' % skin_roughness)
            self.sediment_dl_ks = float(skin_roughness) if skin_roughness else 0.0
            self.sediment_dl_ks_factor = float(skin_roughness_factor)

        self.sediment_erosion_mode = 3 if entrainment == 'de_leeuw' else materials[material]
        self.sediment_tau_crit = float(tau_crit)
        if material == 'partheniades' and K_e is not None:
            # [E-4]'s coefficient is a MASS flux in kg m-2 s-1, not [E-5]'s
            # m3 N-1 s-1. Different quantity, so it lands in its own field.
            self.sediment_K_partheniades = float(K_e)
        # [E-5] Hanson & Simon jet-test erodibility, k_d = 0.2 tau_c^-0.5 in
        # cm3 N-1 s-1; the 1e-6 converts to m3 N-1 s-1.
        self.sediment_K_e = (float(K_e) if K_e is not None
                             else 0.2e-6 / self.sediment_tau_crit**0.5)

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_sediment_friction(self, mode='constant', k_s=None, r_d=2.0,
                              r_br=2.0, sigma_br=None, bed='sand',
                              grain_size=None):
        """Select the friction closure used by the sediment kernel (spec 3.3).

        This affects only `tau_b` in the sediment source term; the
        hydrodynamic friction operator is untouched.

        `mode='constant'` (default)
            `f_c = g n^2 h^(-1/3)` `[T-6]` with `n` from the domain's own
            friction quantity.

        `mode='larsen_lamb'`
            Manning-Strickler `[T-13]`-`[T-15]`: `k_s = r_d r_br sigma_br` and
            `n = k_s^(1/6) / (8.1 sqrt(g))`, uniform in space and time. Pass
            either `k_s` directly or `sigma_br` (one standard deviation of
            bedrock elevation). LL16 measured `sigma_br ~ 5 m` at Moses Coulee,
            giving `k_s = 20 m` and `n = 0.065` -- a **site-measured** quantity,
            not a universal default.

        `mode='wilson'`
            `f_c` from bed type and relative submergence `[T-8]`-`[T-10]`.
            `bed` is 'sand' (uses D50), 'gravel' or 'boulder' (both use D84);
            pass the percentile as `grain_size`.

        Notes
        -----
        **`larsen_lamb` and `wilson` are megaflood/planetary parameterisations,
        not general-purpose flood closures.** W04 is a *Mars outflow channel*
        study under Martian gravity; LL16 is the Channeled Scablands. Neither is
        calibrated for ordinary river or urban flood modelling, which is ANUGA's
        main use. For standard flood work keep `constant` with an `n` from
        conventional tables.
        """
        import math
        from anuga.config import g as _g

        modes = {'constant': 0, 'larsen_lamb': 1, 'wilson': 2}
        if mode not in modes:
            raise ValueError('unknown friction mode %r; expected one of %r'
                             % (mode, sorted(modes)))
        self.sediment_friction_mode = modes[mode]

        if mode == 'larsen_lamb':
            if k_s is None:
                if sigma_br is None:
                    raise ValueError(
                        "larsen_lamb needs either k_s or sigma_br; sigma_br is "
                        "site-measured (LL16 report ~5 m at Moses Coulee) and "
                        "has no universal default")
                k_s = r_d * r_br * sigma_br            # [T-15]
            if k_s <= 0.0:
                raise ValueError('k_s must be > 0, got %g' % k_s)
            self.sediment_manning_ll = k_s**(1.0 / 6.0) / (8.1 * math.sqrt(_g))

        elif mode == 'wilson':
            beds = {'sand': 0, 'gravel': 1, 'boulder': 2}
            if bed not in beds:
                raise ValueError('unknown bed type %r; expected one of %r'
                                 % (bed, sorted(beds)))
            if grain_size is None or grain_size <= 0.0:
                raise ValueError(
                    "wilson needs grain_size > 0 (D50 for sand, D84 for "
                    "gravel/boulder)")

            # `bed` picks the curve; `grain_size` sets the relative submergence
            # h/D it is evaluated at. They are independent, so a boulder curve
            # at a sand-sized D is a well-formed calculation of nothing real --
            # and it fails quietly and hugely. The gravel and boulder relations
            # are logarithmic in h/D, so shrinking D inflates the submergence
            # and drives f_c down: bed='boulder' at D = 2e-4 gives f_c = 0.0016
            # against 0.031 at a plausible D = 0.5, a factor of 19 in f_c and
            # therefore in tau_b. Measured end to end, that under-predicted
            # scour four-fold.
            #
            # Bounds are deliberately loose -- roughly Wentworth, widened so
            # the usual borderline choices pass. This warns, it does not
            # refuse: an unusual bed is the user's call to make.
            plausible = {'sand': (6.0e-5, 2.0e-3),
                         'gravel': (2.0e-3, 0.25),
                         'boulder': (0.05, 10.0)}
            lo, hi = plausible[bed]
            if not lo <= grain_size <= hi:
                warnings.warn(
                    "set_sediment_friction('wilson', bed=%r, grain_size=%g): "
                    "grain_size looks implausible for a %s bed, which the "
                    "relations expect to be roughly %g to %g m. `bed` only "
                    "selects the curve; `grain_size` sets the relative "
                    "submergence h/D it is evaluated at, and nothing "
                    "cross-checks them, so a mismatch runs without error and "
                    "can shift f_c -- and therefore tau_b -- by an order of "
                    "magnitude. Note grain_size is the roughness length scale "
                    "of the BED (D50 for sand, D84 for gravel/boulder), not "
                    "the diameter passed to add_sediment_fraction()."
                    % (bed, grain_size, bed, lo, hi),
                    UserWarning, stacklevel=2)

            self.sediment_wilson_bed = beds[bed]
            self.sediment_wilson_D = float(grain_size)

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    #: Bedload parameter sets for `[K-1]`, keyed by name.
    #: Which of the two Wong & Parker relations FG21 used is still open, so
    #: both are provided and neither is privileged beyond the default.
    BEDLOAD_PARAMETER_SETS = {
        'wong_parker_eq24': (3.97, 1.5, 0.0495),
        'wong_parker_eq23': (4.93, 1.60, 0.0470),
    }

    def set_bedload(self, formula='wong_parker_eq24', K=None, m=None,
                    tau_c_star=None, open_boundaries=None, supply=None,
                    min_depth=0.0):
        """Enable bedload transport `[K-1]`-`[K-4]` and its bed evolution `[G-5]`.

        Parameters
        ----------
        formula : str
            `'wong_parker_eq24'` (default; `K`=3.97, `m`=1.5, `tau_c*`=0.0495),
            `'wong_parker_eq23'` (4.93, 1.60, 0.0470), `'engelund_hansen'`,
            `'grass'`, or `'off'`.
        K, m, tau_c_star : float, optional
            Override the chosen set's values individually. For `'grass'`, `K`
            is the Grass coefficient `A_g` (s^2/m, no default: it is a
            calibration) and `m` the velocity exponent (default 3).
        open_boundaries : sequence of str, optional
            Boundary tags through which bedload passes. Across an open edge
            the flux is the cell's own `q_b . n` (zero gradient), so an
            outflow carries bedload away at the rate it arrives and an inflow
            supplies it at the rate the first cell carries it off: the
            equilibrium-supply condition. Every other boundary edge is closed
            (no bedload across it), which is right for walls and is what
            keeps a closed domain exactly conservative. Give this for the
            inflow and outflow of a reach, or the inflow cell exports and
            never imports and digs a hole that travels downstream. On a
            distributed sub-domain a tag this rank owns no part of is
            ignored; in serial an unknown tag is an error.
        min_depth : float, optional
            `[K-7]` depth (m) below which bedload is ramped to zero: none
            below `min_depth`, full above twice it. Default 0, no ramp.
            Bedload relations assume a flow many grains deep; at a wetting
            front over an erodible bed the thin film cells carry a large
            nominal shear and, unramped, dig steps that collapse the time
            step. A few grain diameters is a physical choice.
        supply : dict, optional
            `{tag: q_b}`: a PRESCRIBED bedload inflow across the edges of
            `tag`, in m2/s volumetric per unit width (mass rate / rho_s),
            in place of the zero-gradient import. The tag is opened if it
            is not already. Use it where the supply is known, as in a flume
            fed at a set rate: the zero-gradient import equals the inflow
            cell's own export, so a cell that aggrades under a fixed inflow
            stage sees its transport and hence its import rise, a feedback
            that runs away in a few hours; a prescribed supply does not
            depend on the cell at all. An outflow edge with a supply set
            still carries the supply in, so give it only to inflow tags.

        Notes
        -----
        **Grass is a TOTAL LOAD relation** `[K-6]`, `q_b = A_g |u|^m`, with no
        threshold and no grain size: the classic law of the Exner test cases
        (Hudson & Sweby 2003), for which the migrating-hump validation has a
        closed form. Like Engelund-Hansen it turns the suspended exchange
        off. A fraction must still be registered for the bedload kernel to
        run.

        **Engelund & Hansen is a TOTAL LOAD relation** `[K-5]`: suspension is
        already inside it, so running the suspended source alongside it double
        counts. Spec 6 calls this out as a critical usage rule, and this method
        enforces it rather than warning -- selecting `'engelund_hansen'` turns
        the suspended exchange off, and it has no threshold by construction.

        Which Wong & Parker relation FG21 used, their Eq 23 or Eq 24, is an open
        item. The two differ enough to matter, so the parameters are exposed:
        resolving it is a change of default, not an edit.
        """
        if min_depth < 0.0:
            raise ValueError('min_depth must be >= 0, got %r' % (min_depth,))
        self.sediment_bedload_h_min = float(min_depth)
        if formula == 'off':
            self.sediment_bedload_mode = 0
        elif formula == 'engelund_hansen':
            self.sediment_bedload_mode = 2
            # [K-5] has no threshold; subtracting one would be a different
            # model. And it already contains suspension:
            self.sediment_bedload_tau_c_star = 0.0
            self._sediment_suspended_enabled = False
        elif formula == 'grass':
            if K is None:
                raise ValueError(
                    "set_bedload('grass') needs K, the Grass coefficient A_g "
                    "(s^2/m): it is a calibration and has no default")
            self.sediment_bedload_mode = 3
            self.sediment_bedload_m = 3.0
            # [K-6] has no threshold, and is total load:
            self.sediment_bedload_tau_c_star = 0.0
            self._sediment_suspended_enabled = False
        elif formula in self.BEDLOAD_PARAMETER_SETS:
            self.sediment_bedload_mode = 1
            (self.sediment_bedload_K, self.sediment_bedload_m,
             self.sediment_bedload_tau_c_star) = \
                self.BEDLOAD_PARAMETER_SETS[formula]
        else:
            raise ValueError(
                'unknown bedload formula %r; expected one of %r'
                % (formula, sorted(self.BEDLOAD_PARAMETER_SETS) +
                   ['engelund_hansen', 'grass', 'off']))

        if formula != 'engelund_hansen':
            if K is not None:
                self.sediment_bedload_K = float(K)
            if m is not None:
                self.sediment_bedload_m = float(m)
            if tau_c_star is not None and formula != 'grass':
                self.sediment_bedload_tau_c_star = float(tau_c_star)

        if open_boundaries is not None:
            if isinstance(open_boundaries, str):
                open_boundaries = (open_boundaries,)
            self._sediment_bedload_open_tags = tuple(open_boundaries)
        if supply is not None:
            for tag, q in supply.items():
                if q < 0.0:
                    raise ValueError('bedload supply for %r must be >= 0 m2/s, got %g' % (tag, q))
                if tag not in self.tag_boundary_cells and not self._is_subdomain():
                    raise ValueError(
                        'no boundary tagged %r on this domain; known tags: %s'
                        % (tag, sorted(self.tag_boundary_cells)))
            self._sediment_bedload_supply = {tag: float(q) for tag, q in supply.items()}
            extra = tuple(t for t in self._sediment_bedload_supply
                          if t not in self._sediment_bedload_open_tags)
            self._sediment_bedload_open_tags = tuple(self._sediment_bedload_open_tags) + extra
        self._build_bedload_open()

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    # ---- [D-4] carried near-bed concentration -----------------------------
    ROUSE_DSTAR_COEFFS = (
        (+1.097192252266e-03, +9.816426876103e-04, +2.816550608693e-04, +2.216981593577e-05),
        (+8.152552738643e-01, +2.984288438662e-01, +4.717126357513e-02, +2.488390592718e-03),
        (-3.858022145865e-02, +7.497943739332e-01, +1.494530687071e-01, +1.016829763599e-02),
        (-1.416163484237e-01, -6.145585548869e-01, -2.181478641118e-01, -1.989085090511e-02),
        (+2.441798567588e-02, +2.477861105262e-01, +1.172562489413e-01, +1.380260943049e-02),
        (+1.714535144604e-02, -4.453006825886e-02, -2.922351673152e-02, -4.300807416335e-03),
        (-4.557991043496e-03, +2.511784955218e-03, +2.810236330103e-03, +5.048472261972e-04))
    ROUSE_ZC_COEFFS = (
        (-1.236620662919e+00, -2.966806356338e-01, -5.451826597743e-02, -3.340785646720e-03),
        (-4.934470965259e-01, +1.986068502121e-01, +2.681666552524e-02, +1.430964808434e-03),
        (+5.693206124955e-01, +4.438916973539e-01, +7.169076944620e-02, +3.074161648449e-03),
        (-8.317934945723e-01, -8.674254433522e-01, -2.464308678880e-01, -1.711508840976e-02),
        (+5.466837296437e-01, +6.103587559357e-01, +2.038814693879e-01, +1.776942243679e-02),
        (-1.671607038599e-01, -1.948379068715e-01, -7.092485415218e-02, -7.151638532153e-03),
        (+1.963501463003e-02, +2.364362131184e-02, +9.087924042213e-03, +1.006461512046e-03))
    ROUSE_LIMITS = (0.01, 2.5, 1.0e-3, 0.15)

    @classmethod
    def _rouse_poly(cls, coeffs, Z, a_h):
        import math
        zlo, zhi, alo, ahi = cls.ROUSE_LIMITS
        Z = min(max(float(Z), zlo), zhi)
        a_h = min(max(float(a_h), alo), ahi)
        L = math.log(a_h)
        P = 0.0
        for row in reversed(coeffs):
            P = P * Z + (row[0] + row[1] * L + row[2] * L * L + row[3] * L ** 3)
        return P, L

    @classmethod
    def rouse_d_star(cls, Z, a_h):
        """The kernel's fitted near-bed ratio `d*(Z, a/h)` of `[S-4]`, for
        scripts that need it on the host (an inflow near-bed concentration,
        say). Clamped to the fitted range like the kernel."""
        import math
        P, L = cls._rouse_poly(cls.ROUSE_DSTAR_COEFFS, Z, a_h)
        return max(math.exp(-Z * L + P), 1.0)

    @classmethod
    def rouse_centroid(cls, Z, a_h):
        """`(z_c - a)/h`, the height of the Rouse profile's centroid above the
        reference level, from the kernel's fit; the settling time of `[D-4]`
        is this times `h / w_s`."""
        import math
        P, _ = cls._rouse_poly(cls.ROUSE_ZC_COEFFS, Z, a_h)
        return math.exp(P)

    def _ensure_nearbed_tracers(self):
        """Register the `[D-4]` near-bed tracers, one per fraction, once every
        fraction exists. Called at the top of evolve; a no-op unless
        `adaptation='carried'` and they are not there yet."""
        if (self.sediment_adaptation_mode not in (3, 4) or self.sediment_nearbed_base >= 0
                or self.n_sediment_classes == 0):
            self._ensure_tracer_speed_factor()
            return
        ncl = self.n_sediment_classes
        # The near-bed tracers go at the END of the tracer list, contiguously,
        # and the kernel reads them as nb_base + s. Any OTHER tracers (a
        # passive salt, say) sit between the fractions and them and are
        # untouched by the sediment kernels, so they are allowed; what is not
        # allowed is adding a tracer after these, which add_tracer refuses.
        base = self.number_of_tracers
        names = list(self.get_sediment_names())
        two_layer = self.sediment_adaptation_mode == 4
        for s in range(ncl):
            # [D-4]: 0 = not set, the kernel takes the local d*; [D-5]: a
            # negative upper-layer mass = not set, the kernel takes the
            # equilibrium partition. Either way the initial state and every
            # inflow start at equilibrium with nothing to set.
            if two_layer:
                self.add_tracer(names[s] + '_upper', initial_value=-1.0)
            else:
                self.add_tracer(names[s] + '_nearbed_ratio', initial_value=0.0)
        self.sediment_nearbed_base = base
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized
        self._ensure_tracer_speed_factor()

    def _ensure_tracer_speed_factor(self):
        """Allocate the `[D-5v]` per-tracer speed factors (all 1: the kernel
        fills them from the first step) once the tracers are final.

        Only where there is something to split. [D-5] is the default closure,
        so a domain with no sediment at all reaches here; allocating for it
        would put a factor on the advection of its passive tracers, which
        mode 2 then reads from a pointer it never mapped."""
        if (not self.sediment_velocity_profile
                or self.sediment_adaptation_mode != 4
                or self.sediment_nearbed_base < 0):
            return
        shape = (self.number_of_tracers, self.number_of_elements)
        if self.tracer_speed_factor is not None and self.tracer_speed_factor.shape == shape:
            return
        self.tracer_speed_factor = num.ones(shape, dtype=num.float64)
        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def _build_bedload_open(self):
        """(Re)build the per-boundary-edge bedload flag from the open tags.

        Allocated with the other bedload scratch so the C bindings always
        find it once a class exists; rebuilt whenever set_bedload changes the
        tags. Zero everywhere means every boundary edge is closed.
        """
        if self.sediment_bedload_open is None:
            self.sediment_bedload_open = num.zeros(self.boundary_length,
                                                   dtype=num.int64)
            self.sediment_bedload_supply = num.full(self.boundary_length, -1.0)
        flags = self.sediment_bedload_open
        flags[:] = 0
        supply = self.sediment_bedload_supply
        supply[:] = -1.0
        for tag in self._sediment_bedload_open_tags:
            if tag not in self.tag_boundary_cells:
                if self._is_subdomain():
                    continue      # this rank owns no part of the tag
                raise ValueError(
                    'no boundary tagged %r on this domain; known tags: %s'
                    % (tag, sorted(self.tag_boundary_cells)))
            edges = num.asarray(self.tag_boundary_cells[tag], dtype=num.intp)
            flags[edges] = 1
            if tag in self._sediment_bedload_supply:
                supply[edges] = self._sediment_bedload_supply[tag]

    def get_sediment_names(self):
        """Return the registered sediment class names, in index order."""
        return list(self._sediment_names)
    def set_tracer_boundary(self, name, tag, value):
        """Concentration that tracer `name` brings in across boundary `tag`.

        Only used where water FLOWS IN. The flux kernel picks the upwind value
        edge by edge from the sign of the water flux:

            inflow  (n.U into the domain)   ->  this boundary value
            outflow (n.U out of the domain) ->  the interior edge value

        so there is nothing to set for an outflow, and nothing that can be set:
        prescribing a concentration on an outflow would over-determine the
        advection, and the kernel would ignore it anyway. A boundary that only
        ever lets water out therefore needs no call at all, and one that
        alternates gets Dirichlet-on-inflow / transmissive-on-outflow
        automatically -- the characteristic condition, without a switch.

        Parameters
        ----------
        name : str
            A registered tracer.
        tag : str
            A boundary tag, as used by `set_boundary`. On a distributed
            sub-domain a tag this rank owns no part of is silently ignored,
            so the same call works on every rank; in serial an unknown tag
            is an error.
        value : float, array-like, or callable
            A scalar applies to every edge carrying `tag`. An array must have
            one value per edge of that tag, ordered as
            `domain.tag_boundary_cells[tag]`. A callable is evaluated as
            `value(t)` at each timestep and must return a scalar or such an
            array, for a time-varying inflow.

        Notes
        -----
        UNSET BOUNDARIES BRING c = 0. The array is zero-filled at
        `add_tracer`, so with no call here an inflow carries clean water. That
        is a modelling assumption, not a neutral default: for salinity or
        suspended sediment it is usually wrong, and it is invisible in the
        output. It is preserved because changing it would silently alter
        existing results, but set it explicitly wherever water enters.
        """
        s = self.get_tracer_index(name)

        if tag not in self.tag_boundary_cells:
            if self._is_subdomain():
                # After distribute() a tag lives only on the ranks that own a
                # piece of it, so 'not here' is the normal case, not an error.
                # Same convention as set_boundary, which ignores a tag the
                # domain does not use. In SERIAL there is no such excuse, so a
                # missing tag is still reported -- it can only be a typo.
                return
            raise ValueError(
                'no boundary tagged %r on this domain; known tags: %s'
                % (tag, sorted(self.tag_boundary_cells)))

        idx = num.asarray(self.tag_boundary_cells[tag], dtype=num.intp)

        if callable(value):
            # Re-evaluated each timestep by update_tracer_boundary_values().
            self._tracer_boundary_functions[(s, tag)] = value
            self._apply_tracer_boundary(s, idx, value(self.get_time()), tag)
        else:
            self._tracer_boundary_functions.pop((s, tag), None)
            self._apply_tracer_boundary(s, idx, value, tag)
        self._push_tracer_boundary_to_device()

    def _push_tracer_boundary_to_device(self):
        """Mirror the host tracer boundary array to the device, if there is one.

        The device copy is made once when the arrays are mapped, and the
        per-step boundary push carries only the hydrodynamic values, so every
        host write to tracer_boundary_values has to be followed by this or the
        flux kernel keeps reading the mapped copy. Same reasoning as the
        external source in set_tracer_source. A no-op off the GPU path.
        """
        if (self.multiprocessor_mode == MULTIPROCESSOR_GPU
                and self.gpu_interface is not None):
            from anuga.shallow_water.sw_domain_gpu_ext import (
                sync_tracer_boundary_to_device)
            sync_tracer_boundary_to_device(self.gpu_interface.gpu_dom)

    def _apply_tracer_boundary(self, s, idx, value, tag):
        """Write one tracer's boundary values for the edges of one tag."""
        v = num.asarray(value, dtype=num.float64)
        if v.ndim == 0:
            self.tracer_boundary_values[s, idx] = float(v)
        elif v.shape == idx.shape:
            self.tracer_boundary_values[s, idx] = v
        else:
            raise ValueError(
                'tracer boundary %r on tag %r: expected a scalar or %d values '
                '(one per edge of that tag), got shape %r'
                % (self._tracer_names[s], tag, idx.size, v.shape))

    def update_tracer_boundary_values(self):
        """Re-evaluate any callable tracer boundary values for this time.

        Called from the evolve loop alongside update_boundary(). A no-op unless
        set_tracer_boundary() was given a callable, so a domain with constant
        boundary concentrations pays nothing.
        """
        if not self._tracer_boundary_functions:
            return
        t = self.get_time()
        for (s, tag), f in self._tracer_boundary_functions.items():
            idx = num.asarray(self.tag_boundary_cells[tag], dtype=num.intp)
            self._apply_tracer_boundary(s, idx, f(t), tag)
        self._push_tracer_boundary_to_device()

    def get_tracer_boundary(self, name, tag):
        """Return the boundary concentrations of `name` on `tag`."""
        s = self.get_tracer_index(name)
        if tag not in self.tag_boundary_cells:
            if self._is_subdomain():
                # This rank owns none of that boundary; it has no values for
                # it, which is not the same thing as the tag being wrong.
                return num.zeros(0, dtype=num.float64)
            raise ValueError('no boundary tagged %r on this domain' % tag)
        idx = num.asarray(self.tag_boundary_cells[tag], dtype=num.intp)
        return self.tracer_boundary_values[s, idx]

    def _is_subdomain(self):
        """True if this domain is one rank's piece of a distributed domain."""
        return getattr(self, 'numproc', 1) > 1

    def _local_tracer_mass(self, s):
        """This rank's share of tracer `s`: owned cells only, NO reduction.

        Deliberately collective-free. The conservation baseline is captured
        during setup, inside add_tracer and set_tracer -- and on a parallel run
        that happens on rank 0 alone, before distribute() has handed the domain
        out. A collective there would block forever waiting for ranks that are
        not in that code at all.

        Ghost cells are excluded, so summing this across ranks gives the
        whole-domain figure exactly once.
        """
        m = self.tracer_conserved_values[s]           # m = h*c per cell
        areas = self.areas
        if self.tri_full_flag is not None:
            return float(num.sum(m * areas * (self.tri_full_flag == 1)))
        return float(num.sum(m * areas))

    def _reduce_over_ranks(self, value):
        """Sum a per-rank scalar over the communicator. Collective."""
        from anuga import numprocs
        if numprocs == 1:
            return value
        from mpi4py import MPI
        return MPI.COMM_WORLD.allreduce(value, op=MPI.SUM)

    def get_tracer_mass(self, name):
        """Total mass of tracer `name` in the domain: the integral of m = h*c.

        The tracer analogue of `get_water_volume`. Summed over owned cells only
        and reduced across ranks, so it is the whole-domain figure in parallel
        as well as in serial.

        COLLECTIVE in parallel, like `get_water_volume`: every rank must call
        it, or the ones that do will block.
        """
        s = self.get_tracer_index(name)
        return self._reduce_over_ranks(self._local_tracer_mass(s))

    def get_tracer_boundary_flux_integral(self, name):
        """Net tracer mass that has crossed the domain boundary, since t=0.

        Positive is INTO the domain, matching the sign of the water balance.
        The tracer analogue of `get_boundary_flux_integral`, and the other half
        of a mass budget:

            mass(t) - mass(0) == boundary_flux_integral(t)

        to within the timestepping error, for a domain with no source terms.
        `check_tracer_conservation` does exactly that comparison.

        The kernel accumulates the per-edge boundary flux into a per-cell array
        each substep; `tracer_flux_integral_operator` applies the timestepping
        method's own weights and zeroes it, exactly as
        `boundary_flux_integral_operator` does for water.

        Reduced across ranks, matching `get_boundary_flux_integral` for water:
        each rank counts only the boundary edges of the cells it owns, so the
        per-rank integrals sum to the whole-domain figure. COLLECTIVE.
        """
        s = self.get_tracer_index(name)
        return self._reduce_over_ranks(float(self._tracer_flux_integral[s]))

    def check_tracer_conservation(self, name):
        """Return (change_in_mass, boundary_flux_integral, discrepancy).

        For a domain with no tracer source terms the first two should agree, so
        the discrepancy is the conservation error. On a CLOSED domain the flux
        integral is zero and this reduces to "did the mass stay constant".

        Returns absolute quantities, not a relative error: for a tracer that is
        zero almost everywhere a relative measure is meaningless, and the caller
        knows the scale that matters.

        COLLECTIVE in parallel: every rank must call it.
        """
        s = self.get_tracer_index(name)
        # The baseline is a per-rank number (see _local_tracer_mass), so the
        # difference is formed locally and reduced once. Reducing the two
        # separately would be wrong after a distribute(), where each rank's
        # baseline covers only its own cells.
        change = self._reduce_over_ranks(
            self._local_tracer_mass(s) - self._tracer_initial_mass[s])
        flux = self.get_tracer_boundary_flux_integral(name)
        return change, flux, change - flux

    def _reorder_tracer_arrays(self, new_order, inv_order,
                               old_boundary_enumeration):
        """Permute the tracer blocks onto a reordered mesh (issue #277).

        Called by `Generic_Domain.reorder` after the mesh has been renumbered
        but while `inv_order` and the pre-reorder boundary enumeration are
        still available. A no-op on a domain with no tracers.

        Permutes IN PLACE. The C domain struct holds raw pointers into these
        buffers, so rebinding the attributes would leave the struct addressing
        the old memory.
        """
        ns = self.number_of_tracers
        if ns == 0:
            return

        N = self.number_of_elements
        new_order = num.asarray(new_order, dtype=int)

        # Belt and braces: a block added to _TRACER_ARRAYS without a kind would
        # otherwise be skipped here and silently keep the old ordering, which
        # is the exact failure #277 was.
        unclassified = [a for a in self._TRACER_ARRAYS
                        if a not in self._TRACER_ARRAY_KINDS]
        if unclassified:
            raise NotImplementedError(
                'tracer array(s) %s have no entry in _TRACER_ARRAY_KINDS, so '
                'reorder() does not know how to permute them; classify them '
                'as cell/edge/boundary' % ', '.join(sorted(unclassified)))

        bperm = None

        for attr in self._TRACER_ARRAYS:
            arr = getattr(self, attr, None)
            if arr is None:
                continue
            kind = self._TRACER_ARRAY_KINDS[attr]

            if kind == 'cell':
                # new[:, i] = old[:, new_order[i]], matching the convention the
                # quantities use. Fancy indexing copies before the assignment,
                # so this does not alias.
                arr[:] = arr[:, new_order]
            elif kind == 'edge':
                # (ns, 3N) -> (ns, N, 3), permute triangles, flatten back. The
                # same reshape-permute-ravel the (3N,) arrays get in reorder().
                arr[:] = arr.reshape(ns, N, 3)[:, new_order, :].reshape(ns, 3 * N)
            elif kind == 'boundary':
                if bperm is None:
                    bperm = self._boundary_permutation(old_boundary_enumeration,
                                                       inv_order)
                arr[:] = arr[:, bperm]
            else:
                raise NotImplementedError(
                    'unknown tracer array kind %r for %r' % (kind, attr))

        # The bedload boundary flag is per boundary edge too, but is rebuilt
        # from its tags rather than permuted: tag_boundary_cells has already
        # been renumbered by the time this runs.
        if self.sediment_bedload_open is not None:
            self._build_bedload_open()

    def _boundary_permutation(self, old_boundary_enumeration, inv_order):
        """Map each NEW boundary index to the old index of the same edge.

        `build_boundary_neighbours` numbers boundary edges by their position in
        the sorted (triangle, edge) order, so renumbering triangles renumbers
        the boundary as well. Quantities do not care -- `update_boundary`
        refills their boundary values from the Boundary objects every timestep
        -- but `tracer_boundary_values` is written once by
        `set_tracer_boundary` and then just read, so it has to be carried
        across by hand.

        Returns an index array `bperm` with `bperm[new_j] = old_j`, so that
        `new_values = old_values[bperm]`.
        """
        new_boundary_enumeration = self.mesh.boundary_enumeration
        M = len(new_boundary_enumeration)
        if len(old_boundary_enumeration) != M:
            raise RuntimeError(
                'boundary edge count changed during reorder (%d -> %d); the '
                'tracer boundary values cannot be carried across'
                % (len(old_boundary_enumeration), M))

        bperm = num.empty(M, dtype=int)
        for (old_id, edge), old_j in old_boundary_enumeration.items():
            # The same physical edge, addressed by the triangle's new number.
            new_j = new_boundary_enumeration[int(inv_order[old_id]), edge]
            bperm[new_j] = old_j
        return bperm

    def get_tracer_index(self, name):
        """Return the row index of tracer `name`."""
        try:
            return self._tracer_names.index(name)
        except ValueError:
            raise ValueError('no tracer named %r; registered tracers are %r'
                             % (name, list(self._tracer_names)))

    def set_tracer_source(self, name, values):
        """Prescribe an external source `S_ms` for tracer `name`, in m/s.

        `[G-3]`'s optional external supply: hillslope yield, tributary load,
        rainfall washoff. Same units as `E` and `D`.

        Applied **after** the `[L-1]`/`[L-2]` limiters. Those bound the bed
        exchange by what the bed and water column can supply; an external
        source is neither, so clipping it there would be wrong -- and would
        also make a manufactured solution impossible to impose exactly.
        """
        s = self.get_tracer_index(name)
        # Writes only. The array is allocated by add_tracer, so the pointer the
        # C struct and the device hold stays valid and this is safe to call
        # from inside a fractional step -- which a manufactured source, or any
        # time-varying supply, does on every timestep. See #288.
        v = num.asarray(values, dtype=num.float64)
        if v.ndim == 0:
            self.tracer_external_source[s] = float(v)
        elif v.shape == (self.number_of_elements,):
            self.tracer_external_source[s] = v
        else:
            raise ValueError('tracer %r: expected a scalar or %d values, got %r'
                             % (name, self.number_of_elements, v.shape))

        # In mode 2 the kernel reads the DEVICE copy, so a host write has to be
        # pushed or the source goes stale -- silently, and stale is worse than
        # absent. One array, so a per-step source costs one small transfer.
        if (self.multiprocessor_mode == MULTIPROCESSOR_GPU
                and self.gpu_interface is not None):
            from anuga.shallow_water.sw_domain_gpu_ext import (
                sync_tracer_source_to_device)
            sync_tracer_source_to_device(self.gpu_interface.gpu_dom)

    def get_tracer_names(self):
        """Return the registered tracer names, in index order."""
        return list(self._tracer_names)

    def _tracer_depth(self):
        """Cell depth h = stage - elevation, floored at zero."""
        h = (self.quantities['stage'].centroid_values
             - self.quantities['elevation'].centroid_values)
        return num.maximum(h, 0.0)

    def get_tracer(self, name):
        """Return the concentration `c` of tracer `name`, one value per cell."""
        return self.tracer_centroid_values[self.get_tracer_index(name)]

    def set_tracer(self, name, values):
        """Set the concentration `c` of tracer `name`.

        Seeds both `c` and the conserved `m = h*c`, which must agree: setting
        `c` alone leaves the conserved variable stale and the first substep
        overwrites `c` from it.

        Because `m = h*c`, this reads the CURRENT depth: set `stage` and
        `elevation` before the tracer, not after.

        Called before the run starts, this also rebases the conservation
        baseline (see `check_tracer_conservation`); called mid-run it does not.
        """
        s = self.get_tracer_index(name)
        # asarray, not ascontiguousarray: the latter promotes a 0-d scalar to
        # shape (1,), which would defeat the scalar branch below.
        c = num.asarray(values, dtype=num.float64)
        if c.ndim == 0:
            c = num.full(self.number_of_elements, float(c))
        elif c.shape != (self.number_of_elements,):
            raise ValueError(
                'tracer %r: expected a scalar or %d values, got shape %r'
                % (name, self.number_of_elements, c.shape))
        self.tracer_centroid_values[s] = c
        self.tracer_conserved_values[s] = self._tracer_depth() * c

        # Setting the field before the run starts IS the initial condition, so
        # rebase the conservation baseline on it. Without this, the common
        # add_tracer(name) / set_tracer(name, field) pairing would leave the
        # baseline at the mass add_tracer seeded and check_tracer_conservation
        # would report the difference as a conservation error.
        #
        # Mid-run a set_tracer is a genuine intervention that really does break
        # the budget, so leave the baseline alone there and let it show up as a
        # discrepancy rather than silently absorbing it.
        if (self._tracer_initial_mass is not None
                and s < len(self._tracer_initial_mass)
                and self.relative_time == self.evolve_starttime):
            self._tracer_initial_mass[s] = self._local_tracer_mass(s)

    def update_domain_c_struct(self):
        """Update the C domain structure from the Python Domain object.
        """
        from .sw_domain_openmp_ext import update_Domain_C_struct
        update_Domain_C_struct(self)

    def _ensure_work_arrays(self):
        """Allocate the work arrays actually needed by the C extension.

        Called by ``setup_Domain_C_struct`` in the Cython extension before the
        C domain struct is populated.  Only arrays that the C computation code
        actually reads are allocated here — confirmed-dead arrays remain None
        (passed as NULL to the C struct) for the full simulation lifetime.
        """
        if self.x_centroid_work is not None:
            return  # already allocated

        N  = self.number_of_elements

        # ---- scratch arrays used by the C extrapolation kernel ----
        self.x_centroid_work   = num.zeros(N)               # (N,)
        self.y_centroid_work   = num.zeros(N)               # (N,)

        # ---- per-element max wave speed (CFL timestep calculation) ----
        self.max_speed         = num.zeros(N, float)        # (N,)

        # Note: edge_flux_type / edge_river_wall_counter stay None unless
        # create_riverwalls() is called; all other struct fields remain NULL.

    #---------------------------------------------------------------
    # Plotting methods
    #---------------------------------------------------------------
    def _plotter_quantity(self, name):
        """A derived quantity of the attached plotter, computed on access."""
        if self.dplotter is None:
            raise AttributeError(
                "domain.%s needs a plotter: call domain.set_plotter() first "
                "(it returns the plotter, so dplotter = domain.set_plotter() "
                "gives you dplotter.%s as well)" % (name, name))
        return getattr(self.dplotter, name)

    @property
    def depth(self):
        """Depth at the plotter's vertices, `stage - elevation`."""
        return self._plotter_quantity('depth')

    @property
    def xvel(self):
        """x velocity at the plotter's vertices, zero below `min_depth`."""
        return self._plotter_quantity('xvel')

    @property
    def yvel(self):
        """y velocity at the plotter's vertices, zero below `min_depth`."""
        return self._plotter_quantity('yvel')

    @property
    def speed(self):
        """Speed at the plotter's vertices, `hypot(xvel, yvel)`."""
        return self._plotter_quantity('speed')

    @property
    def speed_depth(self):
        """Speed times depth at the plotter's vertices."""
        return self._plotter_quantity('speed_depth')

    def set_plotter(self, *args, **kwargs):
        """Attach a `Domain_plotter` to this domain and return it.

        The plotter is created on the first call and kept on the domain as
        `domain.dplotter`; later calls return the one already there and
        ignore their arguments. Its arrays and plotting methods are also
        mirrored onto the domain, so `domain.plot_depth_frame()` and
        `plotter.plot_depth_frame()` are the same call.

        Parameters
        ----------
        *args, **kwargs
            Passed to :class:`~anuga.Domain_plotter` on the first call, e.g.
            `absolute=True` to plot in the domain's absolute coordinates.

        Returns
        -------
        Domain_plotter
            The plotter, so it can be named directly::

                dplotter = domain.set_plotter(absolute=True)

            rather than reaching for `domain.dplotter` afterwards.
        """

        #FIXME SR: Should look into seeing if the triang can use the
        # triangulation from Domain rather than having two copies
        if self.dplotter is None:
            import anuga
            self.dplotter = anuga.Domain_plotter(self, *args, **kwargs)


        self.triang = self.dplotter.triang
        self.stage = self.dplotter.stage
        self.xmom = self.dplotter.xmom
        self.ymom = self.dplotter.ymom
        self.elev = self.dplotter.elev
        self.friction = self.dplotter.friction
        # depth, xvel, yvel, speed and speed_depth are NOT copied here: they
        # are properties on the plotter that recompute from stage and
        # elevation, so copying froze them at the moment set_plotter was
        # called. They are properties of the domain too, just below.
        self.x = self.dplotter.x
        self.y = self.dplotter.y
        self.xc = self.dplotter.xc
        self.yc = self.dplotter.yc

        self.plot_mesh = self.dplotter.plot_mesh

        self.save_depth_frame = self.dplotter.save_depth_frame
        self.plot_depth_frame = self.dplotter.plot_depth_frame
        self.make_depth_animation = self.dplotter.make_depth_animation

        self.save_stage_frame = self.dplotter.save_stage_frame
        self.plot_stage_frame = self.dplotter.plot_stage_frame
        self.make_stage_animation = self.dplotter.make_stage_animation

        self.save_speed_frame = self.dplotter.save_speed_frame
        self.plot_speed_frame = self.dplotter.plot_speed_frame
        self.make_speed_animation = self.dplotter.make_speed_animation

        return self.dplotter

    def triplot(self, *args,  **kwargs):

        self.set_plotter()

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

        lines = ax.triplot(self.triang, *args, **kwargs)

        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')

        return fig, ax, lines


    def tripcolor(self, *args,  **kwargs):

        self.set_plotter()

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

        im = ax.tripcolor(self.triang,  *args, **kwargs)

        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')

        return fig, ax, im

    #==============================================================
    # Methods to set and get domain parameters
    #
    # FIXME SR: These (and other paramters) should be refactored
    # to save the underlying quantities in np.ndarray(s) for
    # efficient access in Cython
    #==============================================================

    @property
    def g(self) -> float:
        """Gravitational acceleration [m/s^2]"""
        return self._g

    @g.setter
    def g(self, value: float):
        """Set gravitational acceleration [m/s^2]"""
        self._g = value

    @property
    def timestep(self) -> float:
        """Current timestep [s]"""
        return self._timestep

    @timestep.setter
    def timestep(self, value: float):
        """Set current timestep [s]"""
        self._timestep = value

    @property
    def flux_timestep(self) -> float:
        """Current flux timestep [s]"""
        return self._flux_timestep

    @flux_timestep.setter
    def flux_timestep(self, value: float):
        """Set current flux timestep [s]"""
        self._flux_timestep = value


    #==============================================================
    # Set config defaults
    #==============================================================
    def _set_config_defaults(self):
        """Set the default values in this routine. That way we can inherit class
        and just redefine the defaults for the new class
        """

        from anuga.config import minimum_storable_height
        from anuga.config import minimum_allowed_height, maximum_allowed_speed
        from anuga.config import g
        from anuga.config import tight_slope_limiters
        from anuga.config import extrapolate_velocity_second_order
        from anuga.config import alpha_balance
        from anuga.config import optimise_dry_cells
        from anuga.config import use_centroid_velocities
        from anuga.config import compute_fluxes_method
        from anuga.config import sloped_mannings_function
        from anuga.config import low_froude


        # Early algorithms need elevation to remain continuous
        self.set_using_discontinuous_elevation(False)

        self.set_minimum_allowed_height(minimum_allowed_height)
        self.maximum_allowed_speed = maximum_allowed_speed

        from anuga.config import negative_volume_warning_fraction
        self.negative_volume_warning_fraction = negative_volume_warning_fraction

        self.minimum_storable_height = minimum_storable_height

        self.g = g

        self.alpha_balance = alpha_balance
        self.tight_slope_limiters = tight_slope_limiters

        self.set_low_froude(low_froude)

        self.set_use_optimise_dry_cells(optimise_dry_cells)
        self.set_extrapolate_velocity(extrapolate_velocity_second_order)

        self.use_centroid_velocities = use_centroid_velocities

        self.set_sloped_mannings_function(sloped_mannings_function)
        self.set_compute_fluxes_method(compute_fluxes_method)

        self.set_store_centroids(False)

    def get_algorithm_parameters(self) -> dict:
        """
        Get the standard parameter that are currently set (as a dictionary)
        """

        parameters = {}

        parameters['minimum_allowed_height']  = self.minimum_allowed_height
        parameters['maximum_allowed_speed']   = self.maximum_allowed_speed
        parameters['minimum_storable_height'] = self.minimum_storable_height
        parameters['g']                       = self.g
        parameters['alpha_balance']           = self.alpha_balance
        parameters['tight_slope_limiters']    = self.tight_slope_limiters
        parameters['optimise_dry_cells']      = self.optimise_dry_cells
        parameters['low_froude']              = self.low_froude
        parameters['use_centroid_velocities'] = self.use_centroid_velocities
        parameters['use_sloped_mannings']     = self.use_sloped_mannings
        parameters['compute_fluxes_method']   = self.get_compute_fluxes_method()
        parameters['flow_algorithm']          = self.get_flow_algorithm()
        parameters['CFL']                     = self.get_cfl()
        parameters['timestepping_method']     = self.get_timestepping_method()

        parameters['extrapolate_velocity_second_order'] = self.extrapolate_velocity_second_order



        return parameters

    def print_algorithm_parameters(self) -> None:
        """
        Print the standard parameters that are curently set (as a dictionary)
        """

        print('#============================')
        print('# Domain Algorithm Parameters ')
        print('#============================')
        from pprint import pprint
        pprint(self.get_algorithm_parameters(),indent=4)

        print('#----------------------------')




    def _set_DE0_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE0"
           A 'discontinuous elevation' method
        """

        self._set_config_defaults()

        self.set_cfl(0.9)
        self.set_use_kinematic_viscosity(False)
        #self.timestepping_method='rk2'#'rk3'#'euler'#'rk2'
        self.set_timestepping_method('euler')

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        # Don't place any restriction on the minimum storable height
        #self.minimum_storable_height=-99999999999.0
        self.minimum_allowed_height=1.0e-12

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        # beta 0.5 is not a preference, it is what keeps forward Euler
        # stable.  Raising it to 1.0 makes DE0 go unstable on deep_wave (the
        # stage reaches 25-29 against a forced amplitude of 1.0), destroys
        # monotonicity at the dry dam break (total variation 112 against the
        # exact 10.0, with a 0.43 m overshoot where beta 0.5 has none) and
        # degrades the dry avalanche (RMS 0.29 against 0.051).  Do not raise
        # it here; DE1 and DE_ader2 pair beta 1.0 with a second-order time
        # integrator that can carry it.
        #
        # The price is damping.  DE0 keeps 66% of the forced amplitude across
        # the deep_wave domain against DE1's 97%, so it is a poor choice for
        # wave propagation.  beta < 1 also cannot reconstruct a linear stage
        # exactly -- the limiter scales the gradient by beta on genuinely
        # linear data -- which flattens the bed the solver works with
        # (z = w - h) and costs steady slope flow about 12% of its discharge.
        # See docs/source/setup_anuga_script/flow_algorithms.rst, and
        # _set_DE_ader2_defaults() for the case where that reconstruction
        # error broke a steady state outright.
        self.beta_w=0.5
        self.beta_w_dry=0.0
        self.beta_uh=0.5
        self.beta_uh_dry=0.0
        self.beta_vh=0.5
        self.beta_vh_dry=0.0


        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 2, 'height':2})
        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 1})
        self.set_store_centroids(True)

        self.optimise_dry_cells=False

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # By default vertex values are NOT stored uniquely
        # for storage efficiency. We may override this (but not so important since
        # centroids are stored anyway
        # self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0

        if self.processor == 0 and self.verbose:
            print('Domain: Using discontinuous elevation solver DE0')
            # print('##########################################################################')
            # print('#')
            # print('# Using discontinuous elevation solver DE0')
            # print('#')
            # print('# First order timestepping')
            # print('#')
            # print('# Make sure you use centroid values when reporting on important output quantities')
            # print('#')
            # print('##########################################################################')


    def _set_DE1_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE1"
           A 'discontinuous elevation' method
        """

        self._set_config_defaults()

        self.set_cfl(1.0)
        self.set_use_kinematic_viscosity(False)
        #self.timestepping_method='rk2'#'rk3'#'euler'#'rk2'
        self.set_timestepping_method(2)

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        # Don't place any restriction on the minimum storable height
        #self.minimum_storable_height=-99999999999.0
        self.minimum_allowed_height=1.0e-5

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        self.beta_w=1.0
        self.beta_w_dry=0.0
        self.beta_uh=1.0
        self.beta_uh_dry=0.0
        self.beta_vh=1.0
        self.beta_vh_dry=0.0


        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 2, 'height':2})
        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 1})
        self.set_store_centroids(True)

        self.optimise_dry_cells=False

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # By default vertex values are NOT stored uniquely
        # for storage efficiency. We may override this (but not so important since
        # centroids are stored anyway
        # self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0


        if self.processor == 0 and self.verbose:
            print('##########################################################################')
            print('#')
            print('# Using discontinuous elevation solver DE1 ')
            print('#')
            print('# Uses rk2 timestepping')
            print('#')
            print('# Make sure you use centroid values when reporting on important output quantities')
            print('#')
            print('##########################################################################')

    def _set_DE_ader2_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE_ader2"
           DE1 settings with ADER-2 (Cauchy-Kovalewski) timestepping.
        """

        self._set_config_defaults()

        self.set_cfl(1.0)
        self.set_use_kinematic_viscosity(False)
        self.set_timestepping_method('ader2')

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        self.minimum_allowed_height = 1.0e-5

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        # beta 1.0, not the 0.5 DE0 uses.  Both of the cases that pin this
        # down want it: on the deep-water sinusoid, beta 0.5 keeps only 32% of
        # the wave amplitude across the domain where beta 1.0 keeps 100%; and
        # on rundown_mild_slope_coarse, beta 0.5 holds 0.880 of the prescribed
        # discharge where beta 1.0 holds 0.9997, matching DE1.  The 0.5 error
        # on the slope is the limiter's, not the ADER step's -- DE0 at beta
        # 0.5 gives 0.877 there too -- and it only became visible once the
        # predictor was made steady-state exact (see
        # core_ader_friction_factor() in core_kernels.c).
        self.beta_w = 1.0
        self.beta_w_dry = 0.0
        self.beta_uh = 1.0
        self.beta_uh_dry = 0.0
        self.beta_vh = 1.0
        self.beta_vh_dry = 0.0

        self.set_store_centroids(True)

        self.optimise_dry_cells = False

        self.edge_coordinates = self.get_edge_midpoint_coordinates()

        self.maximum_allowed_speed = 0.0

        if self.processor == 0 and self.verbose:
            print('##########################################################################')
            print('#')
            print('# Using discontinuous elevation solver DE_ader2')
            print('#')
            print('# Uses ADER-2 (Cauchy-Kovalewski) timestepping')
            print('#')
            print('# Make sure you use centroid values when reporting on important output quantities')
            print('#')
            print('##########################################################################')

    def _set_DE2_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE2"
           A 'discontinuous elevation' method
        """

        self._set_config_defaults()

        self.set_cfl(1.0)
        self.set_use_kinematic_viscosity(False)
        #self.timestepping_method='rk2'#'rk3'#'euler'#'rk2'
        self.set_timestepping_method(3)

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        # Don't place any restriction on the minimum storable height
        #self.minimum_storable_height=-99999999999.0
        self.minimum_allowed_height=1.0e-5

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        self.beta_w=1.0
        self.beta_w_dry=0.0
        self.beta_uh=1.0
        self.beta_uh_dry=0.0
        self.beta_vh=1.0
        self.beta_vh_dry=0.0


        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 2, 'height':2})
        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 1})
        self.set_store_centroids(True)

        self.optimise_dry_cells=False

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # By default vertex values are NOT stored uniquely
        # for storage efficiency. We may override this (but not so important since
        # centroids are stored anyway
        # self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0


        if self.processor == 0 and self.verbose:
            print('##########################################################################')
            print('#')
            print('# Using discontinuous elevation solver DE2')
            print('#')
            print('# Using rk3 timestepping')
            print('#')
            print('# Make sure you use centroid values when reporting on important output quantities')
            print('#')
            print('##########################################################################')

    def _set_DE1_7_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE0_7"
           A 'discontinuous elevation' method
        """

        self._set_config_defaults()

        self.set_cfl(1.0)
        self.set_use_kinematic_viscosity(False)
        #self.timestepping_method='rk2'#'rk3'#'euler'#'rk2'
        self.set_timestepping_method(2)

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        # Don't place any restriction on the minimum storable height
        #self.minimum_storable_height=-99999999999.0
        self.minimum_allowed_height=1.0e-12

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        self.beta_w=0.75
        self.beta_w_dry=0.1
        self.beta_uh=0.75
        self.beta_uh_dry=0.1
        self.beta_vh=0.75
        self.beta_vh_dry=0.1


        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 2, 'height':2})
        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 1})
        self.set_store_centroids(True)

        self.optimise_dry_cells=False

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # By default vertex values are NOT stored uniquely
        # for storage efficiency. We may override this (but not so important since
        # centroids are stored anyway
        # self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0

        if self.processor == 0 and self.verbose:
            print('##########################################################################')
            print('#')
            print('# Using discontinuous elevation solver DE1_7 ')
            print('#')
            print('# A slightly more diffusive version of DE1, does use rk2 timestepping')
            print('#')
            print('# Make sure you use centroid values when reporting on important output quantities')
            print('#')
            print('##########################################################################')


    def _set_DE0_7_defaults(self):
        """Set up the defaults for running the flow_algorithm "DE3"
           A 'discontinuous elevation' method
        """

        self._set_config_defaults()

        self.set_cfl(0.9)
        self.set_use_kinematic_viscosity(False)
        #self.timestepping_method='rk2'#'rk3'#'euler'#'rk2'
        self.set_timestepping_method(1)

        self.set_using_discontinuous_elevation(True)
        self.set_using_centroid_averaging(True)
        self.set_compute_fluxes_method('DE')

        # Don't place any restriction on the minimum storable height
        #self.minimum_storable_height=-99999999999.0
        self.minimum_allowed_height=1.0e-12

        self.set_default_order(2)
        self.set_extrapolate_velocity()

        self.beta_w=0.7
        self.beta_w_dry=0.1
        self.beta_uh=0.7
        self.beta_uh_dry=0.1
        self.beta_vh=0.7
        self.beta_vh_dry=0.1


        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 2, 'height':2})
        #self.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2,
        #         'ymomentum': 2, 'elevation': 1})
        self.set_store_centroids(True)

        self.optimise_dry_cells=False

        # We need the edge_coordinates for the extrapolation
        self.edge_coordinates=self.get_edge_midpoint_coordinates()

        # By default vertex values are NOT stored uniquely
        # for storage efficiency. We may override this (but not so important since
        # centroids are stored anyway
        # self.set_store_vertices_smoothly(False)

        self.maximum_allowed_speed=0.0

        if self.processor == 0 and self.verbose:
            print('Domain: Using discontinuous elevation solver DE0_7')
            # print('##########################################################################')
            # print('#')
            # print('# Using discontinuous elevation solver DE0_7')
            # print('#')
            # print('# A slightly less diffusive version than DE0, uses euler timestepping')
            # print('#')
            # print('# Make sure you use centroid values when reporting on important output quantities')
            # print('#')
            # print('##########################################################################')



    def update_special_conditions(self):

        my_update_special_conditions(self)

    # Note Padarn 06/12/12: The following line decorates
    # the set_quantity function to be profiled individually.
    # Need to uncomment the decorator at top of file.
    #@profileit("set_quantity.profile")
    def set_quantity(self, name: str, *args, **kwargs) -> None:
        """Set values for named quantity

        We have to do something special for 'elevation'
        otherwise pass through to generic set_quantity

        Mode-2 ('unified'): the device holds the authoritative centroid state once
        the GPU interface exists, so a host-only write has to be mirrored to it.
        That is handled one level down, in Quantity.set_values(), via the
        _notify_*_host_quantity_write() hooks below — which also covers callers
        that reach a Quantity directly and bypass this method.
        """

#        if name == 'elevation':
#            stage_c = self.get_quantity('stage').centroid_values
#            elev_c =  self.get_quantity('elevation').centroid_values
#            height_c = stage_c - elev_c
#            Generic_Domain.set_quantity(self, name, *args, **kwargs)
#            stage_c[:] = elev_c + height_c
#        else:
#            Generic_Domain.set_quantity(self, name, *args, **kwargs)

        Generic_Domain.set_quantity(self, name, *args, **kwargs)

    # ------------------------------------------------------------------
    # Mode-2 host->device coherency for quantity writes
    #
    # Once the GPU interface exists the *device* holds the authoritative centroid
    # state.  A write that touches only the host arrays is then silently ignored —
    # the next step reads the stale device values, so the simulation runs with the
    # wrong data.  This bites whenever something builds the interface *before* the
    # quantities are set (distribute_to_vertices_and_edges() and set_boundary()
    # both call _ensure_gpu_interface()), and for any mid-run write.
    #
    # Quantity.set_values() is the single choke point for host-side quantity
    # writes, so the sync hangs off there rather than off Domain.set_quantity() —
    # that way a caller holding a Quantity directly is covered too.
    # ------------------------------------------------------------------

    # Set while apply_fractional_steps() already brackets its operators with a
    # sync_from_device()/sync_to_device() pair, so operator writes inside that
    # region don't re-sync once per call.
    _gpu_host_writes_suppressed = False

    def _gpu_syncs_host_quantity_write(self, name: str) -> bool:
        """True if a host write to quantity `name` has to be mirrored to the device."""

        return (not self._gpu_host_writes_suppressed
                and name in GPU_SYNCED_QUANTITIES
                and self.multiprocessor_mode == MULTIPROCESSOR_GPU
                and getattr(self, 'gpu_interface', None) is not None)

    def _notify_before_host_quantity_write(self, name: str) -> None:
        """Refresh the host centroids from the device, ahead of a host-side write.

        Without this the push in _notify_after_host_quantity_write() would send
        stale host values for the *other* quantities over the device's current ones.
        """

        if self._gpu_syncs_host_quantity_write(name):
            try:
                self.gpu_interface.sync_from_device()
            except Exception:
                pass

    def _notify_after_host_quantity_write(self, name: str) -> None:
        """Push a host-side quantity write out to the device."""

        if self._gpu_syncs_host_quantity_write(name):
            self.gpu_interface.sync_to_device()


    def _sync_riverwall_to_device(self) -> None:
        """Push pending host-side riverwall changes out to the device (mode 2).

        The riverwall crest elevations and hydraulic properties are mapped to the
        device once, when the mode-2 interface is built, and are never written
        there. A host-side change — RiverWall.set_elevation() to operate a gate,
        say — therefore has no effect on the device until it is pushed across,
        which is what this does.

        Called by evolve() at yieldstep boundaries only: that is where a script
        can make such a change (in the body of the evolve loop) and where the
        host and device are already in step, so the crest never changes part-way
        through a timestep or between RK substeps.
        """

        riverwall_data = getattr(self, 'riverwallData', None)
        if riverwall_data is None:
            return
        if not getattr(riverwall_data, 'device_data_dirty', False):
            return

        if (self.multiprocessor_mode == MULTIPROCESSOR_GPU
                and getattr(self, 'gpu_interface', None) is not None):
            self.gpu_interface.sync_riverwall_to_device()

        # Cleared unconditionally: on the CPU paths the kernels read the host
        # arrays directly, and a later mode-2 setup maps their current values.
        riverwall_data.device_data_dirty = False

    def set_timezone(self, tz: str | ZoneInfoType | None = None) -> None:
        """Set timezone for domain

        :param tz: either a timezone object or string

        We recommend using the ZoneInfo provided by zoneinfo.
        Default is ZoneInfo('UTC')

        Example: Set default timezone UTC

        >>> domain.set_timezone()

        Example: Set timezone using tsdata string

        >>> domain.set_timezone('Australia/Syndey')

        Example: Set timezone using ZoneInfo timezone

        >>> from zoneinfo import ZoneInfo
        >>> new_tz = ZoneInfo('Australia/Sydney')
        >>> domain.set_timezone(new_tz)
        """

        try:
            from zoneinfo import ZoneInfo
        except ImportError:
            from backports.zoneinfo import ZoneInfo

        if tz is None:
            new_tz = ZoneInfo('UTC')
        elif isinstance(tz,str):
            new_tz = ZoneInfo(tz)
        elif isinstance(tz, ZoneInfo):
            new_tz = tz
        else:
            msg = "Unknown timezone %s" % tz
            raise Exception(msg)


        self.timezone = new_tz

    def get_timezone(self) -> ZoneInfoType:
        """Retrieve current domain timezone"""

        return self.timezone

    def get_datetime(self, timestamp: float | None = None) -> DateTime:
        """Retrieve datetime corresponding to current timestamp wrt to domain timezone

        param: timestamp: return datetime corresponding to given timestamp"""

        from datetime import datetime

        try:
            from datetime import UTC
        except ImportError:
            from datetime import timezone
            UTC = timezone.utc


        if timestamp is None:
            timestamp = self.get_time()

        #utc_datetime = datetime.utcfromtimestamp(timestamp).replace(tzinfo=ZoneInfo('UTC'))
        utc_datetime = datetime.fromtimestamp(timestamp,UTC)

        current_dt = utc_datetime.astimezone(self.timezone)
        return current_dt

    def set_starttime(self, timestamp: float | DateTime = 0.0) -> None:
        """Set the starttime for the evolution

        :param timestamp: Either a float or a datetime object

        Essentially we use unix time as our absolute time. So
        time = 0 corresponds to Jan 1st 1970 UTC

        Use naive datetime which will be localized to the domain timezone or
        or use zoneinfo.ZoneInfo to set the timezone of datetime.
        Don't use the tzinfo argument of datetime to set timezone as this does not work!

        Example:

            Without setting timezone for the `domain` and the `starttime` then time
            calculations are all based on UTC. Note the timestamp, which is time in seconds
            from 1st Jan 1970 UTC.

        >>> import anuga
        >>> from zoneinfo import ZoneInfo
        >>> from datetime import datetime
        >>>
        >>> domain = anuga.rectangular_cross_domain(10,10)
        >>> dt = datetime(2021,3,21,18,30)
        >>> domain.set_starttime(dt)
        >>> print(domain.get_datetime(), 'TZ', domain.get_timezone(), 'Timestamp: ', domain.get_time())
        2021-03-21 18:30:00+00:00 TZ UTC Timestamp:  1616351400.0

        Example:

            Setting timezone for the `domain`, then naive `datetime` will be localizes to
            the `domain` timezone. Note the timestamp, which is time in seconds
            from 1st Jan 1970 UTC.

        >>> import anuga
        >>> from zoneinfo import ZoneInfo
        >>> from datetime import datetime
        >>>
        >>> domain = anuga.rectangular_cross_domain(10,10)
        >>> AEST = ZoneInfo('Australia/Sydney')
        >>> domain.set_timezone(AEST)
        >>>
        >>> dt = datetime(2021,3,21,18,30)
        >>> domain.set_starttime(dt)
        >>> print(domain.get_datetime(), 'TZ', domain.get_timezone(), 'Timestamp: ', domain.get_time())
        2021-03-21 18:30:00+11:00 TZ Australia/Sydney Timestamp:  1616311800.0

        Example:

            Setting timezone for the `domain`, and setting the timezone for the `datetime`.
            Note the timestamp, which is time in seconds from 1st Jan 1970 UTC is the same
            as the previous example.

        >>> import anuga
        >>> from zoneinfo import ZoneInfo
        >>> from datetime import datetime
        >>>
        >>> domain = anuga.rectangular_cross_domain(10,10)
        >>>
        >>> ACST = ZoneInfo('Australia/Adelaide')
        >>> domain.set_timezone(ACST)
        >>>
        >>> AEST = ZoneInfo('Australia/Sydney')
        >>> dt = datetime(2021,3,21,18,30, tzinfo=AEST)
        >>>
        >>> domain.set_starttime(dt)
        >>> print(domain.get_datetime(), 'TZ', domain.get_timezone(), 'Timestamp: ', domain.get_time())
        2021-03-21 18:00:00+10:30 TZ Australia/Adelaide Timestamp:  1616311800.0
        """


        from datetime import datetime

        if self.evolved_called:
            msg = ('Can\'t change simulation start time once evolve has '
                   'been called')
            raise Exception(msg)

        if isinstance(timestamp, datetime):
            if timestamp.tzinfo is None:
                dt = timestamp.replace(tzinfo=self.timezone)
                time = dt.timestamp()
            else:
                time = timestamp.timestamp()
        else:
            time = float(timestamp)


        self.starttime = time
        # starttime is now the origin for relative_time
        self.set_relative_time(0.0)

    def get_starttime(self, datetime: bool = False) -> float | DateTime:
        """return starttime, either as timestamp, or as a datetime"""

        starttime = self.starttime

        if not datetime:
            return starttime
        else:
            return self.get_datetime(starttime)


    def set_store(self, flag: bool = True) -> None:
        """Set whether data saved to sww file.
        """

        self.store = flag

    def get_store(self) -> bool:
        """Get whether data saved to sww file.
        """

        return self.store


    def set_store_centroids(self, flag: bool = True) -> None:
        """Set whether centroid data is saved to sww file.
        """

        self.store_centroids = flag

    def get_store_centroids(self) -> bool:
        """Get whether data saved to sww file.
        """

        return self.store_centroids

    def set_checkpointing(self, checkpoint: bool = True, checkpoint_dir: str = 'CHECKPOINTS',
                          checkpoint_step: int = 10, checkpoint_time: float | None = None) -> None:
        """Set up checkpointing.

        param checkpoint: Default = True. Set to False will turn off checkpointing
        param checkpoint_dir: Where to store checkpointing files
        param checkpoint_step: Save checkpoint files after this many yieldsteps
        param checkpoint_time: If set, over-rides checkpoint_step. save checkpoint files
        after this amount of walltime
        """



        if checkpoint:

            from anuga import myid
            # On processor 0 create checkpoint directory if necessary
            if myid == 0:
                if True:
                    if not os.path.exists(checkpoint_dir):
                        os.mkdir(checkpoint_dir)

                    assert os.path.exists(checkpoint_dir)

            self.checkpoint_dir = checkpoint_dir
            if checkpoint_time is not None:
                #import time
                self.walltime_prev = time.time()
                self.checkpoint_time = checkpoint_time
                self.checkpoint_step = 0
            else:
                self.checkpoint_step = checkpoint_step
            self.checkpoint = True
            #print(self.checkpoint_dir, self.checkpoint_step)
        else:
            self.checkpoint = False

    VEGETATION_FORMULATIONS = {'off': 0, 'baptist': 1}

    VEGETATION_SEDIMENT_SHEAR = {'ignore': 0, 'bed': 1, 'total': 2}

    def set_vegetation_drag(self, density=None, diameter=None, height=None,
                            formulation='baptist', Cd=1.68, bed_chezy=65.0,
                            sediment_shear='bed'):
        """Add the drag of a field of stems to the friction (spec 8).

        The vegetated cells get a friction slope from the Chezy coefficient
        of Baptist et al. (2007),

            Cv = (Cb^-2 + Cd m D min(h, hv) / (2 g))^-1/2
                 + (sqrt(g) / kappa) ln(max(h, hv) / hv),

        applied semi-implicitly like Manning friction and on top of it. Cells
        with no stems are untouched, so open water keeps its Manning `n` and
        vegetation carries `(m, D, hv)`, usually with `n` = 0 there. Stems
        taller than the flow are emergent; when the flow overtops them the
        logarithmic term adds the overflow.

        Parameters
        ----------
        density, diameter, height : float, array or callable
            Stem density `m` (stems per m^2), stem diameter `D` (m) and stem
            height `hv` (m), given as `set_quantity` accepts them: a constant,
            one value per cell centroid, or a function of x and y. They become
            the quantities `veg_density`, `veg_diameter` and `veg_height`.
            A cell with `m` = 0 or `hv` = 0 has no vegetation.
        formulation : str
            `'baptist'` (the only one so far) or `'off'`.
        Cd : float
            Stem drag coefficient, 1.68 (Delta-X Wax Lake model).
        bed_chezy : float
            Bed Chezy coefficient `Cb` inside the stems, 65.
        sediment_shear : {'bed', 'total', 'ignore'}
            The bed shear the sediment kernel sees in vegetated cells
            (entrainment, the Rouse near-bed profile and bedload alike).
            `'bed'` (default): the bed's share of the resistance, as Baptist
            et al. split it -- the stems take the drag and the bed feels the
            canopy velocity `u_v = U Cv_r / Cv` on the bed roughness,
            `tau_b / rho = g u_v^2 / Cb^2` (`u_v = U` for emergent stems).
            `'total'`: the whole vegetated resistance, `tau_b / rho = g U^2 / Cv^2`,
            which is what the Delta-X Wax Lake sediment model uses.
            `'ignore'`: the cell's own sediment friction closure, which with
            `n` = 0 in the vegetated classes means no shear at all.

        Notes
        -----
        Runs in the friction kernel of both compute modes; on a GPU build
        the three fields are mapped to the device once, so call this before
        the first evolve or the device copy is stale. On a distributed run
        call it after `distribute()`, on every rank, as for `set_quantity`.
        """
        if formulation not in self.VEGETATION_FORMULATIONS:
            raise ValueError('unknown vegetation drag formulation %r; expected one of %r'
                             % (formulation, sorted(self.VEGETATION_FORMULATIONS)))
        if sediment_shear not in self.VEGETATION_SEDIMENT_SHEAR:
            raise ValueError('sediment_shear must be one of %r, got %r'
                             % (sorted(self.VEGETATION_SEDIMENT_SHEAR), sediment_shear))
        self.sediment_vegetation_shear = self.VEGETATION_SEDIMENT_SHEAR[sediment_shear]
        mode = self.VEGETATION_FORMULATIONS[formulation]
        if mode:
            from anuga.abstract_2d_finite_volumes.quantity import Quantity
            for name, value in (('veg_density', density), ('veg_diameter', diameter),
                                ('veg_height', height)):
                if value is None:
                    raise ValueError('set_vegetation_drag needs density, diameter and height')
                if name not in self.quantities:
                    Quantity(self, name=name, register=True)
                self.quantities[name].set_values(value, location='centroids')
            self.vegetation_Cd = float(Cd)
            self.vegetation_bed_chezy = float(bed_chezy)
        self.vegetation_mode = mode

        self._Domain_C_struct = None
        self.gpu_interface = None
        if hasattr(self, '_gpu_boundary_info_initialized'):
            del self._gpu_boundary_info_initialized

    def set_sloped_mannings_function(self, flag: bool = True) -> None:
        """Set mannings friction function to use the sloped
        wetted area.

        The flag is tested in the python wrapper
        mannings_friction_implicit
        """
        if flag:
            self.use_sloped_mannings = True
        else:
            self.use_sloped_mannings = False



    def set_use_active_set(self, flag: bool = True) -> None:
        """Enable active-set stepping: skip cells that are provably
        unchangeable this step (dry, with an all-dry 2-ring neighbourhood).

        Only the GPU/unified compute path honours the flag (rk2/DE1,
        ader2/DE_ader2 and euler/DE0 stepping; DE2/rk3 runs full steps with a
        notice).  It is also ignored, with a notice, on domains with passive
        tracers or sediment classes, riverwalls, or sloped Manning friction:
        the active path uses the single-solve scatter flux kernel, which
        carries none of those.  It also switches the flux kernel to single-solve scatter
        mode, whose results differ from the default cell-based kernel only at
        floating-point roundoff.  Serial domains only for now -- under MPI
        the flag is ignored with a warning.

        On mostly-dry flood domains (dam break, levee breach, storm surge)
        this is the largest single speedup measured for the GPU path.  Water
        added by rate operators (rainfall) is never lost -- newly wetted
        cells activate on the next step -- but widespread rain activates the
        whole mesh and the speedup degrades toward zero, hence the warning
        when both are in play.
        """
        self.use_active_set = bool(flag)
        gi = getattr(self, 'gpu_interface', None)
        if gi is not None and getattr(gi, 'gpu_dom', None) is not None:
            try:
                from anuga.shallow_water.sw_domain_gpu_ext import set_use_active_set_gpu
                set_use_active_set_gpu(gi.gpu_dom, 1 if flag else 0)
            except ImportError:
                pass
        if flag:
            if getattr(self, 'number_of_tracers', 0) > 0:
                log.warning('active-set stepping requested on a domain with '
                            'passive tracers / sediment classes: only the '
                            'cell-based flux kernel advects tracers, so the '
                            'GPU path will run full steps (flag ignored).')
            rate_like = [type(op).__name__ for op in
                         getattr(self, 'fractional_step_operators', [])
                         if 'rate' in type(op).__name__.lower()]
            if rate_like:
                log.warning('active-set stepping enabled with rate operators '
                            'attached (%s): results stay correct, but rained '
                            'regions activate and the speedup degrades under '
                            'widespread rain.' % ', '.join(rate_like))

    def get_use_active_set(self) -> bool:
        """Whether active-set stepping is requested (see set_use_active_set)."""
        return bool(getattr(self, 'use_active_set', False))

    def get_active_set_stats(self):
        """Return (mean_active_fraction, n_rebuilds) once evolution has run.

        mean_active_fraction is 1.0 whenever the active set never engaged
        (flag off, non-GPU compute mode, MPI fallback, or before the first
        rebuild).
        """
        gi = getattr(self, 'gpu_interface', None)
        if gi is not None and getattr(gi, 'gpu_dom', None) is not None:
            try:
                from anuga.shallow_water.sw_domain_gpu_ext import active_set_stats_gpu
                return active_set_stats_gpu(gi.gpu_dom)
            except (ImportError, AttributeError):
                pass
        return (1.0, 0)

    def set_compute_fluxes_method(self, flag: str = 'original') -> None:
        """Set method for computing fluxes.

        Currently
           original
           wb_1
           wb_2
           wb_3
           tsunami
           DE
        """
        compute_fluxes_methods = ['original', 'wb_1', 'wb_2', 'wb_3', 'tsunami', 'DE']

        if flag in compute_fluxes_methods:
            self.compute_fluxes_method = flag
        else:
            msg = 'Unknown compute_fluxes_method. \nPossible choices are:\n'+ \
            ', '.join(compute_fluxes_methods)+'.'
            raise Exception(msg)


    def get_compute_fluxes_method(self) -> str:
        """Get method for computing fluxes.

        See set_compute_fluxes_method for possible choices.
        """

        return self.compute_fluxes_method



    def set_flow_algorithm(self, algorithm: str = 'DE0') -> None:
        """Set combination of slope limiting and time stepping

        Currently
           DE0
           DE1
           DE2
           DE0_7
           DE1_7
        """

        algorithm = str(algorithm)

        # Replace any dots with dashes
        algorithm = algorithm.replace('.', '_')


        flow_algorithms = ['DE0', 'DE1', 'DE2', \
                           'DE0_7', 'DE1_7', 'DE_ader2']

        if algorithm in flow_algorithms:
            self.flow_algorithm = algorithm
        else:
            msg = 'Unknown flow_algorithm. \nPossible choices are:\n'+ \
            ', '.join(flow_algorithms)+'.'
            raise Exception(msg)

        if self.flow_algorithm == 'DE0':
            self._set_DE0_defaults()

        if self.flow_algorithm == 'DE1':
            self._set_DE1_defaults()

        if self.flow_algorithm == 'DE2':
            self._set_DE2_defaults()

        if self.flow_algorithm == 'DE0_7':
            self._set_DE0_7_defaults()

        if self.flow_algorithm == 'DE1_7':
            self._set_DE1_7_defaults()

        if self.flow_algorithm == 'DE_ader2':
            self._set_DE_ader2_defaults()


    def get_flow_algorithm(self) -> str:
        """
        Get method used for timestepping and spatial discretisation

        """

        return self.flow_algorithm


    # def set_gravity_method(self):
    #     """Gravity method is determined by the compute_fluxes_method
    #     This is now not used, as gravity is combine in the compute_fluxes method
    #     """

    #     if  self.get_compute_fluxes_method() == 'original':
    #         self.forcing_terms[0] = gravity

    #     elif self.get_compute_fluxes_method() == 'wb_1':
    #         self.forcing_terms[0] = gravity_wb

    #     elif self.get_compute_fluxes_method() == 'wb_2':
    #         self.forcing_terms[0] = gravity

    #     else:
    #         raise Exception('undefined compute_fluxes method')

    def set_extrapolate_velocity(self, flag: bool = True) -> None:
        """ Extrapolation routine uses momentum by default,
        can change to velocity extrapolation which
        seems to work better.
        """

        if flag is True:
            self.extrapolate_velocity_second_order = True
        elif flag is False:
            self.extrapolate_velocity_second_order = False

    def set_low_froude(self, low_froude: int = 0) -> None:
        """ For low Froude problems the standard flux calculations
        can lead to excessive damping. Set low_froude to 1 or 2 for
        flux calculations which minimize the damping in this case.
        """

        assert low_froude in [LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2]

        self.low_froude = low_froude

    def set_use_optimise_dry_cells(self, flag: bool = True) -> None:
        """ Try to optimize calculations where region is dry
        """

        if flag is True:
            self.optimise_dry_cells = int(True)
        elif flag is False:
            self.optimise_dry_cells = int(False)




    def set_use_kinematic_viscosity(self, flag: bool = True) -> None:

        from anuga.operators.kinematic_viscosity_operator import Kinematic_viscosity_operator

        if flag :
            # Create Operator if necessary
            if self.kv_operator is None:
                self.kv_operator = Kinematic_viscosity_operator(self)
        else:
            if self.kv_operator is None:
                return
            else:
                # Remove operator from fractional_step_operators
                self.fractional_step_operators.remove(self.kv_operator)
                self.kv_operator = None





    def set_collect_max_quantities(self,
                               update_frequency=1,
                               collection_start_time=0.,
                               velocity_zero_height=None,
                               store_to_sww=True) -> Collect_max_quantities_operator:
        """Create (or return existing) Collect_max_quantities_operator on this domain.

        Tracks running maxima of stage, depth, speed, and momentum magnitude
        (||(uh, vh)||) over the simulation.  Call once before domain.evolve().

        Parameters
        ----------
        update_frequency : int
            Update maxima every this many timesteps (default 1).
        collection_start_time : float
            Only collect after this simulation time (default 0).
        velocity_zero_height : float or None
            Zero velocity below this depth; defaults to minimum_allowed_height.
        store_to_sww : bool
            If True (default), write running maxima to the SWW file every yield
            step as centroid quantities max_stage_c, max_depth_c, max_speed_c,
            max_uh_c.

        Returns
        -------
        Collect_max_quantities_operator
        """
        from anuga.operators.collect_max_quantities_operator import \
            Collect_max_quantities_operator

        if self.max_quantities_operator is None:
            self.max_quantities_operator = Collect_max_quantities_operator(
                self,
                update_frequency=update_frequency,
                collection_start_time=collection_start_time,
                velocity_zero_height=velocity_zero_height,
                store_to_sww=store_to_sww,
            )
        return self.max_quantities_operator


    def set_beta(self, beta: float) -> None:
        """Shorthand to assign one constant value [0,2] to all limiters.
        0 Corresponds to first order, where as larger values make use of
        the second order scheme.
        """

        self.beta_w = beta
        self.beta_w_dry = beta
        self.quantities['stage'].beta = beta

        self.beta_uh = beta
        self.beta_uh_dry = beta
        self.quantities['xmomentum'].beta = beta

        self.beta_vh = beta
        self.beta_vh_dry = beta
        self.quantities['ymomentum'].beta = beta


    def set_betas(self, beta_w: float, beta_w_dry: float, beta_uh: float,
                  beta_uh_dry: float, beta_vh: float, beta_vh_dry: float) -> None:
        """Assign beta values in the range  [0,2] to all limiters.
        0 Corresponds to first order, where as larger values make use of
        the second order scheme.
        """

        self.beta_w = beta_w
        self.beta_w_dry = beta_w_dry
        self.quantities['stage'].beta = beta_w

        self.beta_uh = beta_uh
        self.beta_uh_dry = beta_uh_dry
        self.quantities['xmomentum'].beta = beta_uh

        self.beta_vh = beta_vh
        self.beta_vh_dry = beta_vh_dry
        self.quantities['ymomentum'].beta = beta_vh


    def set_store_vertices_uniquely(self, flag: bool = True, reduction: Callable | None = None) -> None:
        """Decide whether vertex values should be stored uniquely as
        computed in the model (True) or whether they should be reduced to one
        value per vertex using self.reduction (False).
        """

        # FIXME (Ole): how about using the word "continuous vertex values" or
        # "continuous stage surface"
        self.smooth = not flag

        # Reduction operation for get_vertex_values
        if reduction is None:
            self.reduction = mean
            #self.reduction = min  #Looks better near steep slopes

    def set_store_vertices_smoothly(self, flag: bool = True, reduction: Callable | None = None) -> None:
        """Decide whether vertex values should be stored smoothly (one value per vertex)
        or uniquely as
        computed in the model (False).
        """

        # FIXME (Ole): how about using the word "continuous vertex values" or
        # "continuous stage surface"
        self.smooth = flag

        # Reduction operation for get_vertex_values
        if reduction is None:
            self.reduction = mean
            #self.reduction = min  #Looks better near steep slopes

    def set_minimum_storable_height(self, minimum_storable_height: float) -> None:
        """Set the minimum depth that will be written to an SWW file.

        minimum_storable_height  minimum allowed SWW depth is in meters

        This is useful for removing thin water layers that seems to be caused
        by friction creep.
        """

        self.minimum_storable_height = minimum_storable_height


    def get_minimum_storable_height(self) -> float:

        return self.minimum_storable_height


    def set_minimum_allowed_height(self, minimum_allowed_height: float) -> None:
        """Set minimum depth that will be recognised in the numerical scheme.

        minimum_allowed_height  minimum allowed depth in meters

        The parameter H0 (Minimal height for flux computation) is also set by
        this function.
        """

        #FIXME (Ole): rename H0 to minimum_allowed_height_in_flux_computation

        #FIXME (Ole): Maybe use histogram to identify isolated extreme speeds
        #and deal with them adaptively similarly to how we used to use 1 order
        #steps to recover.

        self.minimum_allowed_height = minimum_allowed_height
        self.H0 = minimum_allowed_height



    def get_minimum_allowed_height(self) -> float:

        return self.minimum_allowed_height

    def set_negative_volume_warning_fraction(self, fraction: float) -> None:
        """Set the fraction of the total domain water volume that must be added by
        clamping negative-depth cells to zero depth in a single timestep before
        update_conserved_quantities() emits a "possible loss of conservation"
        warning.

        Clamping a few cells by a near-zero depth is normal in wetting/drying and
        involves negligible volume, so the default (see
        anuga.config.negative_volume_warning_fraction) avoids warning on almost
        every step. Set to 0.0 to warn whenever any volume is added.
        """

        self.negative_volume_warning_fraction = fraction

    def get_negative_volume_warning_fraction(self) -> float:

        return self.negative_volume_warning_fraction

    def set_maximum_allowed_speed(self, maximum_allowed_speed: float) -> None:
        """Set the maximum particle speed that is allowed in water shallower
        than minimum_allowed_height.

        maximum_allowed_speed

        This is useful for controlling speeds in very thin layers of water and
        at the same time allow some movement avoiding pooling of water.
        """

        self.maximum_allowed_speed = maximum_allowed_speed

    def set_points_file_block_line_size(self, points_file_block_line_size: int) -> None:
        """
        """

        self.points_file_block_line_size = points_file_block_line_size


    # FIXME: Probably obsolete in its curren form
    def set_quantities_to_be_stored(self, q: dict[str, int] | list[str] | None) -> None:
        """Specify which quantities will be stored in the SWW file.

        q must be either:
          - a dictionary with quantity names
          - a list of quantity names (for backwards compatibility)
          - None

        The format of the dictionary is as follows

        quantity_name: flag where flag must be either 1 or 2.
        If flag is 1, the quantity is considered static and will be
        stored once at the beginning of the simulation in a 1D array.

        If flag is 2, the quantity is considered time dependent and
        it will be stored at each yieldstep by appending it to the
        appropriate 2D array in the sww file.

        If q is None, storage will be switched off altogether.

        Once the simulation has started and thw sww file opened,
        this function will have no effect.

        The format, where q is a list of names is for backwards compatibility
        only.
        It will take the specified quantities to be time dependent and assume
        'elevation' to be static regardless.
        """

        if q is None:
            self.quantities_to_be_stored = {}
            self.store = False
            return

        # Check correctness
        for quantity_name in q:
            msg = ('Quantity %s is not a valid conserved quantity'
                   % quantity_name)
            assert quantity_name in self.quantities, msg

        assert isinstance(q, dict)
        self.quantities_to_be_stored = q

    def get_wet_elements(self, indices: list[int] | num.ndarray | None = None,
                         minimum_height: float | None = None) -> num.ndarray:
        """Return indices for elements where h > minimum_allowed_height

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            indices = get_wet_elements()

        Note, centroid values are used for this operation
        """

        # Water depth below which it is considered to be 0 in the model
        # FIXME (Ole): Allow this to be specified as a keyword argument as well
        from anuga.config import minimum_allowed_height

        if minimum_height is None:
            minimum_height = minimum_allowed_height

        elevation = self.get_quantity('elevation').\
                        get_values(location='centroids', indices=indices)
        stage = self.get_quantity('stage').\
                    get_values(location='centroids', indices=indices)
        depth = stage - elevation

        # Select indices for which depth > 0
        wet_indices = num.compress(depth > minimum_height,
                                   num.arange(len(depth)))
        return wet_indices

    def load_balance_statistics(self, minimum_height: float | None = None) -> dict:
        """Return load balance statistics for this domain (single-rank version).

        For a parallel domain use :meth:`Parallel_domain.load_balance_statistics`
        which gathers across all MPI ranks via Allgather.  This serial version
        returns a dict with length-1 arrays so the interface is identical.

        Parameters
        ----------
        minimum_height : float, optional
            Depth threshold for "wet" classification.  Defaults to
            ``anuga.config.minimum_allowed_height``.

        Returns
        -------
        dict
            Keys and shapes are the same as the parallel version::

                n_full           int[1]   total triangle count
                n_ghost          int[1]   0 (no ghost triangles in serial)
                n_wet_full       int[1]   wet triangle count
                wet_fraction     float[1] n_wet_full / n_full
                ghost_fraction   float[1] 0.0
                wall_time        float[1] total wall time since evolve() started
                comm_time        float[1] 0.0
                reduce_wait_time float[1] 0.0
                compute_time     float[1] same as wall_time
                imbalance_ratio  float    1.0
                wet_compute_corr float    nan
        """
        from time import time as walltime
        from anuga.config import minimum_allowed_height as default_mah

        if minimum_height is None:
            minimum_height = default_mah

        n_full = self.get_number_of_triangles()
        stage_c = self.get_quantity('stage').centroid_values
        elev_c  = self.get_quantity('elevation').centroid_values
        n_wet   = int(num.sum((stage_c - elev_c) > minimum_height))

        w_time = walltime() - self.evolve_start_walltime

        return {
            'n_full':           num.array([n_full], dtype=int),
            'n_ghost':          num.array([0],      dtype=int),
            'n_wet_full':       num.array([n_wet],  dtype=int),
            'wet_fraction':     num.array([n_wet / n_full if n_full > 0 else 0.0]),
            'ghost_fraction':   num.array([0.0]),
            'wall_time':        num.array([w_time]),
            'comm_time':        num.array([0.0]),
            'reduce_wait_time': num.array([0.0]),
            'compute_time':     num.array([w_time]),
            'imbalance_ratio':  1.0,
            'wet_compute_corr': float('nan'),
        }

    def print_load_balance_statistics(self, minimum_height: float | None = None) -> None:
        """Print a load balance summary to stdout.

        For a single-process domain this just reports wet fraction and
        triangle count.  The parallel override prints a per-rank table.

        Parameters
        ----------
        minimum_height : float, optional
            Passed through to :meth:`load_balance_statistics`.
        """
        stats = self.load_balance_statistics(minimum_height=minimum_height)
        n      = stats['n_full'][0]
        n_wet  = stats['n_wet_full'][0]
        wf     = stats['wet_fraction'][0]
        w_time = stats['wall_time'][0]
        print(f"Load balance: triangles={n}, wet={n_wet} ({100*wf:.1f}%), "
              f"wall_time={w_time:.3f}s")

    def get_maximum_inundation_elevation(self, indices: list[int] | num.ndarray | None = None,
                                         minimum_height: float | None = None) -> float:
        """Return highest elevation where h > 0

        Optional argument:
            indices is the set of element ids that the operation applies to.
            minimum_height for testing h > minimum_height
        Usage:
            q = get_maximum_inundation_elevation()

        Note, centroid values are used for this operation
        """

        wet_elements = self.get_wet_elements(indices, minimum_height)
        return self.get_quantity('elevation').\
                   get_maximum_value(indices=wet_elements)

    def get_maximum_inundation_location(self, indices: list[int] | num.ndarray | None = None) -> tuple[float, float]:
        """Return location of highest elevation where h > 0

        Optional argument:
            indices is the set of element ids that the operation applies to.

        Usage:
            q = get_maximum_inundation_location()

        Note, centroid values are used for this operation
        """

        wet_elements = self.get_wet_elements(indices)
        return self.get_quantity('elevation').\
                   get_maximum_location(indices=wet_elements)

    def get_global_wet_element_count(self, indices: list[int] | num.ndarray | None = None,
                                     minimum_height: float | None = None) -> int:
        """Return total number of wet elements across all MPI ranks.

        Optional arguments:
            indices: set of element ids that the operation applies to
            minimum_height: threshold for considering an element wet

        Usage:
            count = get_global_wet_element_count()

        Note: This performs MPI reduction, so all ranks get the same result.
        """
        from anuga import numprocs

        wet_indices = self.get_wet_elements(indices, minimum_height)
        local_count = len(wet_indices)

        if numprocs == 1:
            return local_count

        from mpi4py import MPI
        global_count = MPI.COMM_WORLD.allreduce(local_count, op=MPI.SUM)
        return global_count

    def get_global_max_stage(self, indices: list[int] | num.ndarray | None = None) -> float:
        """Return maximum stage value across all MPI ranks.

        Optional argument:
            indices: set of element ids that the operation applies to

        Usage:
            max_stage = get_global_max_stage()

        Note: This performs MPI reduction, so all ranks get the same result.
        """
        from anuga import numprocs

        stage = self.get_quantity('stage')
        local_max = stage.get_maximum_value(indices)

        if numprocs == 1:
            return local_max

        from mpi4py import MPI
        global_max = MPI.COMM_WORLD.allreduce(local_max, op=MPI.MAX)
        return global_max

    def get_global_max_speed(self) -> float:
        """Return maximum speed across all MPI ranks.

        Usage:
            max_speed = get_global_max_speed()

        Note: This performs MPI reduction, so all ranks get the same result.
        """
        from anuga import numprocs
        import numpy as num

        if self.max_speed is None:
            return 0.0

        local_max = num.max(self.max_speed)

        if numprocs == 1:
            return local_max

        from mpi4py import MPI
        global_max = MPI.COMM_WORLD.allreduce(local_max, op=MPI.MAX)
        return global_max

    def get_water_volume(self, region=None,
                         indices: list[int] | num.ndarray | None = None) -> float:
        """Volume of water (m^3) in the domain, or in part of it.

        Parameters
        ----------
        region : Region, optional
            Restrict the volume to the triangles of an
            :class:`~anuga.Region` (built from a polygon, a circle or
            explicit indices), as :meth:`Quantity.get_integral` does.
        indices : array-like of int, optional
            Restrict to these triangle indices instead. At most one of
            `region` and `indices` may be given.

        Returns
        -------
        float
            The volume, summed across ranks under MPI (every rank must call
            this). Only ghost-free (full) triangles are counted, so a region
            straddling a partition boundary is not double-counted.

        Notes
        -----
        The whole-domain volume is appended to ``volume_history``; a regional
        query is not, since the history tracks the domain total. Under GPU
        offload the whole-domain volume is reduced on the device; a regional
        query syncs the state to the host first.

        Examples
        --------
        >>> pond = anuga.Region(domain, polygon=[[0, 0], [10, 0], [10, 10], [0, 10]])
        >>> domain.get_water_volume(region=pond)
        """
        from anuga import numprocs

        assert region is None or indices is None, \
            'get_water_volume: give region or indices, not both'
        partial = region is not None or indices is not None

        if (self.multiprocessor_mode == MULTIPROCESSOR_GPU
                and self.gpu_interface is not None and not partial):
            # Whole domain: reduce on the device (avoids an expensive D2H sync)
            from anuga.shallow_water.sw_domain_gpu_ext import compute_water_volume_gpu
            volume = compute_water_volume_gpu(self.gpu_interface.gpu_dom)
        else:
            if (self.multiprocessor_mode == MULTIPROCESSOR_GPU
                    and self.gpu_interface is not None):
                self.gpu_interface.sync_from_device()
            volume = self._height_quantity().get_integral(region=region,
                                                          indices=indices)

        if numprocs > 1:
            # Use MPI_Allreduce instead of manual gather-broadcast
            from mpi4py import MPI
            volume = MPI.COMM_WORLD.allreduce(volume, op=MPI.SUM)

        if not partial:
            self.volume_history.append(volume)
        return volume

    def _height_quantity(self):
        """The water depth as a Quantity, from the stage and elevation
        centroids.

        Always rebuilt from stage - elevation rather than taken from the
        'height' quantity: the device sync brings stage and elevation back to
        the host but not height, so under GPU offload the stored height is
        stale between yieldsteps. The two agree wherever height is current.
        """
        from anuga import Quantity
        Stage = self.quantities['stage']
        Elev = self.quantities['elevation']
        Height = Quantity(self)
        Height.set_values(Stage.centroid_values - Elev.centroid_values,
                          location='centroids')
        return Height

    def get_boundary_flux_integral(self) -> float:
        """Compute the boundary flux integral.

        Should work in parallel
        """

        from anuga import numprocs

        if not self.compute_fluxes_method=='DE':
            msg='Boundary flux integral only supported for DE fluxes '+\
                '(because computation of boundary_flux_sum is only implemented there)'
            raise Exception(msg)

        flux_integral = self.boundary_flux_integral.boundary_flux_integral[0]

        if numprocs == 1:
            return flux_integral

        # Use MPI_Allreduce instead of manual gather-broadcast
        from mpi4py import MPI
        flux_integral = MPI.COMM_WORLD.allreduce(flux_integral, op=MPI.SUM)

        return flux_integral

    def get_fractional_step_volume_integral(self) -> float:
        """Compute the integrated flows from fractional steps.

        This requires that the fractional step operators update the fractional_step_volume_integral.

        Should work in parallel
        """

        from anuga import numprocs

        flux_integral = self.fractional_step_volume_integral

        if numprocs == 1:
            return flux_integral

        # Use MPI_Allreduce instead of manual gather-broadcast
        from mpi4py import MPI
        flux_integral = MPI.COMM_WORLD.allreduce(flux_integral, op=MPI.SUM)

        return flux_integral

    def get_flow_through_cross_section(self, polyline: ArrayLike, verbose: bool = False) -> float:
        """Get the total flow through an arbitrary poly line.

        This is a run-time equivalent of the function with same name
        in sww_interrogate.py

        Input:
            polyline: Representation of desired cross section - it may contain
                      multiple sections allowing for complex shapes. Assume
                      absolute UTM coordinates.
                      Format [[x0, y0], [x1, y1], ...]

        Output:
            Q: Total flow [m^3/s] across given segments.
        """


        cross_section = Cross_section(self, polyline, verbose)

        return cross_section.get_flow_through_cross_section()


    def get_energy_through_cross_section(self, polyline: ArrayLike,
                                         kind: str = 'total',
                                         verbose: bool = False) -> float:
        """Obtain average energy head [m] across specified cross section.

        Inputs:
            polyline: Representation of desired cross section - it may contain
                      multiple sections allowing for complex shapes. Assume
                      absolute UTM coordinates.
                      Format [[x0, y0], [x1, y1], ...]
            kind:     Select which energy to compute.
                      Options are 'specific' and 'total' (default)

        Output:
            E: Average energy [m] across given segments for all stored times.

        The average velocity is computed for each triangle intersected by
        the polyline and averaged weighted by segment lengths.

        The typical usage of this function would be to get average energy of
        flow in a channel, and the polyline would then be a cross section
        perpendicular to the flow.

        #FIXME (Ole) - need name for this energy reflecting that its dimension
        is [m].
        """



        cross_section = Cross_section(self, polyline, verbose)

        return cross_section.get_energy_through_cross_section(kind)


    def check_integrity(self) -> None:
        """ Run integrity checks on shallow water domain. """
        Generic_Domain.check_integrity(self)

        #Check that we are solving the shallow water wave equation
        msg = 'First conserved quantity must be "stage"'
        assert self.conserved_quantities[0] == 'stage', msg
        msg = 'Second conserved quantity must be "xmomentum"'
        assert self.conserved_quantities[1] == 'xmomentum', msg
        msg = 'Third conserved quantity must be "ymomentum"'
        assert self.conserved_quantities[2] == 'ymomentum', msg


    #@profile
    def compute_fluxes(self):
        """Compute fluxes and timestep suitable for all volumes in domain.

        Compute total flux for each conserved quantity using "flux_function"

        Fluxes across each edge are scaled by edgelengths and summed up
        Resulting flux is then scaled by area and stored in
        explicit_update for each of the three conserved quantities
        stage, xmomentum and ymomentum

        The maximal allowable speed computed by the flux_function for each volume
        is converted to a timestep that must not be exceeded. The minimum of
        those is computed as the next overall timestep.

        Post conditions:
          domain.explicit_update is reset to computed flux values
          domain.flux_timestep is set to the largest step satisfying all volumes.

        This wrapper calls the underlying C version of compute fluxes
        """

        # Using Gareth Davies discontinuous elevation scheme
        # Flux calculation and gravity incorporated in same
        # procedure

        nvtxRangePush("compute_fluxes")
        self._ensure_gpu_interface()
        # Choose the correct extension module
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import compute_fluxes_ext_central
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            # change over to cuda routines as developed
            # from .sw_domain_simd_ext import compute_fluxes_ext_central
            # FIXME SR: 2023_10_16 currently compute_fluxes and distribute together
            # is producing incorrect results, but work separately!
            compute_fluxes_ext_central = self.gpu_interface.compute_fluxes_ext_central_kernel
        else:
            raise Exception('Not implemented')

        timestep = self.evolve_max_timestep
        self.flux_timestep = compute_fluxes_ext_central(self, timestep)

        nvtxRangePop()

    def _init_gpu_boundary_info(self, gpu_dom) -> bool:
        """Classify this domain's boundaries for the unified (mode-2) path.

        Runs once per boundary configuration (``set_boundary`` clears the
        cache). Fills the per-type lists the device path iterates over,
        decides whether every boundary can be evaluated on the device
        (``_gpu_all_on_gpu``), and if not sets up the host<->device edge sync
        the host path needs.

        Boundary types the device cannot evaluate are reported once with a
        warning, or refused with ``RuntimeError`` in strict mode (see
        :meth:`set_compute_mode`). Returns ``_gpu_all_on_gpu``.
        """
        if getattr(self, '_gpu_boundary_info_initialized', False):
            return self._gpu_all_on_gpu

        self._gpu_cpu_tags = []
        self._gpu_cpu_boundary_types = []
        self._gpu_all_on_gpu = True
        self._gpu_transmissive_n_zero_t_boundaries = []
        self._gpu_time_boundaries = []
        self._gpu_absorbing_wave_boundaries = []
        self._gpu_characteristic_wave_boundaries = []
        self._gpu_flather_boundaries = []

        for tag, B in self.boundary_map.items():
            if B is None:
                continue
            btype = B.__class__.__name__
            if btype not in GPU_BOUNDARY_TYPES:
                self._gpu_cpu_tags.append(tag)
                self._gpu_all_on_gpu = False
                self._gpu_cpu_boundary_types.append((tag, btype))
            elif btype == 'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary':
                self._gpu_transmissive_n_zero_t_boundaries.append(B)
            elif btype == 'Time_boundary':
                self._gpu_time_boundaries.append(B)
            elif btype == 'Absorbing_wave_boundary':
                self._gpu_absorbing_wave_boundaries.append(B)
            elif btype == 'Characteristic_wave_boundary':
                self._gpu_characteristic_wave_boundaries.append(B)
            elif btype == 'Flather_external_stage_zero_velocity_boundary':
                self._gpu_flather_boundaries.append(B)

        if not self._gpu_all_on_gpu:
            msg = (
                "compute mode 'unified': these boundaries are of types the device "
                f"path cannot evaluate: {self._gpu_cpu_boundary_types}. They will be "
                "evaluated on the host every step and their edge values synced to the "
                "device, which is correct but costs a host<->device round trip per "
                f"step. Device-evaluated types: {sorted(GPU_BOUNDARY_TYPES)}.")
            if getattr(self, 'compute_mode_strict', False):
                raise RuntimeError(msg + " (strict mode: refusing the host fallback)")
            import warnings
            warnings.warn(msg, stacklevel=3)
            from .sw_domain_gpu_ext import init_boundary_edge_sync
            boundary_cell_ids = num.unique(self.boundary_cells).astype(num.intc)
            init_boundary_edge_sync(gpu_dom, boundary_cell_ids)

        self._gpu_boundary_info_initialized = True
        return self._gpu_all_on_gpu

    def update_boundary(self):
        """Go through list of boundary objects and update boundary values
        for all conserved quantities on boundary.
        It is assumed that the ordering of conserved quantities is
        consistent between the domain and the boundary object, i.e.
        the jth element of vector q must correspond to the jth conserved
        quantity in domain.
        """

        nvtxRangePush('update_boundary')

        # Time-varying tracer inflow concentrations, if any. Hooked here rather
        # than at update_boundary()'s eight call sites in generic_domain, and a
        # no-op unless set_tracer_boundary() was given a callable.
        if self.number_of_tracers > 0:
            self.update_tracer_boundary_values()

        # GPU mode - use GPU boundary functions if all boundaries are GPU-supported
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            from anuga.shallow_water.sw_domain_gpu_ext import (
                evaluate_reflective_boundary_gpu,
                evaluate_dirichlet_boundary_gpu,
                evaluate_transmissive_boundary_gpu,
                set_transmissive_n_zero_t_stage,
                evaluate_transmissive_n_zero_t_boundary_gpu,
                evaluate_time_boundary_gpu,
                set_file_boundary_values_from_domain,
                evaluate_file_boundary_gpu,
                set_absorbing_wave_value,
                evaluate_absorbing_wave_boundary_gpu,
                set_characteristic_wave_value,
                evaluate_characteristic_wave_boundary_gpu,
                set_flather_value,
                evaluate_flather_boundary_gpu,
                boundary_edge_sync,
                sync_boundary_values,
                init_boundary_edge_sync,
            )
            import numpy as np
            gpu_dom = self.gpu_interface.gpu_dom

            # Lazily initialize GPU boundary info
            self._init_gpu_boundary_info(gpu_dom)

            if self._gpu_all_on_gpu:
                # All boundaries are GPU-supported - evaluate entirely on GPU
                evaluate_reflective_boundary_gpu(gpu_dom)
                evaluate_dirichlet_boundary_gpu(gpu_dom)
                evaluate_transmissive_boundary_gpu(gpu_dom)

                # Handle transmissive_n_zero_t boundaries (need Python function call for stage)
                for B in self._gpu_transmissive_n_zero_t_boundaries:
                    stage_val = B.get_boundary_values()
                    try:
                        stage_val = float(stage_val)
                    except (TypeError, ValueError):
                        stage_val = float(stage_val[0])
                    set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
                evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)

                # Handle Time_boundary (need Python function call for values)
                self._push_gpu_time_boundary_values(gpu_dom)
                evaluate_time_boundary_gpu(gpu_dom)

                # Handle File_boundary / Field_boundary (per-edge values from SWW interpolation)
                set_file_boundary_values_from_domain(gpu_dom, self)
                evaluate_file_boundary_gpu(gpu_dom)

                # Handle Absorbing_wave_boundary (scalar wave stage updated each timestep)
                for B in self._gpu_absorbing_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        wave_val = float(value)
                    except (TypeError, ValueError):
                        wave_val = float(value[0])
                    set_absorbing_wave_value(gpu_dom, wave_val)
                evaluate_absorbing_wave_boundary_gpu(gpu_dom)

                # Handle Characteristic_wave_boundary (scalar perturbation updated each timestep)
                for B in self._gpu_characteristic_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        perturb = float(value)
                    except (TypeError, ValueError):
                        perturb = float(value[0])
                    set_characteristic_wave_value(gpu_dom, perturb)
                evaluate_characteristic_wave_boundary_gpu(gpu_dom)

                for B in self._gpu_flather_boundaries:
                    value = B.get_boundary_values()
                    try:
                        stage_val = float(value)
                    except (TypeError, ValueError):
                        stage_val = float(value[0])
                    set_flather_value(gpu_dom, stage_val)
                evaluate_flather_boundary_gpu(gpu_dom)
            else:
                # Some boundaries need CPU - sync edge values, evaluate on CPU, sync back
                boundary_edge_sync(gpu_dom)
                for tag in self.tag_boundary_cells:
                    B = self.boundary_map[tag]
                    if B is not None:
                        B.evaluate_segment(self, self.tag_boundary_cells[tag])
                sync_boundary_values(gpu_dom)
        else:
            # CPU mode - evaluate boundaries on CPU
            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]

                if B is None:
                    continue

                boundary_segment_edges = self.tag_boundary_cells[tag]

                B.evaluate_segment(self, boundary_segment_edges)

        nvtxRangePop()


    def compute_forcing_terms(self):
        """If there are any forcing functions driving the system
        they should be defined in Domain subclass and appended to
        the list self.forcing_terms
        """

        # The parameter self.flux_timestep should be updated
        # by the forcing_terms to ensure stability but it isn't
        # currently.

        nvtxRangePush('compute_forcing_terms')

        self._ensure_gpu_interface()

        if self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            # GPU mode: use GPU Manning friction, fall back to CPU for others.
            # Forcing terms may be plain functions (with __name__) or callable
            # objects which have none.
            for f in self.forcing_terms:
                if getattr(f, '__name__', None) == 'manning_friction_semi_implicit':
                    self.gpu_interface.manning_friction_kernel(self)
                else:
                    # Other forcing terms (rain, etc.) run on CPU
                    # Need to sync from device, run, sync back
                    self.gpu_interface.sync_from_device()
                    f(self)
                    self.gpu_interface.sync_to_device()
        else:
            for f in self.forcing_terms:
                f(self)

        nvtxRangePop()

    def _warn_unsupported_mode2_forcing(self):
        """Warn (once) if forcing terms other than Manning friction are present
        in mode 2.

        The mode-2 C step loop applies forcing in C and only handles Manning
        friction, so user-supplied Python ``forcing_terms`` are NOT applied —
        they are silently skipped. Use fractional-step operators instead
        (``Rate_operator``, ``Wind_stress_operator``,
        ``Barometric_pressure_operator``, or a custom ``Operator``). This
        converts that silent correctness gap into a loud, actionable message.
        """
        if getattr(self, '_warned_mode2_forcing', False):
            return
        self._warned_mode2_forcing = True
        ignored = [getattr(f, '__name__', None) or f.__class__.__name__
                   for f in self.forcing_terms
                   if getattr(f, '__name__', None) != 'manning_friction_semi_implicit']
        if ignored:
            import warnings
            warnings.warn(
                "multiprocessor_mode=2 ('unified') applies forcing in C and only "
                "handles Manning friction; these Python forcing terms are NOT "
                f"applied and are silently skipped: {ignored}. Use the equivalent "
                "operators instead — Rate_operator.rainfall()/inflow(), "
                "Wind_stress_operator, Barometric_pressure_operator.",
                stacklevel=2)

    def _warn_mode2_degenerate_protection(self):
        """Warn (once) that mode 2 does not run the degenerate-timestep protection.

        `apply_protection_against_isolated_degenerate_timesteps()` damps the
        momentum of triangles whose timestep is anomalously small. It is reached
        only from `update_timestep()`, and mode 2 ('unified') never gets there:
        the C step loops return before it, and the Python-orchestrated GPU loops
        that do call it find a host `max_speed` that the device never syncs back
        (the flux kernel writes the device copy), so the routine's own
        `max(max_speed) < 10` guard returns immediately.

        The feature is default-off (`config.protect_against_isolated_degenerate_timesteps`),
        so the sharp edge is a user who turns it on under GPU offload and gets no
        protection AND no warning. This makes it visible, as mode 2 already does
        for unsupported forcing terms.
        """
        if getattr(self, '_warned_mode2_degenerate_protection', False):
            return
        if not self.protect_against_isolated_degenerate_timesteps:
            return
        if self.multiprocessor_mode != MULTIPROCESSOR_GPU:
            return

        self._warned_mode2_degenerate_protection = True
        import warnings
        warnings.warn(
            "protect_against_isolated_degenerate_timesteps is True, but "
            "multiprocessor_mode=2 ('unified') does not implement it: no "
            "isolated-degenerate-triangle damping will be applied. Use "
            "domain.set_multiprocessor_mode(1) ('legacy') if you need this "
            "protection.",
            stacklevel=2)

    def set_boundary(self, boundary_map):
        """Associate boundary objects with tagged segments (see base class).

        Mode-2 ('unified') captures a device-side boundary classification and
        the per-edge Dirichlet/Time/File/... mappings when the GPU interface is
        first built. A later set_boundary() — e.g. switching a tag from
        Reflective to Dirichlet partway through a simulation — must invalidate
        those caches, otherwise the stale boundary keeps being applied on the
        device and the new condition is silently ignored.
        """
        super().set_boundary(boundary_map)

        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            # Preserve current conserved-quantity state before tearing the
            # interface down. In CPU-no-offload mode host arrays alias the
            # device; the sync makes this correct for the offload case too.
            try:
                self.gpu_interface.sync_from_device()
            except Exception:
                pass
            self.gpu_interface = None
            if hasattr(self, '_gpu_boundary_info_initialized'):
                del self._gpu_boundary_info_initialized
            # Rebuild immediately from the new boundary map so subsequent
            # mode-2 steps see a valid interface.
            self._ensure_gpu_interface()

    def distribute_to_vertices_and_edges(self, distribute_to_vertices=True):
        """ extrapolate centroid values to vertices and edges"""

        nvtxRangePush('distribute_to_vertices_and_edges')

        # Build a deferred mode-2 device interface on demand (a default-'unified'
        # domain may reach here, e.g. from a test, without going through evolve()).
        self._ensure_gpu_interface()

        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            # Output path (yieldsteps / direct reads). This method is NOT on the
            # mode-2 stepping path (the C RK loop extrapolates on the device
            # itself); it exists to make vertex/edge values readable from Python
            # (SWW writer, inundation queries, tests). The gpu edge kernel does
            # NOT compute vertex values, so sync centroids from the device and
            # run the host (openmp) protect + extrapolate, which produces edges
            # AND vertices. The host protect does not affect the device
            # trajectory — the next C step re-protects on the device.
            self.gpu_interface.sync_from_device()
            from .sw_domain_openmp_ext import (protect_new,
                                               extrapolate_second_order_edge_sw)
            protect_new(self)
            extrapolate_second_order_edge_sw(self, distribute_to_vertices=distribute_to_vertices)
            nvtxRangePop()
            return

        # Do protection step
        self.protect_against_infinitesimal_and_negative_heights()

        # Do extrapolation step
        # Choose the correct extension module
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import extrapolate_second_order_edge_sw
        else:
            raise Exception('Not implemented')

        extrapolate_second_order_edge_sw(self, distribute_to_vertices=distribute_to_vertices)

        nvtxRangePop()

    def distribute_to_edges(self):
        """ extrapolate centroid values edges"""

        nvtxRangePush('distribute_to_edges')

        # Do protection step
        self.protect_against_infinitesimal_and_negative_heights()

        # Do extrapolation step
        # Choose the correct extension module
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import distribute_to_edges as extrapolate_second_order_edge_sw
            extrapolate_second_order_edge_sw(self)
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            # change over to cuda routines as developed
            #from .sw_domain_simd_ext import extrapolate_second_order_edge_sw
            extrapolate_second_order_edge_sw = self.gpu_interface.extrapolate_second_order_edge_sw_kernel
            extrapolate_second_order_edge_sw(self)
        else:
            raise Exception('Not implemented')

        nvtxRangePop()

    def distribute_edges_to_vertices(self):
        """Distribute edge values to vertices.

        This is a wrapper for the C implementation of the distribution
        from edges to vertices.
        """

        nvtxRangePush('distribute_edges_to_vertices')

        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            # Using OpenMP extension
            from .sw_domain_openmp_ext import distribute_edges_to_vertices as distribute_edges_to_vertices_ext
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            # Using cupy extension
            # FIXME SR: Not implemented yet so use OpenMP version
            from .sw_domain_openmp_ext import distribute_edges_to_vertices as distribute_edges_to_vertices_ext
            # distribute_edges_to_vertices_ext = self.gpu_interface.distribute_edges_to_vertices_kernel
        else:
            raise Exception('Not implemented')


        distribute_edges_to_vertices_ext(self)
        nvtxRangePop()


    def update_timestep(self, yieldstep, finaltime):
        """Calculate the next timestep to take
        """

        # Protect against degenerate timesteps arising from isolated
        # triangles
        self.apply_protection_against_isolated_degenerate_timesteps()

        # disable variable timestepping
        if self.fixed_flux_timestep is not None:
            self.flux_timestep = self.fixed_flux_timestep
            timestep = self.fixed_flux_timestep
        else:
            # self.timestep is calculated from speed of characteristics
            # Apply CFL condition here
            timestep = min(self.CFL * self.flux_timestep, self.evolve_max_timestep)

        # Record maximal and minimal values of timestep for reporting
        self.recorded_max_timestep = max(timestep, self.recorded_max_timestep)
        self.recorded_min_timestep = min(timestep, self.recorded_min_timestep)

        # Stop if degenerate timestep
        if timestep < self.evolve_min_timestep:
            msg = 'WARNING: Too small timestep %.16f reached ' \
                % timestep
            msg += 'even after %d steps of 1 order scheme' \
                % self.max_smallsteps
            log.info(msg)
            timestep = self.evolve_min_timestep  # Try enforce min_step

            stats = self.timestepping_statistics(track_speeds=True)
            log.info(stats)

            raise Exception(msg)

        # NOTE: Now timestep is redefined. This can lead to a timestep
        #       being smaller than the self.recorded_min_timestep, which
        #       confused me (GD).
        #       The behaviour is good though, since then the
        #       recorded_min_timestep reflects the mathematical constraints on
        #       the timestep, EXCEPT the constraint that we yield at the
        #       required time. Otherwise we would often have very small
        #       recorded_min_timesteps simply because of we have to yield at a
        #       given time

        timestep = self._clip_timestep_to_output_times(timestep, yieldstep, finaltime)

        self.timestep = timestep


    def protect_against_infinitesimal_and_negative_heights(self):
        """ Clean up the stage and momentum values to ensure non-negative heights
        """

        # nvtxRangePush('protect_new')
        # Choose the correct extension module
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import protect_new
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            # change over to cuda routines as developed
            # # from .sw_domain_simd_ext import  protect_new
            #from .sw_domain_openmp_ext import protect_new
            protect_new = self.gpu_interface.protect_against_infinitesimal_and_negative_heights_kernel
        else:
            raise Exception('Not implemented')


        mass_error = protect_new(self)
        # nvtxRangePop()

        if mass_error > 0.0 and self.verbose :
            #print('Cumulative mass protection: ' + str(mass_error) + ' m^3 ')
            # From https://stackoverflow.com/questions/22397261/cant-convert-float-object-to-str-implicitly
            print('Cumulative mass protection: {0} m^3'.format(mass_error))


    def apply_protection_against_isolated_degenerate_timesteps(self):

        if self.protect_against_isolated_degenerate_timesteps is False:
            return

        # Not implemented in mode 2: max_speed is computed on the device and
        # never synced back, so the histogram below would be built from a stale
        # host array. Say so and skip, rather than silently damping nothing (or
        # damping on the strength of stale values) — see
        # _warn_mode2_degenerate_protection().
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            self._warn_mode2_degenerate_protection()
            return

        # FIXME (Ole): Make this configurable
        if num.max(self.max_speed) < 10.0:
            return

        # Setup 10 bins for speed histogram
        from anuga.utilities.numerical_tools import histogram, create_bins

        bins = create_bins(self.max_speed, 10)
        hist = histogram(self.max_speed, bins)

        # Look for characteristic signature
        if len(hist) > 1 and hist[-1] > 0 and \
           hist[4] == hist[5] == hist[6] == hist[7] == hist[8] == 0:
            # Danger of isolated degenerate triangles

            # Find triangles in last bin
            # FIXME - speed up using numeric package
            d = 0
            for i in range(self.number_of_triangles):
                if self.max_speed[i] > bins[-1]:
                    msg = 'Time=%f: Ignoring isolated high ' % self.get_time()
                    msg += 'speed triangle '
                    msg += '#%d of %d with max speed = %f' \
                        % (i, self.number_of_triangles, self.max_speed[i])

                    self.get_quantity('xmomentum').set_values(0.0, indices=[i])
                    self.get_quantity('ymomentum').set_values(0.0, indices=[i])
                    self.max_speed[i] = 0.0
                    d += 1


    def update_conserved_quantities(self):
        """Update vectors of conserved quantities using previously
        computed fluxes and specified forcing functions.
        """

        nvtxRangePush('update_conserved_quantities')

        timestep = self.timestep

        # Update height based on discontinuous elevation
        assert self.get_using_discontinuous_elevation()

        # Build (or defer to legacy) the mode-2 device interface, matching the
        # other mode-2 entry points. A direct call outside evolve() on a
        # default-'unified' domain would otherwise hit a None gpu_interface.
        self._ensure_gpu_interface()

        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import update_conserved_quantities
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            update_conserved_quantities = self.gpu_interface.update_conserved_quantities_kernel
        else:
            raise Exception('Not implemented')

        num_negative_ids, negative_volume = update_conserved_quantities(self, timestep)

        # Clamping negative-depth cells to zero depth adds water (a conservation
        # error). A few cells clamped by a near-zero depth is normal in
        # wetting/drying and involves negligible volume, so warn on the *volume*
        # added rather than the cell count: only when it is a large enough fraction
        # of the total water volume (threshold via
        # set_negative_volume_warning_fraction; 0.0 warns on any added volume).
        # The absolute floor rejects pure floating-point noise.
        #
        # SERIAL ONLY. "Loss of conservation" is a GLOBAL property, so the ratio
        # must use the whole-domain volume. In parallel this rank sees only its
        # partition: a nearly-dry sub-domain holds femtolitre-scale noise, so a
        # local ratio warns spuriously. Getting the global volume would need a
        # per-substep collective inside this hot function, which is not viable —
        # update_conserved_quantities is not called in guaranteed lock-step
        # across ranks (structure operators, euler vs rk2, small/empty
        # partitions), so any collective here deadlocks (two separate hangs were
        # traced to exactly this). In parallel, use the periodic global
        # report_water_volume_statistics() (e.g. the TOML runner's per-yieldstep
        # water balance) to check conservation instead.
        from anuga import numprocs
        if numprocs == 1 and num_negative_ids > 0 \
                and negative_volume > _negative_volume_noise_floor:
            total_volume = self.get_water_volume()
            if total_volume > 0.0 and \
                    negative_volume > self.negative_volume_warning_fraction * total_volume:
                import warnings
                fraction = negative_volume / total_volume
                msg = (
                    f'{num_negative_ids} negative cells set to zero depth, adding '
                    f'{negative_volume:.3g} m^3 ({100.0 * fraction:.3g}% of the '
                    f'{total_volume:.3g} m^3 in the domain): possible loss of '
                    'conservation. \nConsider using '
                    'domain.report_water_volume_statistics() to check the extent '
                    'of the problem'
                )
                warnings.warn(msg)

        nvtxRangePop()

    def update_other_quantities(self):
        """ There may be a need to calculate some of the other quantities
        based on the new values of conserved quantities
        """

        return


    def update_centroids_of_velocities_and_height(self):
        """Calculate the centroid values of velocities and height based
        on the values of the quantities stage and x and y momentum

        Assumes that stage and momentum are up to date

        Useful for kinematic viscosity calculations
        """

        # For shallow water we need to update height xvelocity and yvelocity

        #Shortcuts
        W  = self.quantities['stage']
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']

        #print num.min(W.centroid_values)

        # Make sure boundary values of conserved quantites
        # are consistent with value of functions at centroids
        Z.set_boundary_values_from_edges()

        # FIXME(Ole): Why are these not necessary?
        #W.set_boundary_values_from_edges()
        #UH.set_boundary_values_from_edges()
        #VH.set_boundary_values_from_edges()


        #Aliases
        w_C   = W.centroid_values
        z_C   = Z.centroid_values
        uh_C  = UH.centroid_values
        vh_C  = VH.centroid_values
        u_C   = U.centroid_values
        v_C   = V.centroid_values
        h_C   = H.centroid_values

        w_B   = W.boundary_values
        z_B   = Z.boundary_values
        uh_B  = UH.boundary_values
        vh_B  = VH.boundary_values
        u_B   = U.boundary_values
        v_B   = V.boundary_values
        h_B   = H.boundary_values

        h_C[:] = w_C-z_C
        h_C[:] = num.where(h_C >= 0, h_C , 0.0)

        h_B[:] = w_B-z_B
        h_B[:] = num.where(h_B >=0, h_B, 0.0)

        # Update height values
        #H.set_values( num.where(W.centroid_values-Z.centroid_values>=0,
        #                        W.centroid_values-Z.centroid_values, 0.0), location='centroids')
        #H.set_boundary_values( num.where(W.boundary_values-Z.boundary_values>=0,
        #                                 W.boundary_values-Z.boundary_values, 0.0))



        #assert num.min(h_C) >= 0
        #assert num.min(h_B) >= 0


        H0 = 1.0e-8

        #U.set_values(uh_C/(h_C + H0/h_C), location='centroids')
        #V.set_values(vh_C/(h_C + H0/h_C), location='centroids')

        factor = h_C/(h_C*h_C + H0)
        u_C[:]  = uh_C*factor
        v_C[:]  = vh_C*factor

        #U.set_boundary_values(uh_B/(h_B + H0/h_B))
        #V.set_boundary_values(vh_B/(h_B + H0/h_B))

        factor = h_B/(h_B*h_B + H0)
        u_B[:]  = uh_B*factor
        v_B[:]  = vh_B*factor



    def update_centroids_of_momentum_from_velocity(self):
        """
        Calculate the centroid value of x and y momentum from height and velocities.

        This method computes the centroid values of x and y momentum (xmomentum and
        ymomentum) by multiplying the centroid velocities by the centroid height values.
        The method assumes that the centroids of height and velocities are already
        up to date.

        This is particularly useful for kinematic viscosity calculations where momentum
        values at cell centroids are required.

        The method updates:
        - xmomentum.centroid_values: product of xvelocity and height at centroids
        - ymomentum.centroid_values: product of yvelocity and height at centroids

        After updating centroid values, the method distributes these values to vertices
        and edges via distribute_to_vertices_and_edges().

        Notes
        -----
        This method modifies the centroid_values arrays in-place for both xmomentum
        and ymomentum quantities.

        See Also
        --------
        distribute_to_vertices_and_edges : Distribute centroid values to vertices and edges
        """

        # For shallow water we need to update height xvelocity and yvelocity

        #Shortcuts
        UH = self.quantities['xmomentum']
        VH = self.quantities['ymomentum']
        H  = self.quantities['height']
        Z  = self.quantities['elevation']
        U  = self.quantities['xvelocity']
        V  = self.quantities['yvelocity']


        #Arrays
        u_C  = U.centroid_values
        v_C  = V.centroid_values
        uh_C = UH.centroid_values
        vh_C = VH.centroid_values
        h_C  = H.centroid_values

        u_B  = U.boundary_values
        v_B  = V.boundary_values
        uh_B = UH.boundary_values
        vh_B = VH.boundary_values
        h_B  = H.boundary_values

        uh_C[:] = u_C*h_C
        vh_C[:] = v_C*h_C
        #UH.set_values(u_C*h_C , location='centroids')
        #VH.set_values(v_C*h_C , location='centroids')

        self.distribute_to_vertices_and_edges()




    def evolve(self,
               yieldstep: float | None = None,
               outputstep: float | None = None,
               finaltime: float | DateTime | None = None,
               duration: float | None = None,
               skip_initial_step: bool = False) -> Iterator[float]:
        """Evolve method from Domain class.

        Parameters
        ----------
        yieldstep : float, optional
            Yield every yieldstep time period
        outputstep : float, optional
            Output to sww file every outputstep time period. outputstep should be
            an integer multiple of yieldstep.
        finaltime : float or datetime, optional
            Evolve until finaltime (can be a float in seconds or a datetime object)
        duration : float, optional
            Evolve for a time of length duration (seconds)
        skip_initial_step : bool, optional
            Can be used to restart a simulation (not often used).

        Notes
        -----
        If outputstep is None, the output to sww file happens every yieldstep.
        If yieldstep is None then simply evolve to finaltime or for a duration.
        """

        # Call check integrity here rather than from user scripts
        # self.check_integrity()

        from time import time as walltime
        self.evolve_start_walltime = walltime()
        self.last_walltime = self.evolve_start_walltime

        from datetime import datetime
        if finaltime is not None:
            if isinstance(finaltime, datetime):
                if finaltime.tzinfo is None:
                    dt = finaltime.replace(tzinfo=self.timezone)
                else:
                    dt = finaltime
                finaltime = dt.timestamp()
            else:
                finaltime = float(finaltime)

        if outputstep is None:
            outputstep = yieldstep

        if yieldstep is None:
            self.output_frequency = 1
        else:
            msg = f'outputstep ({outputstep}) should be an integer multiple of yieldstep ({yieldstep})'
            output_frequency = outputstep/yieldstep
            assert float(output_frequency).is_integer(), msg
            self.output_frequency = int(output_frequency)

        msg = 'Attribute self.beta_w must be in the interval [0, 2]'
        assert 0 <= self.beta_w <= 2.0, msg

        # [D-4]: the near-bed tracers are registered once every fraction
        # exists and before the C struct and the device mapping are built.
        self._ensure_nearbed_tracers()

        # Build the mode-2 device interface lazily if it was deferred (a
        # default-'unified' domain constructed before boundaries were set). Must
        # happen before distribute_to_vertices_and_edges(), which uses the
        # interface in mode 2.
        self._ensure_gpu_interface()

        # Initial update of vertex and edge values before any STORAGE
        # and or visualisation.
        # This is done again in the initialisation of the Generic_Domain
        # evolve loop but we do it here to ensure the values are ok for storage.

        self.distribute_to_vertices_and_edges()

        if self.store is True and (self.get_relative_time() == 0.0 or self.evolved_called is False):
            self.initialise_storage()

        # Eagerly run the mode-2 fractional-step setup (Boyd culvert registration
        # and CPU-only-operator detection) so its log lines print before the first
        # yielded step instead of after the t=0 yield. Cached, so the first
        # apply_fractional_steps() does not repeat it.
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            self._has_cpu_only_fractional_operators()
            self._warn_unsupported_mode2_forcing()
            self._warn_mode2_degenerate_protection()

        # Any riverwall change made before evolve() (or between two evolve()
        # calls) reaches the device here, before the first step.
        self._sync_riverwall_to_device()

        #nvtx marker
        nvtxRangePush('_evolve_base')

        # Call basic machinery from parent class
        for t in self._evolve_base(yieldstep=yieldstep,
                                   finaltime=finaltime, duration=duration,
                                   skip_initial_step=skip_initial_step):


            walltime = time.time()

            #print t , self.get_time()
            # Store model data, e.g. for subsequent visualisation
            if self.store:
                if self.yieldstep_counter%self.output_frequency == 0:
                    self.store_timestep()

            if self.checkpoint:
                save_checkpoint=False
                if self.checkpoint_step == 0:
                    if rank() == 0:
                        if walltime - self.walltime_prev > self.checkpoint_time:

                            save_checkpoint = True
                        for cpu in range(size()):
                            if cpu != rank():
                                send(save_checkpoint, cpu)
                    else:
                        save_checkpoint = receive(0)

                elif self.yieldstep_counter%self.checkpoint_step == 0:
                        save_checkpoint = True

                if save_checkpoint:
                    pickle_name = os.path.join(self.checkpoint_dir,self.get_name())+'_'+str(self.get_time())+'.pickle'
                    pickle.dump(self, open(pickle_name, 'wb'))

                    barrier()
                    self.walltime_prev = time.time()

                    #print 'Stored Checkpoint File '+pickle_name

            # Pass control on to outer loop for more specific actions
            yield(t)

            # The outer loop may have operated a riverwall (set_elevation() and
            # friends). In mode 2 the device copy is stale until pushed; do it
            # here so the change takes effect from the next timestep on.
            self._sync_riverwall_to_device()

            self.yieldstep_counter += 1

        #nvtx marker
        nvtxRangePop()


    def initialise_storage(self) -> None:
        """Create and initialise self.writer object for storing data.
        Also, save x,y and bed elevation
        """

        nvtxRangePush('SWW_file')

        # Erosion operators promote elevation to time-varying storage when they
        # are created. If it has since been reset to static (flag != 2), the
        # eroded bed will not be recorded — warn. Skip when the user
        # deliberately chose static (they were already told; e.g. the TOML
        # scenario warns at parse time).
        if getattr(self, '_erosion_present', False) \
                and self.quantities_to_be_stored.get('elevation') != 2 \
                and not getattr(self, '_elevation_static_by_user', False):
            import warnings
            warnings.warn(
                'An erosion operator is active but elevation is stored '
                'statically, so the eroded bed will not appear in the SWW '
                "output. Set domain.quantities_to_be_stored['elevation'] = 2 "
                'to store it time-varying.',
                stacklevel=2)

        # Initialise writer
        self.writer = SWW_file(self)

        # Store vertices and connectivity
        self.writer.store_connectivity()

        nvtxRangePop()


    def store_timestep(self) -> None:
        """Store time dependent quantities and time.

        Precondition:
           self.writer has been initialised
        """

        nvtxRangePush('store_timestep')
        self.writer.store_timestep()
        nvtxRangePop()


    def sww_merge(self, *args, **kwargs) -> None:
        """Dummy function for sequential algorithms where the sww produced is the final products.

        For parallel runs, a similarly named routine in parallel_shallow_water will merge all the
        sub domain sww files into a global sww file

        :param bool verbose: Flag to produce more output
        :param bool delete_old: Flag to delete sub domain sww files after
            creating global sww file

        """

        pass


    # Boundary types whose values are produced by a Python callback each step
    # (Time/File/Field, transmissive-set-stage, and the wave/Flather boundaries).
    # The single-call C RK loop (_evolve_one_rk*_step_c) sets these on the device
    # once per step, so with a multi-substep method (RK2/RK3) they are NOT
    # refreshed between substeps — unlike the legacy (mode-1) solver, which calls
    # update_boundary() before every substep. For time-varying boundaries this
    # gives an O(dt) boundary-forcing error (see issue #170). Until the
    # C RK loop evaluates them per substep, domains using such boundaries are
    # routed through the Python-orchestrated GPU loop, which refreshes them each
    # substep and so bit-matches mode-1 (at negligible GPU cost, ~<=4%).
    _PYTHON_EVALUATED_GPU_BOUNDARY_TYPES = frozenset((
        'Time_boundary', 'File_boundary', 'Field_boundary',
        'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary',
        'Absorbing_wave_boundary', 'Characteristic_wave_boundary',
        'Flather_external_stage_zero_velocity_boundary',
    ))

    def _has_python_evaluated_gpu_boundaries(self):
        """True if any boundary's value is set from a Python callback each step.

        Single-substep methods (Euler, ADER2) are unaffected — they impose the
        boundary once per step in both paths — so this is only consulted by the
        multi-substep RK2/RK3 dispatch.
        """
        bmap = getattr(self, 'boundary_map', None) or {}
        return any(
            B is not None
            and B.__class__.__name__ in self._PYTHON_EVALUATED_GPU_BOUNDARY_TYPES
            for B in bmap.values()
        )



    def _push_gpu_time_boundary_values(self, gpu_dom):
        """Push per-edge Time_boundary values to the device (mode 2).

        Each Time_boundary object carries one spatially-uniform [stage, xmom,
        ymom] from its time function. The edges of all Time_boundary tags are
        concatenated in boundary_map order (matching init_time_boundary in the
        GPU extension), so a single per-edge array addresses every time-boundary
        edge and multiple Time_boundary objects with different values no longer
        clobber one another (previously a single global scalar was shared).
        """
        import numpy as num
        from anuga.shallow_water.sw_domain_gpu_ext import set_time_boundary_values
        stage_vals = []
        xmom_vals = []
        ymom_vals = []
        for tag, B in self.boundary_map.items():
            if B is not None and B.__class__.__name__ == 'Time_boundary':
                edges = self.tag_boundary_cells.get(tag, None)
                if edges is None or len(edges) == 0:
                    continue
                ne = len(edges)
                q = B.get_boundary_values()
                stage_vals.extend([float(q[0])] * ne)
                xmom_vals.extend([float(q[1])] * ne)
                ymom_vals.extend([float(q[2])] * ne)
        if stage_vals:
            set_time_boundary_values(
                gpu_dom,
                num.ascontiguousarray(stage_vals, dtype=float),
                num.ascontiguousarray(xmom_vals, dtype=float),
                num.ascontiguousarray(ymom_vals, dtype=float))
    def evolve_one_euler_step(self, yieldstep, finaltime):
        """One Euler Time Step
        Q^{n+1} = E(h) Q^n

        Does not assume that centroid values have been extrapolated to
        vertices and edges
        """

        # GPU mode: dispatch to C Euler loop (eliminates Python round-trip overhead)
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            self._evolve_one_euler_step_c(yieldstep, finaltime)
            return

        # From centroid values calculate edge
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)

        # Apply boundary conditions
        self.update_boundary()

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        self.update_conserved_quantities()


    def evolve_one_rk2_step(self, yieldstep, finaltime):
        """One 2nd order RK timestep
        Q^{n+1} = 0.5 Q^n + 0.5 E(h)^2 Q^n

        Does not assume that centroid values have been extrapolated to
        vertices and edges
        """

        # GPU mode: use C RK loop (faster) or Python-orchestrated GPU loop.
        # Fall back to the Python-orchestrated loop when a Python-evaluated
        # (possibly time-varying) boundary is present, so it is refreshed every
        # substep and matches mode-1 (the C RK loop only sets it once per step).
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            if self.use_c_rk_loop and not self._has_python_evaluated_gpu_boundaries():
                self._evolve_one_rk2_step_c(yieldstep, finaltime)
            else:
                self._evolve_one_rk2_step_gpu(yieldstep, finaltime)
            return

        # Fractional-step operators (applied by the evolve loop *after* this step,
        # before it advances relative_time to t+dt) must see the pre-step time t —
        # consistent with DE0/DE2/DE_ader2 and the mode-2 GPU loops. The mid-step
        # set_relative_time() below advances time to t+dt for the substep-2
        # boundary evaluation, so capture t here and restore it at the end;
        # otherwise time-varying operators evaluate forcing "one step too far".
        initial_relative_time = self.get_relative_time()

        # Save initial initial conserved quantities values
        self.backup_conserved_quantities() # has C, ported to GPU

        #==========================================
        # First euler step
        #==========================================

        # From centroid values calculate edge values
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False) # has C, ported to GPU

        # Apply boundary conditions
        self.update_boundary() # has C, ported to GPU

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime) #needs C

        # Update centroid values of conserved quantities
        self.update_conserved_quantities() # has C, ported to GPU

        #===========================
        # End of first euler step
        #===========================

        # Update time
        self.set_relative_time(self.get_relative_time() + self.timestep) # needs C

        # Update ghosts
        if self.ghost_layer_width < 4:
            self.update_ghosts() # needs C

        #=========================================
        # Second Euler step using the same timestep
        # calculated in the first step. Might lead to
        # stability problems but we have not seen any
        # example.
        #=========================================

        # Update edge values
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)

        # Update boundary values
        self.update_boundary()

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update conserved quantities using timestep from first step
        self.update_conserved_quantities()

        #========================================
        # Combine initial and final values
        # of conserved quantities and cleanup
        #========================================

        # Combine steps
        self.saxpy_conserved_quantities(0.5, 0.5) # has C, not ported

        # Restore the pre-step time so fractional-step operators evaluate forcing
        # at t (not t+dt); the evolve loop advances relative_time to t+dt after
        # apply_fractional_steps(). Fixes an operator-timing mismatch where DE1
        # evaluated time-varying operators one step too far, unlike DE0/DE2 and
        # the mode-2 GPU loops.
        self.set_relative_time(initial_relative_time)

    def evolve_one_ader2_step(self, yieldstep, finaltime):
        """One ADER-2 timestep using the local Cauchy-Kovalewski predictor.

        Q^{n+1} = Q^n + dt * R(Q^{n+1/2})

        Uses the fused edge predictor: edge values are shifted to Q^{n+1/2}
        in-place while centroid values remain at Q^n, eliminating the second
        extrapolation pass and the backup/saxpy restore pattern.

        Single-flux-call variant: the previous step's CFL timestep is reused
        for the predictor half-advance.  The first step bootstraps with dt=0
        (Euler) to establish the initial CFL timestep.

        Cost: 1 flux call + 1 extrapolation + 1 edge C-K predictor (after step 1).
        Accuracy: 2nd-order in space and time.
        """
        # GPU mode: use C loop (faster) or Python-orchestrated GPU loop
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            if self.use_c_rk_loop:
                self._evolve_one_ader2_step_c(yieldstep, finaltime)
            else:
                self._evolve_one_ader2_step_gpu(yieldstep, finaltime)
            return

        from .sw_domain_openmp_ext import ader_ck_predictor_edge

        # Bootstrap: prev_dt=0 on first call → Euler step to establish CFL dt
        if not hasattr(self, '_ader2_prev_dt'):
            self._ader2_prev_dt = 0.0

        prev_dt = self._ader2_prev_dt

        # Extrapolate Q^n centroids → edges (1 pass; centroids untouched)
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)
        self.update_boundary()

        if prev_dt > 0.0:
            # Fused edge predictor: shift edge values to Q^{n+1/2}, Q^n in centroids
            ader_ck_predictor_edge(self, prev_dt * 0.5)
            # Re-apply boundary conditions to boundary edges
            self.update_boundary()

        # Single flux call from Q^{n+1/2} edge values (or Q^n on bootstrap step)
        self.compute_fluxes()
        self.compute_forcing_terms()

        # Clip to yieldstep / finaltime / evolve_max_timestep.
        # In MPI parallel mode update_timestep does an Allreduce so that
        # self.timestep is the global-minimum CFL dt across all ranks.
        # _ader2_prev_dt must equal the timestep actually taken so that
        # all ranks use the same predictor half-step next iteration;
        # saving the local flux_timestep before the Allreduce would give
        # each rank a different prev_dt, corrupting ghost-boundary edges.
        self.update_timestep(yieldstep, finaltime)

        # Record the timestep actually taken as prev_dt for the next predictor.
        self._ader2_prev_dt = self.timestep

        # Q^n centroids are untouched — update directly (no saxpy restore needed)
        self.update_conserved_quantities()           # Q^{n+1} = Q^n + dt*R(Q^{n+1/2})

    def _evolve_one_ader2_step_gpu(self, yieldstep, finaltime):
        """Python-orchestrated GPU implementation of ADER-2 step.

        Fallback when the C ADER-2 loop cannot be used (e.g., unsupported
        boundary types). Prefer _evolve_one_ader2_step_c() for better performance.

        Uses the fused edge predictor with _ader2_prev_dt to match the CPU path
        exactly (single flux call, no backup/restore of centroid values needed).
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            extrapolate_second_order_gpu,
            protect_gpu,
            compute_fluxes_gpu,
            update_conserved_quantities_gpu,
            ader_ck_predictor_edge_gpu,
            sync_boundary_values,
            init_boundary_edge_sync,
            boundary_edge_sync,
            exchange_ghosts,
            evaluate_reflective_boundary_gpu,
            evaluate_dirichlet_boundary_gpu,
            evaluate_transmissive_boundary_gpu,
            set_transmissive_n_zero_t_stage,
            evaluate_transmissive_n_zero_t_boundary_gpu,
            evaluate_time_boundary_gpu,
            set_file_boundary_values_from_domain,
            evaluate_file_boundary_gpu,
            set_absorbing_wave_value,
            evaluate_absorbing_wave_boundary_gpu,
            set_characteristic_wave_value,
            evaluate_characteristic_wave_boundary_gpu,
            set_flather_value,
            evaluate_flather_boundary_gpu,
        )
        import numpy as np

        gpu_dom = self.gpu_interface.gpu_dom

        self._init_gpu_boundary_info(gpu_dom)

        # Bootstrap: prev_dt=0 on first call → plain Euler step
        if not hasattr(self, '_ader2_prev_dt'):
            self._ader2_prev_dt = 0.0
        prev_dt = self._ader2_prev_dt

        def _eval_boundaries():
            if self._gpu_all_on_gpu:
                evaluate_reflective_boundary_gpu(gpu_dom)
                evaluate_dirichlet_boundary_gpu(gpu_dom)
                evaluate_transmissive_boundary_gpu(gpu_dom)
                for B in self._gpu_transmissive_n_zero_t_boundaries:
                    stage_val = B.get_boundary_values()
                    try:
                        stage_val = float(stage_val)
                    except (TypeError, ValueError):
                        stage_val = float(stage_val[0])
                    set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
                evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)
                self._push_gpu_time_boundary_values(gpu_dom)
                evaluate_time_boundary_gpu(gpu_dom)
                set_file_boundary_values_from_domain(gpu_dom, self)
                evaluate_file_boundary_gpu(gpu_dom)
                for B in self._gpu_absorbing_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        wave_val = float(value)
                    except (TypeError, ValueError):
                        wave_val = float(value[0])
                    set_absorbing_wave_value(gpu_dom, wave_val)
                evaluate_absorbing_wave_boundary_gpu(gpu_dom)
                for B in self._gpu_characteristic_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        perturb = float(value)
                    except (TypeError, ValueError):
                        perturb = float(value[0])
                    set_characteristic_wave_value(gpu_dom, perturb)
                evaluate_characteristic_wave_boundary_gpu(gpu_dom)

                for B in self._gpu_flather_boundaries:
                    value = B.get_boundary_values()
                    try:
                        stage_val = float(value)
                    except (TypeError, ValueError):
                        stage_val = float(value[0])
                    set_flather_value(gpu_dom, stage_val)
                evaluate_flather_boundary_gpu(gpu_dom)
            else:
                boundary_edge_sync(gpu_dom)
                for tag in self.tag_boundary_cells:
                    B = self.boundary_map[tag]
                    if B is not None:
                        B.evaluate_segment(self, self.tag_boundary_cells[tag])
                sync_boundary_values(gpu_dom)

        # --- Step 1: extrapolate Q^n → edges + evaluate boundaries ---
        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        _eval_boundaries()

        if prev_dt > 0.0:
            # --- Step 2: fused edge C-K predictor → Q^{n+1/2} edges in-place ---
            # (centroid values unchanged — no backup/restore needed)
            ader_ck_predictor_edge_gpu(gpu_dom, prev_dt * 0.5)
            # Re-apply boundaries to boundary edges
            _eval_boundaries()

        # --- Step 3: single flux call from Q^{n+1/2} edges ---
        self.flux_timestep = compute_fluxes_gpu(gpu_dom)
        self.compute_forcing_terms()

        # --- Step 4: Allreduce + yieldstep/finaltime clip ---
        self.update_timestep(yieldstep, finaltime)
        self._ader2_prev_dt = self.timestep

        # --- Step 5: update Q^{n+1} = Q^n + timestep * R (no restore needed) ---
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        self.set_relative_time(self.get_relative_time() + self.timestep)
        # Record the CFL-constrained step (pre yield/final cap), matching legacy
        # update_timestep(), rather than the yield-limited step actually taken.
        cfl_dt = min(self.CFL * self.flux_timestep, self.evolve_max_timestep)
        self.recorded_max_timestep = max(cfl_dt, self.recorded_max_timestep)
        self.recorded_min_timestep = min(cfl_dt, self.recorded_min_timestep)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)

    def _evolve_one_ader2_step_c(self, yieldstep, finaltime):
        """ADER-2 step executed entirely in C - eliminates Python round-trip overhead.

        Prefer this over _evolve_one_ader2_step_gpu() for better performance.
        Falls back to _evolve_one_ader2_step_gpu() if any boundary requires CPU.
        Rate_operators must be applied separately (after this call).
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            evolve_one_ader2_step_gpu,
            set_transmissive_n_zero_t_stage,
            set_file_boundary_values_from_domain,
            set_absorbing_wave_value,
            set_characteristic_wave_value,
            set_flather_value,
        )

        gpu_dom = self.gpu_interface.gpu_dom

        self._init_gpu_boundary_info(gpu_dom)

        if not self._gpu_all_on_gpu:
            return self._evolve_one_ader2_step_gpu(yieldstep, finaltime)

        for B in self._gpu_transmissive_n_zero_t_boundaries:
            stage_val = B.get_boundary_values()
            try:
                stage_val = float(stage_val)
            except (TypeError, ValueError):
                stage_val = float(stage_val[0])
            set_transmissive_n_zero_t_stage(gpu_dom, stage_val)

        self._push_gpu_time_boundary_values(gpu_dom)

        set_file_boundary_values_from_domain(gpu_dom, self)

        for B in self._gpu_absorbing_wave_boundaries:
            value = B.get_boundary_values()
            try:
                wave_val = float(value)
            except (TypeError, ValueError):
                wave_val = float(value[0])
            set_absorbing_wave_value(gpu_dom, wave_val)

        for B in self._gpu_characteristic_wave_boundaries:
            value = B.get_boundary_values()
            try:
                perturb = float(value)
            except (TypeError, ValueError):
                perturb = float(value[0])
            set_characteristic_wave_value(gpu_dom, perturb)

        for B in self._gpu_flather_boundaries:
            value = B.get_boundary_values()
            try:
                stage_val = float(value)
            except (TypeError, ValueError):
                stage_val = float(value[0])
            set_flather_value(gpu_dom, stage_val)

        # Time-varying tracer inflow concentrations. The C step evaluates the
        # hydrodynamic boundaries on the device and never calls
        # update_boundary(), so without this a callable given to
        # set_tracer_boundary() is only re-evaluated at yield points and the
        # inflow carries a concentration up to a whole yieldstep stale.
        # A no-op without callables (the common case).
        if self.number_of_tracers > 0 and self._tracer_boundary_functions:
            self.update_tracer_boundary_values()

        max_timestep = self._get_max_timestep_to_output_times(yieldstep, finaltime)


        if not hasattr(self, '_ader2_prev_dt'):
            self._ader2_prev_dt = 0.0

        self.timestep = evolve_one_ader2_step_gpu(gpu_dom, max_timestep, 1, self._ader2_prev_dt)
        self._ader2_prev_dt = self.timestep

        # Do NOT advance relative_time here — the evolve loop does it after
        # apply_fractional_steps(); see _evolve_one_euler_step_c for why.
        # Record the CFL-constrained step (pre yield/final cap), matching legacy
        # update_timestep(), rather than the yield-limited step actually taken.
        cfl_dt = gpu_dom.recorded_flux_timestep
        self.recorded_max_timestep = max(cfl_dt, self.recorded_max_timestep)
        self.recorded_min_timestep = min(cfl_dt, self.recorded_min_timestep)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            from anuga.shallow_water.sw_domain_gpu_ext import exchange_ghosts
            exchange_ghosts(gpu_dom)

    def _evolve_one_rk2_step_gpu(self, yieldstep, finaltime):
        """Python-orchestrated GPU implementation of RK2 step.

        This is a fallback for when the C RK2 loop cannot be used (e.g., unsupported
        boundary types). Prefer _evolve_one_rk2_step_c() for better performance.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            backup_conserved_quantities_gpu,
            extrapolate_second_order_gpu,
            protect_gpu,
            compute_fluxes_gpu,
            update_conserved_quantities_gpu,
            saxpy_conserved_quantities_gpu,
            sync_boundary_values,
            init_boundary_edge_sync,
            boundary_edge_sync,
            exchange_ghosts,
            evaluate_reflective_boundary_gpu,
            evaluate_dirichlet_boundary_gpu,
            evaluate_transmissive_boundary_gpu,
            set_transmissive_n_zero_t_stage,
            evaluate_transmissive_n_zero_t_boundary_gpu,
            evaluate_time_boundary_gpu,
            set_file_boundary_values_from_domain,
            evaluate_file_boundary_gpu,
            set_absorbing_wave_value,
            evaluate_absorbing_wave_boundary_gpu,
            set_characteristic_wave_value,
            evaluate_characteristic_wave_boundary_gpu,
            set_flather_value,
            evaluate_flather_boundary_gpu,
        )
        import numpy as np

        gpu_dom = self.gpu_interface.gpu_dom

        # Supported GPU boundary types
        self._init_gpu_boundary_info(gpu_dom)

        # Backup for RK2
        backup_conserved_quantities_gpu(gpu_dom)

        # ==========================================
        # First Euler step
        # ==========================================

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)

        # Evaluate boundaries
        if self._gpu_all_on_gpu:
            evaluate_reflective_boundary_gpu(gpu_dom)
            evaluate_dirichlet_boundary_gpu(gpu_dom)
            evaluate_transmissive_boundary_gpu(gpu_dom)
            for B in self._gpu_transmissive_n_zero_t_boundaries:
                stage_val = B.get_boundary_values()
                try:
                    stage_val = float(stage_val)
                except (TypeError, ValueError):
                    stage_val = float(stage_val[0])
                set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
            evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)
            self._push_gpu_time_boundary_values(gpu_dom)
            evaluate_time_boundary_gpu(gpu_dom)
            set_file_boundary_values_from_domain(gpu_dom, self)
            evaluate_file_boundary_gpu(gpu_dom)
            for B in self._gpu_absorbing_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    wave_val = float(value)
                except (TypeError, ValueError):
                    wave_val = float(value[0])
                set_absorbing_wave_value(gpu_dom, wave_val)
            evaluate_absorbing_wave_boundary_gpu(gpu_dom)
            for B in self._gpu_characteristic_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    perturb = float(value)
                except (TypeError, ValueError):
                    perturb = float(value[0])
                set_characteristic_wave_value(gpu_dom, perturb)
            evaluate_characteristic_wave_boundary_gpu(gpu_dom)

            for B in self._gpu_flather_boundaries:
                value = B.get_boundary_values()
                try:
                    stage_val = float(value)
                except (TypeError, ValueError):
                    stage_val = float(value[0])
                set_flather_value(gpu_dom, stage_val)
            evaluate_flather_boundary_gpu(gpu_dom)
        else:
            boundary_edge_sync(gpu_dom)
            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]
                if B is not None:
                    B.evaluate_segment(self, self.tag_boundary_cells[tag])
            sync_boundary_values(gpu_dom)

        # Compute fluxes
        self.flux_timestep = compute_fluxes_gpu(gpu_dom, 0, 2)

        # Forcing terms
        self.compute_forcing_terms()

        # Update timestep
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        # End of first Euler step
        self.set_relative_time(self.get_relative_time() + self.timestep)

        # Ghost exchange (MPI)
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)

        # =========================================
        # Second Euler step
        # =========================================

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)

        # Evaluate boundaries
        if self._gpu_all_on_gpu:
            evaluate_reflective_boundary_gpu(gpu_dom)
            evaluate_dirichlet_boundary_gpu(gpu_dom)
            evaluate_transmissive_boundary_gpu(gpu_dom)
            for B in self._gpu_transmissive_n_zero_t_boundaries:
                stage_val = B.get_boundary_values()
                try:
                    stage_val = float(stage_val)
                except (TypeError, ValueError):
                    stage_val = float(stage_val[0])
                set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
            evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)
            self._push_gpu_time_boundary_values(gpu_dom)
            evaluate_time_boundary_gpu(gpu_dom)
            set_file_boundary_values_from_domain(gpu_dom, self)
            evaluate_file_boundary_gpu(gpu_dom)
            for B in self._gpu_absorbing_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    wave_val = float(value)
                except (TypeError, ValueError):
                    wave_val = float(value[0])
                set_absorbing_wave_value(gpu_dom, wave_val)
            evaluate_absorbing_wave_boundary_gpu(gpu_dom)
            for B in self._gpu_characteristic_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    perturb = float(value)
                except (TypeError, ValueError):
                    perturb = float(value[0])
                set_characteristic_wave_value(gpu_dom, perturb)
            evaluate_characteristic_wave_boundary_gpu(gpu_dom)

            for B in self._gpu_flather_boundaries:
                value = B.get_boundary_values()
                try:
                    stage_val = float(value)
                except (TypeError, ValueError):
                    stage_val = float(value[0])
                set_flather_value(gpu_dom, stage_val)
            evaluate_flather_boundary_gpu(gpu_dom)
        else:
            boundary_edge_sync(gpu_dom)
            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]
                if B is not None:
                    B.evaluate_segment(self, self.tag_boundary_cells[tag])
            sync_boundary_values(gpu_dom)

        # Compute fluxes (ignore timestep from second step)
        compute_fluxes_gpu(gpu_dom, 1, 2)

        # Forcing terms
        self.compute_forcing_terms()

        # Update conserved quantities
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        # RK2 averaging
        saxpy_conserved_quantities_gpu(gpu_dom, 0.5, 0.5)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)


    def _evolve_one_euler_step_gpu(self, yieldstep, finaltime):
        """Python-orchestrated GPU Euler (DE0) step.

        Fallback for _evolve_one_euler_step_c() when the boundary map contains a
        type the C Euler loop cannot evaluate on the device (e.g.
        Transmissive_momentum_set_stage_boundary). GPU-supported boundaries are
        evaluated on the device; any others are evaluated on the host via
        evaluate_segment() and synced back. This keeps DE0 results correct
        (identical to legacy) for every boundary type — the same fallback that
        rk2/rk3/ader2 already perform. Prefer _evolve_one_euler_step_c() (faster)
        when all boundaries are GPU-supported.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            extrapolate_second_order_gpu,
            protect_gpu,
            compute_fluxes_gpu,
            update_conserved_quantities_gpu,
            sync_boundary_values,
            init_boundary_edge_sync,
            boundary_edge_sync,
            exchange_ghosts,
            evaluate_reflective_boundary_gpu,
            evaluate_dirichlet_boundary_gpu,
            evaluate_transmissive_boundary_gpu,
            set_transmissive_n_zero_t_stage,
            evaluate_transmissive_n_zero_t_boundary_gpu,
            evaluate_time_boundary_gpu,
            set_file_boundary_values_from_domain,
            evaluate_file_boundary_gpu,
            set_absorbing_wave_value,
            evaluate_absorbing_wave_boundary_gpu,
            set_characteristic_wave_value,
            evaluate_characteristic_wave_boundary_gpu,
            set_flather_value,
            evaluate_flather_boundary_gpu,
        )
        import numpy as np

        gpu_dom = self.gpu_interface.gpu_dom

        self._init_gpu_boundary_info(gpu_dom)

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)

        if self._gpu_all_on_gpu:
            evaluate_reflective_boundary_gpu(gpu_dom)
            evaluate_dirichlet_boundary_gpu(gpu_dom)
            evaluate_transmissive_boundary_gpu(gpu_dom)
            for B in self._gpu_transmissive_n_zero_t_boundaries:
                stage_val = B.get_boundary_values()
                try:
                    stage_val = float(stage_val)
                except (TypeError, ValueError):
                    stage_val = float(stage_val[0])
                set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
            evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)
            self._push_gpu_time_boundary_values(gpu_dom)
            evaluate_time_boundary_gpu(gpu_dom)
            set_file_boundary_values_from_domain(gpu_dom, self)
            evaluate_file_boundary_gpu(gpu_dom)
            for B in self._gpu_absorbing_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    wave_val = float(value)
                except (TypeError, ValueError):
                    wave_val = float(value[0])
                set_absorbing_wave_value(gpu_dom, wave_val)
            evaluate_absorbing_wave_boundary_gpu(gpu_dom)
            for B in self._gpu_characteristic_wave_boundaries:
                value = B.get_boundary_values()
                try:
                    perturb = float(value)
                except (TypeError, ValueError):
                    perturb = float(value[0])
                set_characteristic_wave_value(gpu_dom, perturb)
            evaluate_characteristic_wave_boundary_gpu(gpu_dom)
            for B in self._gpu_flather_boundaries:
                value = B.get_boundary_values()
                try:
                    stage_val = float(value)
                except (TypeError, ValueError):
                    stage_val = float(value[0])
                set_flather_value(gpu_dom, stage_val)
            evaluate_flather_boundary_gpu(gpu_dom)
        else:
            # Host evaluation of any non-GPU boundary type (e.g.
            # Transmissive_momentum_set_stage_boundary), then sync edge values
            # back to the device for the flux kernel.
            boundary_edge_sync(gpu_dom)
            for tag in self.tag_boundary_cells:
                B = self.boundary_map[tag]
                if B is not None:
                    B.evaluate_segment(self, self.tag_boundary_cells[tag])
            sync_boundary_values(gpu_dom)

        # Compute fluxes (sets flux_timestep = CFL flux step)
        self.flux_timestep = compute_fluxes_gpu(gpu_dom)

        # Forcing terms (friction)
        self.compute_forcing_terms()

        # Update timestep to fit yieldstep/finaltime (also records
        # recorded_min/max_timestep from the CFL constraint).
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        # Do NOT advance relative_time here — the evolve loop advances it after
        # apply_fractional_steps(); see _evolve_one_euler_step_c for why.

        # Post-step ghost exchange
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)


    def _evolve_one_euler_step_c(self, yieldstep, finaltime):
        """Euler step executed entirely in C - eliminates Python round-trip overhead.

        All kernel calls (extrapolation, boundaries, fluxes, forcing, update) happen
        in C with MPI timestep reduction also in C.  Only one Python->C call per step.

        Limitations:
        - Only supports GPU-evaluated boundary types
        - Rate_operators must be applied separately (after this call)
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            evolve_one_euler_step_gpu,
            set_transmissive_n_zero_t_stage,
            set_file_boundary_values_from_domain,
            set_absorbing_wave_value,
            set_characteristic_wave_value,
            set_flather_value,
            exchange_ghosts,
        )

        gpu_dom = self.gpu_interface.gpu_dom

        self._init_gpu_boundary_info(gpu_dom)

        # Boundaries the C loop cannot evaluate on the device (e.g.
        # Transmissive_momentum_set_stage_boundary) are handled by the
        # Python-orchestrated fallback, which evaluates them on the host and
        # syncs — matching rk2/rk3/ader2. Without this, those boundaries are
        # silently ignored in DE0 and results diverge from legacy.
        if not self._gpu_all_on_gpu:
            return self._evolve_one_euler_step_gpu(yieldstep, finaltime)

        # Set time-dependent boundary values before calling C function
        for B in self._gpu_transmissive_n_zero_t_boundaries:
            stage_val = B.get_boundary_values()
            try:
                stage_val = float(stage_val)
            except (TypeError, ValueError):
                stage_val = float(stage_val[0])
            set_transmissive_n_zero_t_stage(gpu_dom, stage_val)

        self._push_gpu_time_boundary_values(gpu_dom)

        set_file_boundary_values_from_domain(gpu_dom, self)

        for B in self._gpu_absorbing_wave_boundaries:
            value = B.get_boundary_values()
            try:
                wave_val = float(value)
            except (TypeError, ValueError):
                wave_val = float(value[0])
            set_absorbing_wave_value(gpu_dom, wave_val)

        for B in self._gpu_characteristic_wave_boundaries:
            value = B.get_boundary_values()
            try:
                perturb = float(value)
            except (TypeError, ValueError):
                perturb = float(value[0])
            set_characteristic_wave_value(gpu_dom, perturb)

        for B in self._gpu_flather_boundaries:
            value = B.get_boundary_values()
            try:
                stage_val = float(value)
            except (TypeError, ValueError):
                stage_val = float(value[0])
            set_flather_value(gpu_dom, stage_val)

        # Time-varying tracer inflow concentrations. The C step evaluates the
        # hydrodynamic boundaries on the device and never calls
        # update_boundary(), so without this a callable given to
        # set_tracer_boundary() is only re-evaluated at yield points and the
        # inflow carries a concentration up to a whole yieldstep stale.
        # A no-op without callables (the common case).
        if self.number_of_tracers > 0 and self._tracer_boundary_functions:
            self.update_tracer_boundary_values()

        max_timestep = self._get_max_timestep_to_output_times(yieldstep, finaltime)


        # Execute full Euler step in C (includes MPI timestep reduction)
        self.timestep = evolve_one_euler_step_gpu(gpu_dom, max_timestep, 1)

        # NOTE: do NOT advance relative_time here. The evolve loop advances it
        # (relative_time = initial_relative_time + timestep) AFTER
        # apply_fractional_steps(), so advancing it here makes time-dependent
        # operators (variable-Q inlet, time-varying rate, ...) see the time one
        # step too far — they would evaluate forcing at t+dt instead of t.

        # Record the CFL-constrained step (pre yield/final cap), matching legacy
        # update_timestep(), rather than the yield-limited step actually taken.
        cfl_dt = gpu_dom.recorded_flux_timestep
        self.recorded_max_timestep = max(cfl_dt, self.recorded_max_timestep)
        self.recorded_min_timestep = min(cfl_dt, self.recorded_min_timestep)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)


    def _evolve_one_rk2_step_c(self, yieldstep, finaltime):
        """RK2 step executed entirely in C - eliminates Python round-trip overhead.

        This is faster than _evolve_one_rk2_step_gpu() because:
        - All kernel calls happen in C without Python round-trips
        - MPI reduction for timestep happens in C
        - Only one Python->C call per RK2 step

        Limitations:
        - Only supports GPU-evaluated boundary types
        - Rate_operators must be applied separately (after this call)
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            evolve_one_rk2_step_gpu,
            set_transmissive_n_zero_t_stage,
            set_file_boundary_values_from_domain,
            set_absorbing_wave_value,
            set_characteristic_wave_value,
            set_flather_value,
        )

        gpu_dom = self.gpu_interface.gpu_dom

        # Supported GPU boundary types
        self._init_gpu_boundary_info(gpu_dom)

        # If any boundary requires CPU, fall back to Python-orchestrated loop
        if not self._gpu_all_on_gpu:
            return self._evolve_one_rk2_step_gpu(yieldstep, finaltime)

        # Set time-dependent boundary values BEFORE calling C function
        # Python function calls are cheap (~microseconds)
        for B in self._gpu_transmissive_n_zero_t_boundaries:
            stage_val = B.get_boundary_values()
            try:
                stage_val = float(stage_val)
            except (TypeError, ValueError):
                stage_val = float(stage_val[0])
            set_transmissive_n_zero_t_stage(gpu_dom, stage_val)

        self._push_gpu_time_boundary_values(gpu_dom)

        set_file_boundary_values_from_domain(gpu_dom, self)

        for B in self._gpu_absorbing_wave_boundaries:
            value = B.get_boundary_values()
            try:
                wave_val = float(value)
            except (TypeError, ValueError):
                wave_val = float(value[0])
            set_absorbing_wave_value(gpu_dom, wave_val)

        for B in self._gpu_characteristic_wave_boundaries:
            value = B.get_boundary_values()
            try:
                perturb = float(value)
            except (TypeError, ValueError):
                perturb = float(value[0])
            set_characteristic_wave_value(gpu_dom, perturb)

        for B in self._gpu_flather_boundaries:
            value = B.get_boundary_values()
            try:
                stage_val = float(value)
            except (TypeError, ValueError):
                stage_val = float(value[0])
            set_flather_value(gpu_dom, stage_val)

        # Time-varying tracer inflow concentrations. The C step evaluates the
        # hydrodynamic boundaries on the device and never calls
        # update_boundary(), so without this a callable given to
        # set_tracer_boundary() is only re-evaluated at yield points and the
        # inflow carries a concentration up to a whole yieldstep stale.
        # A no-op without callables (the common case).
        if self.number_of_tracers > 0 and self._tracer_boundary_functions:
            self.update_tracer_boundary_values()

        max_timestep = self._get_max_timestep_to_output_times(yieldstep, finaltime)


        # Execute full RK2 step in C (includes MPI timestep reduction)
        # apply_forcing=1 enables Manning friction on GPU
        self.timestep = evolve_one_rk2_step_gpu(gpu_dom, max_timestep, 1)

        # Do NOT advance relative_time here — the evolve loop does it after
        # apply_fractional_steps(); see _evolve_one_euler_step_c for why.

        # Record the CFL-constrained step (pre yield/final cap), matching legacy
        # update_timestep(), rather than the yield-limited step actually taken.
        cfl_dt = gpu_dom.recorded_flux_timestep
        self.recorded_max_timestep = max(cfl_dt, self.recorded_max_timestep)
        self.recorded_min_timestep = min(cfl_dt, self.recorded_min_timestep)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            from anuga.shallow_water.sw_domain_gpu_ext import exchange_ghosts
            exchange_ghosts(gpu_dom)


    def _evolve_one_rk3_step_gpu(self, yieldstep, finaltime):
        """Python-orchestrated GPU implementation of SSP-RK3 step.

        This is a fallback for when the C RK3 loop cannot be used (e.g., unsupported
        boundary types). Prefer _evolve_one_rk3_step_c() for better performance.
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            backup_conserved_quantities_gpu,
            extrapolate_second_order_gpu,
            protect_gpu,
            compute_fluxes_gpu,
            update_conserved_quantities_gpu,
            saxpy_conserved_quantities_gpu,
            saxpy3_conserved_quantities_gpu,
            sync_boundary_values,
            init_boundary_edge_sync,
            boundary_edge_sync,
            exchange_ghosts,
            evaluate_reflective_boundary_gpu,
            evaluate_dirichlet_boundary_gpu,
            evaluate_transmissive_boundary_gpu,
            set_transmissive_n_zero_t_stage,
            evaluate_transmissive_n_zero_t_boundary_gpu,
            evaluate_time_boundary_gpu,
            set_file_boundary_values_from_domain,
            evaluate_file_boundary_gpu,
            set_absorbing_wave_value,
            evaluate_absorbing_wave_boundary_gpu,
            set_characteristic_wave_value,
            evaluate_characteristic_wave_boundary_gpu,
            set_flather_value,
            evaluate_flather_boundary_gpu,
        )
        import numpy as np

        gpu_dom = self.gpu_interface.gpu_dom

        # Supported GPU boundary types
        self._init_gpu_boundary_info(gpu_dom)

        def _eval_boundaries():
            """Evaluate all boundary conditions on GPU (or CPU fallback)."""
            if self._gpu_all_on_gpu:
                evaluate_reflective_boundary_gpu(gpu_dom)
                evaluate_dirichlet_boundary_gpu(gpu_dom)
                evaluate_transmissive_boundary_gpu(gpu_dom)
                for B in self._gpu_transmissive_n_zero_t_boundaries:
                    stage_val = B.get_boundary_values()
                    try:
                        stage_val = float(stage_val)
                    except (TypeError, ValueError):
                        stage_val = float(stage_val[0])
                    set_transmissive_n_zero_t_stage(gpu_dom, stage_val)
                evaluate_transmissive_n_zero_t_boundary_gpu(gpu_dom)
                self._push_gpu_time_boundary_values(gpu_dom)
                evaluate_time_boundary_gpu(gpu_dom)
                set_file_boundary_values_from_domain(gpu_dom, self)
                evaluate_file_boundary_gpu(gpu_dom)
                for B in self._gpu_absorbing_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        wave_val = float(value)
                    except (TypeError, ValueError):
                        wave_val = float(value[0])
                    set_absorbing_wave_value(gpu_dom, wave_val)
                evaluate_absorbing_wave_boundary_gpu(gpu_dom)
                for B in self._gpu_characteristic_wave_boundaries:
                    value = B.get_boundary_values()
                    try:
                        perturb = float(value)
                    except (TypeError, ValueError):
                        perturb = float(value[0])
                    set_characteristic_wave_value(gpu_dom, perturb)
                evaluate_characteristic_wave_boundary_gpu(gpu_dom)

                for B in self._gpu_flather_boundaries:
                    value = B.get_boundary_values()
                    try:
                        stage_val = float(value)
                    except (TypeError, ValueError):
                        stage_val = float(value[0])
                    set_flather_value(gpu_dom, stage_val)
                evaluate_flather_boundary_gpu(gpu_dom)
            else:
                boundary_edge_sync(gpu_dom)
                for tag in self.tag_boundary_cells:
                    B = self.boundary_map[tag]
                    if B is not None:
                        B.evaluate_segment(self, self.tag_boundary_cells[tag])
                sync_boundary_values(gpu_dom)

        initial_relative_time = self.get_relative_time()

        # Backup Q^n
        backup_conserved_quantities_gpu(gpu_dom)

        # ==========================================
        # Stage 1: Q^(1) = Q^n + h*L(Q^n)
        # ==========================================

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        _eval_boundaries()

        self.flux_timestep = compute_fluxes_gpu(gpu_dom, 0, 3)

        # Forcing terms
        self.compute_forcing_terms()

        # Update timestep
        self.update_timestep(yieldstep, finaltime)

        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        self.set_relative_time(initial_relative_time + self.timestep)

        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)

        # ==========================================
        # Stage 2: Q^(2) = Q^(1) + h*L(Q^(1))
        # ==========================================

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        _eval_boundaries()

        compute_fluxes_gpu(gpu_dom, 1, 3)
        self.compute_forcing_terms()
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        # Intermediate: Q = 0.25*Q^(2) + 0.75*Q^n
        saxpy_conserved_quantities_gpu(gpu_dom, 0.25, 0.75)

        self.set_relative_time(initial_relative_time + self.timestep * 0.5)

        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)

        # ==========================================
        # Stage 3: Q^(3) = Q^(1)_mid + h*L(Q^(1)_mid)
        # ==========================================

        protect_gpu(gpu_dom)
        extrapolate_second_order_gpu(gpu_dom)
        _eval_boundaries()

        compute_fluxes_gpu(gpu_dom, 2, 3)
        self.compute_forcing_terms()
        update_conserved_quantities_gpu(gpu_dom, self.timestep)

        # Final: Q^{n+1} = (2*Q^(3) + Q^n) / 3
        saxpy3_conserved_quantities_gpu(gpu_dom, 2.0, 1.0, 3.0)

        # Restore the pre-step time so fractional-step operators evaluate forcing
        # at t (not t+dt); see evolve_one_rk3_step. The evolve loop advances
        # relative_time to t+dt after apply_fractional_steps().
        self.set_relative_time(initial_relative_time)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            exchange_ghosts(gpu_dom)


    def _evolve_one_rk3_step_c(self, yieldstep, finaltime):
        """RK3 step executed entirely in C - eliminates Python round-trip overhead.

        This is faster than _evolve_one_rk3_step_gpu() because:
        - All kernel calls happen in C without Python round-trips
        - MPI reduction for timestep happens in C
        - Only one Python->C call per RK3 step

        Limitations:
        - Only supports GPU-evaluated boundary types
        - Rate_operators must be applied separately (after this call)
        """
        from anuga.shallow_water.sw_domain_gpu_ext import (
            evolve_one_rk3_step_gpu,
            set_transmissive_n_zero_t_stage,
            set_file_boundary_values_from_domain,
            set_absorbing_wave_value,
            set_characteristic_wave_value,
            set_flather_value,
        )

        gpu_dom = self.gpu_interface.gpu_dom

        # Supported GPU boundary types
        self._init_gpu_boundary_info(gpu_dom)

        # If any boundary requires CPU, fall back to Python-orchestrated loop
        if not self._gpu_all_on_gpu:
            return self._evolve_one_rk3_step_gpu(yieldstep, finaltime)

        # Set time-dependent boundary values BEFORE calling C function
        for B in self._gpu_transmissive_n_zero_t_boundaries:
            stage_val = B.get_boundary_values()
            try:
                stage_val = float(stage_val)
            except (TypeError, ValueError):
                stage_val = float(stage_val[0])
            set_transmissive_n_zero_t_stage(gpu_dom, stage_val)

        self._push_gpu_time_boundary_values(gpu_dom)

        set_file_boundary_values_from_domain(gpu_dom, self)

        for B in self._gpu_absorbing_wave_boundaries:
            value = B.get_boundary_values()
            try:
                wave_val = float(value)
            except (TypeError, ValueError):
                wave_val = float(value[0])
            set_absorbing_wave_value(gpu_dom, wave_val)

        for B in self._gpu_characteristic_wave_boundaries:
            value = B.get_boundary_values()
            try:
                perturb = float(value)
            except (TypeError, ValueError):
                perturb = float(value[0])
            set_characteristic_wave_value(gpu_dom, perturb)

        for B in self._gpu_flather_boundaries:
            value = B.get_boundary_values()
            try:
                stage_val = float(value)
            except (TypeError, ValueError):
                stage_val = float(value[0])
            set_flather_value(gpu_dom, stage_val)

        # Time-varying tracer inflow concentrations. The C step evaluates the
        # hydrodynamic boundaries on the device and never calls
        # update_boundary(), so without this a callable given to
        # set_tracer_boundary() is only re-evaluated at yield points and the
        # inflow carries a concentration up to a whole yieldstep stale.
        # A no-op without callables (the common case).
        if self.number_of_tracers > 0 and self._tracer_boundary_functions:
            self.update_tracer_boundary_values()

        max_timestep = self._get_max_timestep_to_output_times(yieldstep, finaltime)


        # Execute full RK3 step in C (includes MPI timestep reduction)
        # apply_forcing=1 enables Manning friction on GPU
        self.timestep = evolve_one_rk3_step_gpu(gpu_dom, max_timestep, 1)

        # Do NOT advance relative_time here — the evolve loop advances it to t+dt
        # after apply_fractional_steps(), so fractional-step operators evaluate
        # forcing at the pre-step time t (matching the rk2 C loop, DE0, DE_ader2).

        # Record the CFL-constrained step (pre yield/final cap), matching legacy
        # update_timestep(), rather than the yield-limited step actually taken.
        cfl_dt = gpu_dom.recorded_flux_timestep
        self.recorded_max_timestep = max(cfl_dt, self.recorded_max_timestep)
        self.recorded_min_timestep = min(cfl_dt, self.recorded_min_timestep)

        # Post-step ghost exchange — update_ghosts() is a no-op in GPU mode
        if self.ghost_layer_width < 4:
            from anuga.shallow_water.sw_domain_gpu_ext import exchange_ghosts
            exchange_ghosts(gpu_dom)


    def evolve_one_rk3_step(self, yieldstep, finaltime):
        """One 3rd order RK timestep
        Q^(1) = 3/4 Q^n + 1/4 E(h)^2 Q^n  (at time t^n + h/2)
        Q^{n+1} = 1/3 Q^n + 2/3 E(h) Q^(1) (at time t^{n+1})

        Does not assume that centroid values have been extrapolated to
        vertices and edges
        """

        # GPU mode: use C RK loop (faster) or Python-orchestrated GPU loop.
        # Fall back to the Python-orchestrated loop when a Python-evaluated
        # (possibly time-varying) boundary is present, so it is refreshed every
        # substep and matches mode-1 (the C RK loop only sets it once per step).
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            if self.use_c_rk_loop and not self._has_python_evaluated_gpu_boundaries():
                self._evolve_one_rk3_step_c(yieldstep, finaltime)
            else:
                self._evolve_one_rk3_step_gpu(yieldstep, finaltime)
            return

        # Save initial initial conserved quantities values
        self.backup_conserved_quantities()

        initial_relative_time = self.get_relative_time()

        ######
        # First euler step
        ######

        # From centroid values calculate edge values
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)

        # Apply boundary conditions
        self.update_boundary()

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update timestep to fit yieldstep and finaltime
        self.update_timestep(yieldstep, finaltime)

        # Update conserved quantities
        self.update_conserved_quantities()

        #====================================
        # End of first euler step
        #====================================

        # Update time
        self.set_relative_time(self.relative_time+ self.timestep)

        # Update ghosts
        self.update_ghosts()

        #============================================
        # Second Euler step using the same timestep
        # calculated in the first step. Might lead to
        # stability problems but we have not seen any
        # example.
        #============================================

        # Update edge values
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)

        # Update boundary values
        self.update_boundary()

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update conserved quantities using timestep from first step
        self.update_conserved_quantities()

        #============================================
        # End of second euler step
        #============================================

        #============================================
        # Combine steps to obtain intermediate
        # solution at time t^n + 0.5 h
        #============================================

        # Combine steps
        self.saxpy_conserved_quantities(0.25, 0.75)

        # Set substep time
        self.set_relative_time(initial_relative_time + self.timestep * 0.5)

        # Update ghosts
        self.update_ghosts()

        ######
        # Third Euler step
        ######

        # Update edge values
        self.distribute_to_vertices_and_edges(distribute_to_vertices=False)

        # Update boundary values
        self.update_boundary()

        # Compute fluxes across each element edge
        # In MPI parallel mode this involves an allreduce to find global minimal timestep
        self.compute_fluxes()

        # Compute forcing terms (friction)
        self.compute_forcing_terms()

        # Update conserved quantities using timestep from first step
        self.update_conserved_quantities()

        #=======================================
        # Combine final and initial values
        # and cleanup
        #=======================================

        # self.saxpy_conserved_quantities(2.0/3.0, 1.0/3.0)
        # This caused a roundoff error that created negative water heights

        # So do this instead!
        self.saxpy_conserved_quantities(2.0, 1.0, 3.0)

        # Restore the pre-step time so fractional-step operators (applied by the
        # evolve loop before it advances relative_time to t+dt) evaluate forcing
        # at t, not t+dt — consistent with rk2 (DE1), DE0 and DE_ader2. The
        # mid-step set_relative_time() calls above advance time for the substep
        # boundary evaluations; the evolve loop sets the final t+dt after
        # apply_fractional_steps().
        self.set_relative_time(initial_relative_time)


    def backup_conserved_quantities(self):

        # Backup conserved_quantities centroid values
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            from .sw_domain_openmp_ext import backup_conserved_quantities
            backup_conserved_quantities(self)
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            self.gpu_interface.backup_conserved_quantities_kernel(self)
        else:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                Q.backup_centroid_values()

    def saxpy_conserved_quantities(self, a, b, c=None):

        # saxpy conserved_quantities centroid values with backup values
        if self.multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            if c is None:
                c = 1.0
            from .sw_domain_openmp_ext import saxpy_conserved_quantities
            saxpy_conserved_quantities(self, a, b, c)
        elif self.multiprocessor_mode == MULTIPROCESSOR_GPU:
            if c is not None:
                from anuga.shallow_water.sw_domain_gpu_ext import saxpy3_conserved_quantities_gpu
                saxpy3_conserved_quantities_gpu(self.gpu_interface.gpu_dom, a, b, c)
            else:
                self.gpu_interface.saxpy_conserved_quantities_kernel(self, a, b)
        else:
            for name in self.conserved_quantities:
                Q = self.quantities[name]
                Q.saxpy_centroid_values(a, b)
                if c is not None:
                    Q.centroid_values[:] = Q.centroid_values / c

    def _has_cpu_only_fractional_operators(self):
        """Check if any fractional step operators require CPU execution.

        Rate_operators with GPU support don't need CPU sync.
        boundary_flux_integral_operator is GPU-safe (only reads boundary_flux_sum).
        Sediment_transport_operator is GPU-safe (device-resident kernel, updates in place).
        Boyd_box_operator/Boyd_pipe_operator are GPU-safe via GPUCulvertManager.
        Inlet_operator with GPU support doesn't need CPU sync.

        Result is cached after first call since operators don't change during simulation.
        """
        # Check cache first
        if hasattr(self, '_cached_has_cpu_only_ops'):
            return self._cached_has_cpu_only_ops

        from anuga.operators.rate_operators import Rate_operator
        from anuga.operators.boundary_flux_integral_operator import boundary_flux_integral_operator
        from anuga.operators.tracer_flux_integral_operator import (
            tracer_flux_integral_operator)
        from anuga.structures.inlet_operator import Inlet_operator
        from anuga.structures.gpu_culvert_manager import GPUCulvertManager
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.operators.sediment_operator import Sediment_transport_operator

        # Initialize GPU culvert manager for Boyd operators if needed
        has_boyd_ops = any(GPUCulvertManager.is_boyd_operator(op)
                          for op in self.fractional_step_operators)
        if has_boyd_ops and self.gpu_culvert_manager is None:
            self.gpu_culvert_manager = GPUCulvertManager(self)
            self.gpu_culvert_manager.register_all()

        result = False
        cpu_only_ops = []
        for op in self.fractional_step_operators:
            op_name = op.__class__.__name__
            if isinstance(op, Rate_operator):
                # Force GPU initialization if not already done (lazy init causes race with caching)
                if hasattr(op, '_init_gpu') and not getattr(op, '_gpu_initialized', False):
                    op._init_gpu()
                # Rate_operator with GPU support doesn't need CPU sync
                if hasattr(op, '_gpu_initialized') and op._gpu_initialized:
                    continue  # GPU-accelerated, no sync needed
            elif isinstance(op, boundary_flux_integral_operator):
                # boundary_flux_integral_operator only reads boundary_flux_sum (small array)
                # and accumulates a scalar - doesn't need full centroid sync
                continue
            elif isinstance(op, tracer_flux_integral_operator):
                # Same shape as the water one: reads the small per-substep
                # tracer_boundary_flux_sum, which the kernel writes on the HOST,
                # and accumulates into a per-tracer scalar. The per-cell scratch
                # it totals lives and is zeroed on the device, so no centroid
                # sync is needed here either.
                continue
            elif isinstance(op, Inlet_operator):
                # Force GPU initialization if not already done
                if hasattr(op, '_init_gpu') and not getattr(op, '_gpu_initialized', False):
                    op._init_gpu()
                if hasattr(op, '_gpu_initialized') and op._gpu_initialized:
                    continue  # GPU-accelerated, no sync needed
            elif isinstance(op, Collect_max_quantities_operator):
                # GPU kernel reads device-resident quantities and updates device-resident max
                # arrays; no host centroid sync needed.
                if hasattr(op, '_init_gpu') and not getattr(op, '_gpu_initialized', False):
                    op._init_gpu()
                if hasattr(op, '_gpu_initialized') and op._gpu_initialized:
                    continue  # GPU-accelerated, no sync needed
            elif isinstance(op, Sediment_transport_operator):
                # The sediment kernel runs on the device in mode 2 and updates
                # the tracer and bed arrays in place, so no host sync is needed.
                continue
            elif GPUCulvertManager.is_boyd_operator(op):
                # Handled by GPUCulvertManager (local + cross-boundary via MPI in C)
                if (self.gpu_culvert_manager is not None
                        and op in self.gpu_culvert_manager.operators):
                    continue
            # All other operators need CPU sync
            cpu_only_ops.append(op_name)
            result = True

        rank = getattr(self, 'processor', 0)
        if result:
            # Only print GPU sync warning if GPU offload is actually active
            if getattr(self, 'gpu_offload_active', True):
                print(f"[Rank {rank}] WARNING: CPU-only fractional operators detected (GPU<->CPU sync every RK2 step):")
                for name in cpu_only_ops:
                    print(f"  - {name}")
        else:
            # Only print GPU-safe message if GPU offload is actually active
            if rank == 0 and getattr(self, 'gpu_offload_active', True):
                print(f"[Rank {rank}] All fractional operators are GPU-safe, no GPU<->CPU sync needed")
        import sys; sys.stdout.flush()

        self._cached_has_cpu_only_ops = result
        return result

    def apply_fractional_steps(self):
        """Override to sync GPU data before fractional step operators run.

        Boyd culvert operators are handled via GPUCulvertManager (batched, only 2 GPU
        sync points) instead of the per-operator Python loop.

        Operator ORDER matters: fractional-step operators mutate `stage` in sequence, so
        running them in a different order gives different answers.  Mode 1 runs them in
        registration order.  Mode 2 must too, or the same script silently disagrees with
        itself across compute modes — and mode-1-vs-mode-2 agreement is the oracle we use
        to validate GPU work, so it must not depend on the order a user happened to
        construct their operators in.

        The batch is therefore fired **at the position of the first culvert in the list**,
        not before the loop.  When the culverts are registered together — which is what
        every real script does (towradgi registers 22 in a row) — that is *exactly*
        registration order, and the batching is untouched, so this costs nothing.

        If culverts are interleaved with other operators, the later ones are still pulled
        forward to the first culvert's slot (the C batches all registered culverts in one
        gather/compute/scatter cycle, and there is no unbatched mode-2 culvert path).  We
        warn in that case rather than diverge silently.  See issue #192.
        """
        gpu_mode = (self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None)
        gpu_culverts_active = (gpu_mode and self.gpu_culvert_manager is not None
                               and self.gpu_culvert_manager.num_culverts > 0)
        needs_cpu_sync = False

        if gpu_mode:
            needs_cpu_sync = self._has_cpu_only_fractional_operators()
            if needs_cpu_sync:
                self.gpu_interface.sync_from_device()

            if gpu_culverts_active:
                self._warn_if_culverts_interleaved()

        # The host copy is already in sync here and gets pushed back below, so an
        # operator that writes a quantity must not trigger its own round-trip per
        # call — that would be one full host<->device transfer pair per operator
        # per timestep.
        self._gpu_host_writes_suppressed = (gpu_mode and needs_cpu_sync)
        culverts_done = False
        try:
            for operator in self.fractional_step_operators:
                if gpu_culverts_active and operator in self.gpu_culvert_manager.operators:
                    # Fire the whole batched culvert cycle in the slot of the FIRST
                    # culvert, then skip the rest — they were handled by that batch.
                    if not culverts_done:
                        self.gpu_culvert_manager.apply_all()
                        culverts_done = True
                    continue
                operator()
        finally:
            self._gpu_host_writes_suppressed = False

        # Push the host-side work of any CPU-only operator back to the device. Mode-2
        # correctness depends on this: without it a CPU-only operator (wind stress,
        # say) writes only the host arrays and the device never sees it.
        # Re-check the interface: gpu_mode was decided BEFORE the operators ran,
        # and an operator is allowed to invalidate it mid-loop. set_tracer_source
        # does exactly that -- it rebinds the source array, so the C struct and
        # the device mapping have to be rebuilt -- and an operator that calls it
        # every step (a manufactured source, say) leaves nothing here to sync to.
        # Skipping is correct rather than merely safe: a rebuilt interface is
        # populated from the host arrays, so the host-side writes are picked up
        # anyway.
        if gpu_mode and needs_cpu_sync and self.gpu_interface is not None:
            self.gpu_interface.sync_to_device()

    def _warn_if_culverts_interleaved(self):
        """Warn once if the Boyd culverts are not contiguous in the operator list.

        The mode-2 culvert batch is a single gather/compute/scatter cycle over every
        registered culvert, so it can only be fired at one point in the sequence. We fire
        it where the first culvert sits, which reproduces registration order exactly when
        the culverts are contiguous. When they are not, the later ones get pulled forward
        and mode 2 will disagree with mode 1 — say so instead of diverging in silence.
        """

        if getattr(self, '_culvert_order_checked', False):
            return
        self._culvert_order_checked = True

        culvert_ops = self.gpu_culvert_manager.operators
        positions = [i for i, op in enumerate(self.fractional_step_operators)
                     if op in culvert_ops]
        if not positions:
            return

        contiguous = (positions[-1] - positions[0] + 1) == len(positions)
        if contiguous:
            return

        interleaved = [self.fractional_step_operators[i].__class__.__name__
                       for i in range(positions[0], positions[-1] + 1)
                       if self.fractional_step_operators[i] not in culvert_ops]
        import warnings
        warnings.warn(
            'mode 2: Boyd culvert operators are not registered contiguously — '
            f'{sorted(set(interleaved))} sit between them. The GPU applies all culverts '
            'in one batch at the position of the first, so those operators will run AFTER '
            'the culverts here but BEFORE some of them in mode 1 (legacy), and the two '
            'compute modes will not agree. Register the culverts together to avoid this. '
            '(issue #192)',
            UserWarning, stacklevel=3)

    def update_ghosts(self, quantities=None):
        """Override to use GPU ghost exchange when in GPU mode."""
        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is not None:
            # GPU RK2 loop handles ghost exchange internally via exchange_ghosts(gpu_dom)
            # Don't do another exchange here - it would cause MPI message conflicts
            pass
        else:
            # Fall back to parent implementation
            super().update_ghosts(quantities)

    def timestepping_statistics(self,
                                track_speeds: bool = False,
                                triangle_id: int | None = None,
                                relative_time: bool = False,
                                time_unit: str = 'sec',
                                datetime: bool = False) -> str:
        """Return string with time stepping statistics for printing or logging

        Parameters
        ----------
        time_units : str, optional
            Time units for reporting. Options are 'sec', 'min', 'hr', 'day'.
        datetime : bool, optional
            Flag to use timestamp or datetime.
        track_speeds : bool, optional
            Optional boolean keyword that decides whether to report location of
            smallest timestep as well as a histogram and percentile report.
        relative_time : bool, optional
            Flag to report relative time instead of absolute time.
        triangle_id : int, optional
            Can be used to specify a particular triangle rather than the one with
            the largest speed.

        Returns
        -------
        str
            Formatted string with time stepping statistics.
        """

        from anuga.config import epsilon, g

        # Call basic machinery from parent class
        msg = Generic_Domain.timestepping_statistics(self,
                                                     track_speeds=track_speeds,
                                                     triangle_id=triangle_id,
                                                     relative_time=relative_time,
                                                     time_unit=time_unit,
                                                     datetime=datetime)

        if track_speeds is True and self.max_speed is not None:
            # qwidth determines the text field used for quantities
            qwidth = self.qwidth

            # Selected triangle
            k = self.k

            # Report some derived quantities at vertices, edges and centroid
            # specific to the shallow water wave equation
            z = self.quantities['elevation']
            w = self.quantities['stage']

            Vw = w.get_values(location='vertices', indices=[k])[0]
            Ew = w.get_values(location='edges', indices=[k])[0]
            Cw = w.get_values(location='centroids', indices=[k])

            Vz = z.get_values(location='vertices', indices=[k])[0]
            Ez = z.get_values(location='edges', indices=[k])[0]
            Cz = z.get_values(location='centroids', indices=[k])

            name = 'depth'
            Vh = Vw-Vz
            Eh = Ew-Ez
            Ch = Cw-Cz

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vh[0], Vh[1], Vh[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Eh[0], Eh[1], Eh[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Ch[0])

            msg += message

            uh = self.quantities['xmomentum']
            vh = self.quantities['ymomentum']

            Vuh = uh.get_values(location='vertices', indices=[k])[0]
            Euh = uh.get_values(location='edges', indices=[k])[0]
            Cuh = uh.get_values(location='centroids', indices=[k])

            Vvh = vh.get_values(location='vertices', indices=[k])[0]
            Evh = vh.get_values(location='edges', indices=[k])[0]
            Cvh = vh.get_values(location='centroids', indices=[k])

            # Speeds in each direction
            Vu = Vuh/(Vh + epsilon)
            Eu = Euh/(Eh + epsilon)
            Cu = Cuh/(Ch + epsilon)
            name = 'U'
            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Vu[0], Vu[1], Vu[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Eu[0], Eu[1], Eu[2])

            message += '    %s: centroid_value = %.4f\n' \
                 % (name.ljust(qwidth), Cu[0])

            msg += message

            Vv = Vvh/(Vh + epsilon)
            Ev = Evh/(Eh + epsilon)
            Cv = Cvh/(Ch + epsilon)
            name = 'V'
            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Vv[0], Vv[1], Vv[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n' \
                 % (name.ljust(qwidth), Ev[0], Ev[1], Ev[2])

            message += '    %s: centroid_value = %.4f\n'\
                 %(name.ljust(qwidth), Cv[0])

            msg += message

            # Froude number in each direction
            name = 'Froude (x)'
            Vfx = Vu/(num.sqrt(g*Vh + epsilon))
            Efx = Eu/(num.sqrt(g*Eh + epsilon))
            Cfx = Cu/(num.sqrt(g*Ch + epsilon))

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vfx[0], Vfx[1], Vfx[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Efx[0], Efx[1], Efx[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Cfx[0])

            msg += message

            name = 'Froude (y)'
            Vfy = Vv/(num.sqrt(g*Vh + epsilon))
            Efy = Ev/(num.sqrt(g*Eh + epsilon))
            Cfy = Cv/(num.sqrt(g*Ch + epsilon))

            message  = '    %s: vertex_values =  %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Vfy[0], Vfy[1], Vfy[2])

            message += '    %s: edge_values =    %.4f,\t %.4f,\t %.4f\n'\
                 % (name.ljust(qwidth), Efy[0], Efy[1], Efy[2])

            message += '    %s: centroid_value = %.4f\n'\
                 % (name.ljust(qwidth), Cfy[0])

            msg += message

        return msg

    def print_timestepping_statistics(self, *args, **kwargs) -> None:
        """Print time stepping statistics.

        Parameters
        ----------
        time_units : str, optional
            Time units for reporting. Options are 'sec', 'min', 'hr', 'day'.
        datetime : bool, optional
            Flag to use timestamp or datetime.
        track_speed : bool, optional
            Optional boolean keyword that decides whether to report location of
            smallest timestep as well as a histogram and percentile report.
        relative_time : bool, optional
            Flag to report relative time instead of absolute time.
        triangle_id : int, optional
            Can be used to specify a particular triangle rather than the one with
            the largest speed.
        """

        msg = self.timestepping_statistics(*args, **kwargs)

        print(msg, flush=True)


    def compute_boundary_flows(self) -> tuple[dict[str, float], float, float]:
        """Compute boundary flows at current timestep.

        Computes the total inflow and outflow across the domain boundary,
        as well as the flow across each tagged boundary segment.

        Returns
        -------
        boundary_flows : dict
            Flow rates [m^3/s] for each boundary tag
        total_boundary_inflow : float
            Total inflow across boundary [m^3/s]
        total_boundary_outflow : float
            Total outflow across boundary [m^3/s]

        Notes
        -----
        These calculations are only approximate since they don't use the
        flux calculation used in evolve. For exact computation, see
        get_boundary_flux_integral.
        """

        # Run through boundary array and compute for each segment
        # the normal momentum ((uh, vh) dot normal) times segment length.
        # Based on sign accumulate this into boundary_inflow and
        # boundary_outflow.

        # Compute flows along boundary

        uh = self.get_quantity('xmomentum').get_values(location='edges')
        vh = self.get_quantity('ymomentum').get_values(location='edges')

        # Loop through edges that lie on the boundary and calculate
        # flows
        boundary_flows = {}
        total_boundary_inflow = 0.0
        total_boundary_outflow = 0.0
        for vol_id, edge_id in self.boundary:
            # Compute normal flow across edge. Since normal vector points
            # away from triangle, a positive sign means that water
            # flows *out* from this triangle.

            momentum = [uh[vol_id, edge_id], vh[vol_id, edge_id]]
            normal = self.mesh.get_normal(vol_id, edge_id)
            length = self.mesh.get_edgelength(vol_id, edge_id)
            normal_flow = num.dot(momentum, normal)*length

            # Reverse sign so that + is taken to mean inflow
            # and - means outflow. This is more intuitive.
            edge_flow = -normal_flow

            # Tally up inflows and outflows separately
            if edge_flow > 0:
                # Flow is inflow
                total_boundary_inflow += edge_flow
            else:
                # Flow is outflow
                total_boundary_outflow += edge_flow

            # Tally up flows by boundary tag
            tag = self.boundary[(vol_id, edge_id)]

            if tag not in boundary_flows:
                boundary_flows[tag] = 0.0
            boundary_flows[tag] += edge_flow


        return boundary_flows, total_boundary_inflow, total_boundary_outflow


    def compute_total_volume(self, region=None,
                             indices: list[int] | num.ndarray | None = None) -> float:
        """Total volume (m^3) of water in the domain, or in a region of it.

        An alias of :meth:`get_water_volume`; see there for the arguments.
        """

        return self.get_water_volume(region=region, indices=indices)


    def volumetric_balance_statistics(self) -> str:
        """Create volumetric balance report suitable for printing or logging.
        """

        (boundary_flows, total_boundary_inflow,
         total_boundary_outflow) = self.compute_boundary_flows()

        message = '---------------------------\n'
        message += 'Volumetric balance report:\n'
        message += 'Note: Boundary fluxes are not exact\n'
        message += 'See get_boundary_flux_integral for exact computation\n'
        message += '--------------------------\n'
        message += 'Total boundary inflow [m^3/s]: %.2f\n' % total_boundary_inflow
        message += 'Total boundary outflow [m^3/s]: %.2f\n' % total_boundary_outflow
        message += 'Net boundary flow by tags [m^3/s]\n'
        for tag in boundary_flows:
            message += '    %s [m^3/s]: %.2f\n' % (tag, boundary_flows[tag])

        message += 'Total net boundary flow [m^3/s]: %.2f\n' % \
                    (total_boundary_inflow + total_boundary_outflow)
        message += 'Total volume in domain [m^3]: %.2f\n' % \
                    self.compute_total_volume()

        # The go through explicit forcing update and record the rate of change
        # for stage and
        # record into forcing_inflow and forcing_outflow. Finally compute
        # integral of depth to obtain total volume of domain.

        # FIXME(Ole): This part is not yet done.

        return message

    def print_volumetric_balance_statistics(self) -> None:

        print (self.volumetric_balance_statistics())

    def report_water_volume_statistics(self, verbose: bool = True,
                                       returnStats: bool = False) -> list[float] | None:
        """
        Compute the volume, boundary flux integral, fractional step volume integral, and their difference

        If verbose, print a summary
        If returnStats, return a list with the volume statistics
        """
        from anuga import myid

        if(self.compute_fluxes_method != 'DE'):
            if(myid == 0):
                print('Water_volume_statistics only supported for DE algorithm ')
            return

        # Compute the volume
        Vol = self.get_water_volume()

        # Compute the boundary flux integral
        fluxIntegral=self.get_boundary_flux_integral()
        fracIntegral=self.get_fractional_step_volume_integral()

        if(verbose and myid==0):
            print(' ')
            print(f'    Volume V at time {self.get_time()}:', Vol)
            print('    Boundary Flux integral BF: ', fluxIntegral)
            print('    (rate + inlet) Fractional Step volume integral FS: ', fracIntegral)
            print('    V - BF - FS - InitialVolume :',  Vol- fluxIntegral -fracIntegral - self.volume_history[0])
            print(' ')

        if returnStats:
            return [Vol, fluxIntegral, fracIntegral]
        else:
            return

    def report_cells_with_small_local_timestep(self, threshold_depth: float | None = None) -> None:
        """
        Convenience function to print the locations of cells
        with a small local timestep.

        Computations are at cell centroids

        Useful in models
        with complex meshes, to find ways to speed up the model
        """

        from anuga.parallel import myid, numprocs
        from anuga.config import g, epsilon

        if threshold_depth is None:
            threshold_depth=self.minimum_allowed_height

        uh = self.quantities['xmomentum'].centroid_values
        vh = self.quantities['ymomentum'].centroid_values
        d =  self.quantities['stage'].centroid_values - self.quantities['elevation'].centroid_values
        d = num.maximum(d, threshold_depth)
        v = ((uh)**2 + (vh)**2)**0.5/d
        v = v * (d > threshold_depth)

        for i in range(numprocs):
            if myid == i:
                print('    Processor ', myid)
                gravSpeed = (g * d)**0.5
                waveSpeed = abs(v) + gravSpeed
                localTS = self.radii / num.maximum(waveSpeed, epsilon)
                controlling_pt_ind = localTS.argmin()
                print(f'     * Smallest LocalTS at time {self.get_time()} is approximately: ', localTS[controlling_pt_ind])
                print('     -- Location: ', round(self.centroid_coordinates[controlling_pt_ind,0]+self.geo_reference.xllcorner,2),\
                                        round(self.centroid_coordinates[controlling_pt_ind,1]+self.geo_reference.yllcorner,2))
                print('     -+ Speed: ', v[controlling_pt_ind])
                print('     -* Gravity_wave_speed', gravSpeed[controlling_pt_ind])
                print(' ')

            barrier()
        return

    def diagnose_timestep(self,
                          threshold_dt: float | None = None,
                          top_n: int = 5,
                          threshold_depth: float | None = None,
                          verbose: bool = True) -> dict:
        """Diagnose CFL-limited timestep and NaN/inf in conserved quantities.

        Computes the CFL-limiting dt per centroid (radius / wave_speed) and
        reports the cells driving the global minimum. Also scans stage,
        xmomentum and ymomentum for NaN / inf, which is the usual signature
        of an ADER-2 (or any) timestepping blow-up. Intended to be called
        from inside an evolve loop (e.g. each yieldstep) to localise where
        and when the solver is about to fail.

        Parameters
        ----------
        threshold_dt : float, optional
            If given, a WARNING is emitted whenever the global min local-dt
            falls below this value. Useful to catch the moment dt collapses.
        top_n : int
            Number of worst (smallest-dt) cells per rank to include in the
            returned dict and printed report.
        threshold_depth : float, optional
            Depths below this are clamped when computing speed, to avoid
            division-by-tiny. Defaults to self.minimum_allowed_height.
        verbose : bool
            If True, print a per-rank report.

        Returns
        -------
        dict
            Keys::

              time              relative simulation time
              current_timestep  self.timestep (last accepted global dt)
              cfl               self.CFL setting
              local_min_dt      smallest cell-local CFL dt on this rank
              global_min_dt     smallest cell-local CFL dt across all ranks
              local_max_speed   max wave_speed on this rank
              global_max_speed  max wave_speed across all ranks
              nan_counts        {'stage': n, 'xmomentum': n, 'ymomentum': n}
              first_nan_cell    centroid xy of the first NaN/inf cell on
                                this rank, or None
              worst_cells       list of dicts (top_n) with keys
                                {triangle, x, y, depth, speed, local_dt}
        """
        from anuga import numprocs
        from anuga.config import g, epsilon

        if threshold_depth is None:
            threshold_depth = self.minimum_allowed_height

        stage_c = self.quantities['stage'].centroid_values
        elev_c = self.quantities['elevation'].centroid_values
        uh_c = self.quantities['xmomentum'].centroid_values
        vh_c = self.quantities['ymomentum'].centroid_values

        bad_stage = ~num.isfinite(stage_c)
        bad_uh = ~num.isfinite(uh_c)
        bad_vh = ~num.isfinite(vh_c)
        nan_counts = {
            'stage':     int(bad_stage.sum()),
            'xmomentum': int(bad_uh.sum()),
            'ymomentum': int(bad_vh.sum()),
        }

        first_nan_cell = None
        any_bad = bad_stage | bad_uh | bad_vh
        if any_bad.any():
            i = int(num.argmax(any_bad))
            x = self.centroid_coordinates[i, 0] + self.geo_reference.xllcorner
            y = self.centroid_coordinates[i, 1] + self.geo_reference.yllcorner
            first_nan_cell = (float(x), float(y), i)

        # Local CFL dt per cell. Mask NaNs out so argmin/min are meaningful.
        d = stage_c - elev_c
        d_safe = num.where(num.isfinite(d), num.maximum(d, threshold_depth),
                           threshold_depth)
        speed_h = num.sqrt(uh_c * uh_c + vh_c * vh_c) / d_safe
        speed_h = num.where(num.isfinite(speed_h), speed_h, 0.0)
        speed_h = speed_h * (d > threshold_depth)
        grav_speed = num.sqrt(g * d_safe)
        wave_speed = num.abs(speed_h) + grav_speed

        local_dt = self.radii / num.maximum(wave_speed, epsilon)
        # Don't let NaN cells dominate the argmin
        local_dt = num.where(num.isfinite(local_dt), local_dt, num.inf)

        local_min_dt = float(local_dt.min())
        local_max_speed = float(num.nanmax(wave_speed)) if wave_speed.size else 0.0

        order = num.argsort(local_dt)[:max(top_n, 1)]
        worst_cells = []
        for idx in order:
            i = int(idx)
            worst_cells.append({
                'triangle': i,
                'x': float(self.centroid_coordinates[i, 0]
                           + self.geo_reference.xllcorner),
                'y': float(self.centroid_coordinates[i, 1]
                           + self.geo_reference.yllcorner),
                'depth': float(d[i]),
                'speed': float(wave_speed[i]),
                'local_dt': float(local_dt[i]),
            })

        if numprocs > 1:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            global_min_dt = comm.allreduce(local_min_dt, op=MPI.MIN)
            global_max_speed = comm.allreduce(local_max_speed, op=MPI.MAX)
            global_nan = {
                k: comm.allreduce(v, op=MPI.SUM) for k, v in nan_counts.items()
            }
        else:
            global_min_dt = local_min_dt
            global_max_speed = local_max_speed
            global_nan = dict(nan_counts)

        result = {
            'time':             self.get_time(),
            'current_timestep': self.timestep,
            'cfl':              self.CFL,
            'local_min_dt':     local_min_dt,
            'global_min_dt':    global_min_dt,
            'local_max_speed':  local_max_speed,
            'global_max_speed': global_max_speed,
            'nan_counts':       nan_counts,
            'global_nan_counts': global_nan,
            'first_nan_cell':   first_nan_cell,
            'worst_cells':      worst_cells,
        }

        if verbose:
            from anuga.parallel import myid
            for i in range(numprocs):
                if myid == i:
                    print(f'  [diagnose_timestep] rank {myid} t={result["time"]:.3f}')
                    print(f'    domain.timestep={self.timestep:.6e}  CFL={self.CFL}')
                    print(f'    local  min cell-dt = {local_min_dt:.6e}   max wave_speed = {local_max_speed:.4f}')
                    if numprocs > 1:
                        print(f'    global min cell-dt = {global_min_dt:.6e}   max wave_speed = {global_max_speed:.4f}')
                    total_bad = sum(nan_counts.values())
                    if total_bad:
                        print(f'    NaN/inf cells: stage={nan_counts["stage"]} '
                              f'xmom={nan_counts["xmomentum"]} '
                              f'ymom={nan_counts["ymomentum"]}')
                        if first_nan_cell is not None:
                            fx, fy, fi = first_nan_cell
                            print(f'    first bad cell: tri={fi} at ({fx:.2f}, {fy:.2f})')
                    for c in worst_cells:
                        print(f'    tri {c["triangle"]:>7}  '
                              f'xy=({c["x"]:.2f}, {c["y"]:.2f})  '
                              f'depth={c["depth"]:.4f}  speed={c["speed"]:.4f}  '
                              f'local_dt={c["local_dt"]:.6e}')
                    if threshold_dt is not None and global_min_dt < threshold_dt:
                        print(f'    WARNING: global min cell-dt {global_min_dt:.3e} '
                              f'< threshold {threshold_dt:.3e}')
                    sys.stdout.flush()
                barrier()

        return result

# =======================================================================
# PETE: NEW METHODS FOR FOR PARALLEL STRUCTURES. Note that we assume the
# first "number_of_full_[nodes|triangles]" are full [nodes|triangles]
# For full triangles it is possible to enquire self.tri_full_flag == True
# =======================================================================

    def get_number_of_full_triangles(self, *args, **kwargs) -> int:
        return self.number_of_full_triangles

    def get_full_centroid_coordinates(self, *args, **kwargs) -> num.ndarray:
        C = self.mesh.get_centroid_coordinates(*args, **kwargs)
        return C[:self.number_of_full_triangles, :]

    def get_full_vertex_coordinates(self, *args, **kwargs) -> num.ndarray:
        V = self.mesh.get_vertex_coordinates(*args, **kwargs)
        return V[:3*self.number_of_full_triangles,:]

    def get_full_triangles(self, *args, **kwargs) -> num.ndarray:
        T = self.mesh.get_triangles(*args, **kwargs)
        return T[:self.number_of_full_triangles,:]

    def get_full_nodes(self, *args, **kwargs) -> num.ndarray:
        N = self.mesh.get_nodes(*args, **kwargs)
        return N[:self.number_of_full_nodes,:]

    def get_tri_map(self) -> num.ndarray | None:
        return self.tri_map

    def get_inv_tri_map(self) -> num.ndarray | None:
        return self.inv_tri_map

# ==============================================================================
# Multiprocessor Mode (1=openmp, 2=cupy (in development))
# ==============================================================================

    # User-facing *per-domain* compute mode. This is the only genuinely
    # per-domain choice — it selects the internal ``multiprocessor_mode``:
    #   'legacy'  -> mode 1: sw_domain_openmp_ext solver + serial-Python operators
    #   'unified' -> mode 2: the unified gpu_ext C kernels (solver + operators)
    # Whether 'unified' runs on CPU or offloads to a GPU is NOT a per-domain
    # property: it is a *process-global* OpenMP offload setting controlled by
    # the module function ``set_gpu_offload()`` (and only possible on a GPU
    # build). On a CPU-only build, 'unified' simply runs CPU-multicore.
    COMPUTE_MODES = ('legacy', 'unified')

    def _mode2_mpi_available(self) -> bool:
        """True if ``sw_domain_gpu_ext`` was built with real C MPI support.

        Mode 2 ('unified') performs its halo exchange at the C level
        (``exchange_ghosts``). Without an MPI-enabled build that exchange is a
        silent no-op, so a multi-rank 'unified' run would compute wrong results —
        the selector falls back to 'legacy' (Python MPI exchange) in that case.
        Serial (single-rank) 'unified' needs no MPI and is unaffected.
        """
        try:
            from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext
            return bool(gpu_ext.gpu_has_mpi())
        except Exception:
            return False

    def compute_capabilities(self) -> dict:
        """Report which compute backends this build/run supports.

        Returns a dict with::

            'gpu_offload'     : bool - process can offload mode-2 to a GPU device
                                       (build supports it, device present, offload
                                       not disabled); see set_gpu_offload()
            'num_gpu_devices' : int  - number of offload devices visible
            'mpi'             : bool - gpu_ext built with C MPI ('unified' parallel ok)
            'modes'           : list - per-domain modes available ('unified' only
                                       when the gpu_ext extension is importable)
        """
        try:
            from anuga.shallow_water import sw_domain_gpu_ext as gpu_ext  # noqa: F401
            unified = True
            ndev = int(gpu_ext.get_num_gpu_devices())
            mpi = bool(gpu_ext.gpu_has_mpi())
        except Exception:
            unified, ndev, mpi = False, 0, False
        modes = ['legacy'] + (['unified'] if unified else [])
        return {'gpu_offload': gpu_offload_enabled(), 'num_gpu_devices': ndev,
                'mpi': mpi, 'modes': modes}

    def set_compute_mode(self, mode: str = 'unified', verbose: bool = False,
                         strict: bool | None = None) -> None:
        """Select this domain's compute mode (per-domain).

        This is a per-domain setting — different domains in one script may use
        different modes. Whether 'unified' uses a GPU is a separate, process-wide
        decision (see :func:`set_gpu_offload`), because OpenMP target offload is
        a process-level runtime setting, not a per-domain one.

        Under MPI, 'unified' requires a gpu_ext built with MPI; otherwise this
        falls back to 'legacy' (whose Python MPI exchange is correct in parallel)
        with a rank-0 warning. The active mode is recorded in
        ``self.compute_mode``; the original request in
        ``self.requested_compute_mode``.

        Parameters
        ----------
        mode : {'legacy', 'unified'}
            - ``'legacy'`` — mode 1: the ``sw_domain_openmp_ext`` solver with
              serial-Python fractional-step operators.
            - ``'unified'`` — mode 2: the unified ``sw_domain_gpu_ext`` C kernels
              (solver and operators). Runs CPU-multicore by default; offloads to
              a GPU only when GPU offload is enabled process-wide via
              :func:`anuga.set_gpu_offload` on a GPU-capable build.
        strict : bool, optional
            Fail instead of falling back. With ``strict=True`` every silent
            downgrade of the requested mode becomes a ``RuntimeError``: the
            MPI fallback from 'unified' to 'legacy', running on the CPU when
            GPU offload is enabled but no device was found, and boundary types
            the device path cannot evaluate (which otherwise run on the host
            every step). Defaults to the domain's current setting, initially
            ``ANUGA_STRICT_COMPUTE_MODE`` from the environment. Meant for
            batch jobs, where a run that quietly does something slower than
            what was asked for is the expensive failure.
        """
        import warnings

        if mode not in self.COMPUTE_MODES:
            raise ValueError(
                f"Invalid compute mode {mode!r}. Must be one of {self.COMPUTE_MODES}.")

        if strict is None:
            strict = getattr(self, 'compute_mode_strict', _strict_compute_mode_from_env())
        self.compute_mode_strict = bool(strict)

        requested = mode
        if mode == 'unified':
            # Parallel guard: 'unified' exchanges ghosts at the C level, a silent
            # no-op without an MPI-enabled gpu_ext build. Under MPI that gives
            # wrong results, so fall back to 'legacy' (Python MPI exchange).
            try:
                from anuga import numprocs
            except Exception:
                numprocs = 1
            if numprocs > 1 and not self._mode2_mpi_available():
                try:
                    from anuga import myid
                except Exception:
                    myid = 0
                msg = (f"compute mode 'unified' selected under MPI ({numprocs} ranks) but "
                       "this ANUGA build's sw_domain_gpu_ext was compiled without MPI; the "
                       "C-level ghost exchange would be a silent no-op and give wrong "
                       "parallel results. Falling back to 'legacy' (mode 1, Python MPI "
                       "exchange). Rebuild with MPI to run 'unified' in parallel.")
                if self.compute_mode_strict:
                    raise RuntimeError(msg + " (strict mode: refusing to fall back)")
                if myid == 0:
                    warnings.warn(msg, stacklevel=2)
                mode = 'legacy'

        self.requested_compute_mode = requested
        self.compute_mode = mode

        if mode == 'legacy':
            self.multiprocessor_mode = MULTIPROCESSOR_OPENMP
            self.use_c_rk_loop = False
        else:  # 'unified' -> mode 2 (unified gpu_ext kernels)
            self.multiprocessor_mode = MULTIPROCESSOR_GPU
            self.use_c_rk_loop = True
            # Build the device interface now if boundaries are ready; otherwise
            # defer to the first evolve(). Boundaries are typically set AFTER
            # construction, so a default-'unified' domain must not require them
            # at __init__ time.
            if self._boundaries_ready():
                self.set_gpu_interface()

        if verbose:
            print(f"Compute mode: requested {requested!r} -> active {self.compute_mode!r} "
                  f"(multiprocessor_mode={self.multiprocessor_mode}, "
                  f"gpu_offload={'on' if gpu_offload_enabled() else 'off'})")

    def _boundaries_ready(self) -> bool:
        """True if real boundary objects are set — required to build the gpu_ext
        device interface. After distribute() boundary_map may be
        ``{'exterior': None, 'ghost': None}`` (not None but no real boundaries).
        """
        bmap = getattr(self, 'boundary_map', None)
        return bool(bmap) and any(b is not None for b in bmap.values())

    def _ensure_gpu_interface(self) -> None:
        """Build a deferred mode-2 device interface on demand.

        For a default-'unified' domain the interface is built lazily (boundaries
        are set after construction). Build it the first time a mode-2 path needs
        it. If boundaries are still not set — e.g. a test that pokes the domain
        without a full boundary setup — fall back to 'legacy' so the operation
        can proceed rather than hit a None gpu_interface.
        """
        if self.multiprocessor_mode != MULTIPROCESSOR_GPU or self.gpu_interface is not None:
            return
        if self._boundaries_ready():
            self.set_gpu_interface()
        else:
            self.set_compute_mode('legacy')

    def get_compute_mode(self) -> str:
        """Return the active per-domain compute mode: 'legacy' or 'unified'."""
        return getattr(self, 'compute_mode', 'legacy')

    def set_multiprocessor_mode(self, multiprocessor_mode: int = 1) -> None:
        """
        Set multiprocessor mode (legacy integer API).

        1. openmp - Python RK loop (use_c_rk_loop=False)
        2. gpu/mpi - C RK loop (use_c_rk_loop=True)

        Thin wrapper over :meth:`set_compute_mode`: 1 maps to ``'legacy'``, 2 to
        ``'unified'``. Whether 'unified' offloads to a GPU is a separate,
        process-wide choice — see :func:`anuga.set_gpu_offload`. New code should
        prefer :meth:`set_compute_mode`.
        """
        if multiprocessor_mode not in [MULTIPROCESSOR_OPENMP, MULTIPROCESSOR_GPU]:
            raise ValueError('Invalid multiprocessor mode. Must be one of [1,2] (openmp, gpu/mpi)')

        if multiprocessor_mode == MULTIPROCESSOR_OPENMP:
            self.set_compute_mode('legacy')
        else:
            self.set_compute_mode('unified')

    @property
    def use_c_rk2_loop(self):
        """Deprecated: use use_c_rk_loop instead."""
        import warnings
        warnings.warn(
            "use_c_rk2_loop is deprecated; use use_c_rk_loop instead.",
            DeprecationWarning, stacklevel=2)
        return self.use_c_rk_loop

    @use_c_rk2_loop.setter
    def use_c_rk2_loop(self, value):
        import warnings
        warnings.warn(
            "use_c_rk2_loop is deprecated; use use_c_rk_loop instead.",
            DeprecationWarning, stacklevel=2)
        self.use_c_rk_loop = value

    def get_multiprocessor_mode(self) -> int:
        """
        Get multiprocessor mode

        1. openmp (in development)
        2. gpu/mpi (in development)
        """
        return self.multiprocessor_mode

    @property
    def omp_num_threads(self) -> int:
        """The process-wide OpenMP thread count (read-only view).

        OpenMP thread count is a process-level setting, not per-domain, so this
        reflects the live value set by :func:`anuga.set_omp_num_threads` for
        *every* domain in the session — including domains constructed before the
        call. Assigning to it (``domain.omp_num_threads = n``) is kept for
        backward compatibility and sets the count process-wide.
        """
        return get_omp_num_threads()

    @omp_num_threads.setter
    def omp_num_threads(self, value: int) -> None:
        set_omp_num_threads(value, verbose=False)

    def set_omp_num_threads(self, omp_num_threads: int | None = None, verbose: bool = True) -> None:
        """Set the OpenMP thread count (process-wide).

        OpenMP thread count is a process-level setting, not per-domain. This is a
        thin wrapper that delegates to the module-level
        :func:`anuga.set_omp_num_threads`; prefer that in new code. Kept for
        backward compatibility.
        """
        set_omp_num_threads(omp_num_threads, verbose=verbose)


    @property
    def gpu(self):
        """
        Shortcut to access the GPU interface.

        Returns the gpu_interface object which provides FLOP counters and kernel wrappers.
        Raises AttributeError if GPU mode is not enabled.
        """
        if self.gpu_interface is None:
            raise AttributeError(
                "GPU interface not initialized. Call domain.set_multiprocessor_mode(2) first."
            )
        return self.gpu_interface

    def set_gpu_interface(self):

        if self.multiprocessor_mode == MULTIPROCESSOR_GPU and self.gpu_interface is None:

            # Check that boundaries are properly set before GPU initialization
            # After distribute(), boundary_map may be {'exterior': None, 'ghost': None}
            # which is not None but has no actual boundary objects - this causes silent failures
            if self.boundary_map is None:
                raise RuntimeError(
                    "GPU mode requires boundaries to be set before calling set_multiprocessor_mode(2).\n"
                    "Please call domain.set_boundary({...}) BEFORE domain.set_multiprocessor_mode(2)."
                )

            has_real_boundary = any(b is not None for b in self.boundary_map.values())
            if not has_real_boundary:
                raise RuntimeError(
                    "GPU mode requires boundaries to be set before calling set_multiprocessor_mode(2).\n"
                    "Please call domain.set_boundary({...}) BEFORE domain.set_multiprocessor_mode(2).\n"
                    f"Current boundary_map has no boundary objects: {list(self.boundary_map.keys())}"
                )

            # Reconcile deeply-dry cells (stage < bed) to stage = bed before the initial
            # state is synced to the device. A dry cell should carry stage = bed
            # (depth 0); mode 1 reaches that via its per-step protect on the very first
            # step, but the mode-2 device path only converges to it gradually (halving
            # the deficit each step, ~a dozen steps to close a large stage<<bed gap).
            # While it converges, any forcing applied to those cells — an Inlet_operator,
            # rainfall — is absorbed into raising the sub-bed stage rather than making
            # water depth, and is lost, leaving a permanent startup mass deficit
            # (issue #200). Clamping stage up to bed here is mass-neutral (depth stays 0)
            # and a no-op for wet cells, and makes mode 2 start already reconciled.
            stage_c = self.quantities['stage'].centroid_values
            bed_c = self.quantities['elevation'].centroid_values
            if (stage_c < bed_c).any():
                self.set_quantity('stage', num.maximum(stage_c, bed_c),
                                  location='centroids')

            # Try OpenMP target offloading interface first. Only the build and
            # setup are guarded: the checks that follow must not be swallowed
            # by this except, which exists for "the unified extension is not
            # available here", not for a strict-mode refusal.
            omp_ok = False
            try:
                from .sw_domain_gpu_omp import GPU_OMP_interface
                self.gpu_interface = GPU_OMP_interface(self)
                self.gpu_interface.setup()
                omp_ok = True
            except Exception as e:
                self.gpu_interface = None
                print(f'OpenMP GPU interface not available: {e}')

            if omp_ok:
                from anuga import myid, numprocs
                omp_num_threads = os.environ.get('OMP_NUM_THREADS', '1')
                # Offload is a process-wide decision (set_gpu_offload), resolved
                # against the build. False on a CPU-only build: 'unified' is CPU
                # multicore, not GPU.
                self.gpu_offload_active = gpu_offload_enabled()
                device_id = self.gpu_interface.gpu_dom.device_id
                if self.compute_mode_strict and self.gpu_offload_active and device_id < 0:
                    raise RuntimeError(
                        "GPU offload is enabled but no GPU device was found, so this "
                        "'unified' domain would run on the CPU via the OpenMP target host "
                        "fallback (strict mode: refusing to fall back). Check nvidia-smi, "
                        "the CUDA driver, and that this is a GPU-offload build.")
                # Only print from rank 0
                if myid == 0:
                    # The number of GPUs the runtime can actually see — NOT numprocs.
                    try:
                        from anuga.shallow_water.sw_domain_gpu_ext import get_num_gpu_devices
                        num_devices = get_num_gpu_devices()
                    except Exception:
                        num_devices = -1

                    for line in gpu_startup_banner(numprocs, num_devices, device_id,
                                                   self.gpu_offload_active, omp_num_threads):
                        print(line)

                # Classify the boundaries now rather than on the first evolve
                # step, so an unsupported type is reported (or, in strict mode,
                # refused) at setup, before any batch time is spent.
                self._init_gpu_boundary_info(self.gpu_interface.gpu_dom)
                return

            # Fall back to CUDA/CuPy interface
            try:
                import cupy as cp
                test_cupy_array = cp.array([1,2,3])
                from .sw_domain_cuda import GPU_interface
                self.gpu_interface = GPU_interface(self)
                self.gpu_interface.allocate_gpu_arrays()
                self.gpu_interface.compile_gpu_kernels()
                return
            except Exception:
                pass

            # No unified interface available at all
            if self.compute_mode_strict:
                raise RuntimeError(
                    "compute mode 'unified' was requested but no unified (mode-2) "
                    "interface could be set up in this build, so the domain would "
                    "fall back to 'legacy' (strict mode: refusing to fall back).")
            from anuga import myid
            if myid == 0:
                print('+==============================================================================+')
                print('| WARNING: GPU not available, falling back to multiprocessor_mode 1 (OpenMP)  |')
                print('+==============================================================================+')
            self.set_multiprocessor_mode(1)



################################################################################
# End of class Shallow Water Domain
################################################################################



# FIXME (Ole): Does this do anything?

def my_update_special_conditions(domain):

    pass



if __name__ == "__main__":
    pass
