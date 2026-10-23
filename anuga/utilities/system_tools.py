"""Implementation of tools to do with system administration made as platform independent as possible.
"""

import sys
import os
import platform
import warnings

# Record Python version
major_version = int(platform.python_version_tuple()[0])
version = platform.python_version()

import hashlib


def log_to_file(filename, s, verbose=False, mode='a'):
    """Log string to file name
    """

    with open(filename, mode) as fid:
        if verbose: print(s)
        fid.write(s + '\n')


def get_user_name():
    """Get user name provide by operating system
    """

    import getpass
    user = getpass.getuser()

    return user

def get_host_name():
    """Get host name provide by operating system
    """

    if sys.platform == 'win32':
        host = os.getenv('COMPUTERNAME')
    else:
        host = os.uname()[1]


    return host


def get_version():
    """Get anuga version number as stored in anuga.__version__
    """

    import anuga
    return anuga.__version__


def get_revision_number():
    """Get the (git) sha of this repository copy.
    """
    from anuga import __git_sha__ as revision
    return revision


def get_revision_date():
    """Get the (git) revision date of this repository copy.
    """

    from anuga import __git_committed_datetime__ as revision_date
    return revision_date



def safe_crc(string):
    """64 bit safe crc computation.

    See http://docs.python.org/library/zlib.html#zlib.crc32:

        To generate the same numeric value across all Python versions
        and platforms use crc32(data) & 0xffffffff.
    """

    from zlib import crc32

    return crc32(string) & 0xffffffff


def compute_checksum(filename, max_length=2**20):
    """Compute the CRC32 checksum for specified file

    Optional parameter max_length sets the maximum number
    of bytes used to limit time used with large files.
    Default = 2**20 (1MB)
    """

    fid = open(filename, 'rb') # Use binary for portability
    crcval = safe_crc(fid.read(max_length))
    fid.close()

    return crcval


def get_anuga_pathname():
    """Get pathname of anuga install location

    Typically, this is required in unit tests depending
    on external files.

    """

    import anuga

    return os.path.dirname(anuga.__file__)


def get_pathname_from_package(package):
    """Get pathname of given package (provided as string)

    This is useful for reading files residing in the same directory as
    a particular module. Typically, this is required in unit tests depending
    on external files.

    The given module must start from a directory on the pythonpath
    and be importable using the import statement.

    Example
    path = get_pathname_from_package('anuga.utilities')

    """

    import importlib

    module = importlib.import_module(package)
    return os.path.dirname(module.__file__)


def clean_line(str, delimiter):
    """Split a string into 'clean' fields.

    str        the string to process
    delimiter  the delimiter string to split 'line' with

    Returns a list of 'cleaned' field strings.

    Any fields that were initially zero length will be removed.
    If a field contains '\n' it isn't zero length.
    """

    return [x.strip() for x in str.strip().split(delimiter) if x != '']


################################################################################
# The following two functions are used to get around a problem with numpy and
# NetCDF files.  Previously, using Numeric, we could take a list of strings and
# convert to a Numeric array resulting in this:
#     Numeric.array(['abc', 'xy']) -> [['a', 'b', 'c'],
#                                      ['x', 'y', ' ']]
#
# However, under numpy we get:
#     numpy.array(['abc', 'xy']) -> ['abc',
#                                    'xy']
#
# And writing *strings* to a NetCDF file is problematic.
#
# The solution is to use these two routines to convert a 1-D list of strings
# to the 2-D list of chars form and back.  The 2-D form can be written to a
# NetCDF file as before.
#
# The other option, of inverting a list of tag strings into a dictionary with
# keys being the unique tag strings and the key value a list of indices of where
# the tag string was in the original list was rejected because:
#    1. It's a lot of work
#    2. We'd have to rewite the I/O code a bit (extra variables instead of one)
#    3. The code below is fast enough in an I/O scenario
################################################################################

def string_to_char(l):
    """Convert 1-D list of strings to 2-D list of chars."""

    if not l:
        return []

    if l == ['']:
        l = [' ']


    maxlen = max(len(x) for x in l)
    ll = [x.ljust(maxlen) for x in l]
    result = []
    for s in ll:
        result.append([x for x in s])
    return result



def char_to_string(ll):
    """Convert 2-D list of chars to 1-D list of strings."""

    # https://stackoverflow.com/questions/23618218/numpy-bytes-to-plain-string
    # bytes_string.decode('UTF-8')

    # We might be able to do this a bit more shorthand as we did in Python2.x
    # i.e return [''.join(x).strip() for x in ll]

    # But this works for now.

    result = []
    for i in range(len(ll)):
        x = ll[i]
        string = ''
        for j in range(len(x)):
            c = x[j]
            if type(c) == str:
                string += c
            else:
                string += c.decode()

        result.append(string.strip())

    return result


################################################################################

def get_vars_in_expression(source):
    """Get list of variable names in a python expression."""

    # https://stackoverflow.com/questions/37993137/how-do-i-detect-variables-in-a-python-eval-expression

    import ast

    variables = {}
    syntax_tree = ast.parse(source)
    for node in ast.walk(syntax_tree):
        if type(node) is ast.Name:
            variables[node.id] = 0  # Keep first one, but not duplicates

    # Only return keys
    result = list(variables.keys()) # Only return keys i.e. the variable names
    result.sort() # Sort for uniqueness
    return result



def file_length(in_file):
    """Function to return the length of a file."""

    fid = open(in_file)
    data = fid.readlines()
    fid.close()
    return len(data)


#### Memory functions
_proc_status = '/proc/%d/status' % os.getpid()
_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    """private method"""
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except OSError:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to MB
    return (float(v[1]) * _scale[v[2]]) // (1024.0 * 1024.0)



def _get_rss_mb():
    """Return current resident set size in MB.

    Uses psutil if available, otherwise falls back to /proc/self/status on
    Linux.  Returns 0.0 if neither source is available.
    """
    try:
        import psutil
        return psutil.Process(os.getpid()).memory_info().rss / 1024**2
    except ImportError:
        rss = _VmB('VmRSS:')
        return rss  # already in MB from _VmB


def memory_stats():
    """Return a formatted string showing the current process memory usage.

    Uses the resident set size (RSS) — physical RAM currently occupied by
    the process.  Values below 1 GB are reported in MB; values at or above
    1 GB are reported in GB.

    Returns
    -------
    str
        Memory usage string, e.g. ``'mem=342MB'`` or ``'mem=1.50GB'``.

    Examples
    --------
    >>> from anuga import memory_stats
    >>> print(memory_stats())
    mem=342MB
    """
    rss_mb = _get_rss_mb()
    if rss_mb >= 1024:
        return f'mem={rss_mb / 1024:.2f}GB'
    return f'mem={rss_mb:.0f}MB'


def print_memory_stats():
    """Print the current process memory usage to stdout.

    Convenience wrapper around :func:`memory_stats`.

    Examples
    --------
    >>> from anuga import print_memory_stats
    >>> print_memory_stats()
    mem=342MB
    """
    print(memory_stats())


def quantity_memory_stats(domain):
    """Return a detailed breakdown of memory used by all quantities on *domain*.

    For each registered quantity the function reports which backing arrays are
    allocated and their individual sizes, then prints a per-quantity subtotal
    and a grand total.  Arrays that have been kept lazy (``None``) are shown
    as ``--`` with zero cost.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain whose quantities are to be inspected.

    Returns
    -------
    str
        A multi-line table suitable for printing or logging.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> print(anuga.quantity_memory_stats(domain))
    Quantity memory breakdown  (N=400 triangles, L=80 boundary edges)
    ...
    """
    # Abbreviated column headers (attr name on the Quantity object).
    _ARRAYS = [
        ('cntrd',  'centroid_values'),
        ('edge',   'edge_values'),
        ('vert',   '_vertex_values'),
        ('bndry',  'boundary_values'),
        ('ex_up',  'explicit_update'),
        ('si_up',  'semi_implicit_update'),
        ('ct_bk',  'centroid_backup_values'),
        ('x_grd',  '_x_gradient'),
        ('y_grd',  '_y_gradient'),
        ('phi',    '_phi'),
    ]

    # Abbreviated quantity type names.
    _TYPE_ABBREV = {
        'evolved':          'evol',
        'edge_diagnostic':  'e_diag',
        'centroid_only':    'c_only',
        'coordinate':       'coord',
    }

    N = domain.number_of_elements
    try:
        L = domain.boundary_length
    except AttributeError:
        L = 0

    col_w  = max(len(label) for label, _ in _ARRAYS) + 2   # 7
    name_w = max(len(n) for n in domain.quantities) + 1
    type_w = max(len(_TYPE_ABBREV.get(getattr(q, '_qty_type', '?'), '?'))
                 for q in domain.quantities.values()) + 1

    headers = [label for label, _ in _ARRAYS]
    header_row = (f"{'Quantity':<{name_w}} {'type':<{type_w}} " +
                  ''.join(f"{h:>{col_w}}" for h in headers) +
                  f"{'TOTAL':>{col_w}}")
    sep = '-' * len(header_row)

    lines = [
        f"Quantity memory (N={N:,} tri, L={L:,} bndry) — kB, '--'=lazy",
        "  type: evol=evolved  e_diag=edge_diagnostic  c_only=centroid_only  coord=coordinate",
        "  cols: cntrd=centroid  edge=edge  vert=vertex  bndry=boundary  ex_up=explicit_update",
        "        si_up=semi_implicit_update  ct_bk=centroid_backup  x/y_grd=gradients",
        sep,
        header_row,
        sep,
    ]

    grand_total_bytes = 0

    for qty_name, qty in sorted(domain.quantities.items()):
        raw_type = getattr(qty, '_qty_type', '?')
        qty_type = _TYPE_ABBREV.get(raw_type, raw_type)
        row_bytes = 0
        cells = []

        for _label, attr in _ARRAYS:
            arr = getattr(qty, attr, None)
            if arr is None:
                cells.append('--')
            else:
                b = arr.nbytes
                row_bytes += b
                cells.append(f'{b / 1024:.1f}')

        grand_total_bytes += row_bytes
        total_str = f'{row_bytes / 1024:.1f}'
        line = (f"{qty_name:<{name_w}} {qty_type:<{type_w}} " +
                ''.join(f"{c:>{col_w}}" for c in cells) +
                f"{total_str:>{col_w}}")
        lines.append(line)

    lines.append(sep)
    grand_kb = grand_total_bytes / 1024
    grand_mb = grand_total_bytes / 1024**2
    lines.append(f"{'GRAND TOTAL':<{name_w + type_w + 1}}"
                 f"{grand_kb:>{col_w * (len(_ARRAYS) + 1)},.1f} kB = {grand_mb:.2f} MB")
    lines.append(sep)

    return '\n'.join(lines)


def print_quantity_memory_stats(domain):
    """Print a detailed breakdown of memory used by all quantities on *domain*.

    Convenience wrapper around :func:`quantity_memory_stats`.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain whose quantities are to be inspected.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> anuga.print_quantity_memory_stats(domain)
    Quantity memory breakdown  (N=400 triangles, L=80 boundary edges)
    ...
    """
    print(quantity_memory_stats(domain))


def domain_memory_stats(domain):
    """Return a breakdown of memory used by the domain and all its quantities.

    Reports numpy array memory grouped by category (geometry, connectivity,
    work arrays, quantities, river-wall, flags/parallel), the process RSS, and
    the unaccounted gap (C structs, Python interpreter, imported libraries).

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain to inspect.

    Returns
    -------
    str
        A multi-line table suitable for printing or logging.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> print(anuga.domain_memory_stats(domain))
    Domain memory (N=400 tri) — MB
    ...
    """
    import numpy as np

    # ------------------------------------------------------------------ #
    # 1.  Categorise domain-level numpy arrays by name                    #
    # ------------------------------------------------------------------ #
    _CATEGORIES = {
        'geometry': [
            'vertex_coordinates', 'edge_coordinates', 'centroid_coordinates',
            'nodes', 'normals', 'areas', 'radii', 'edgelengths',
        ],
        'connectivity': [
            'triangles', 'neighbours', 'surrogate_neighbours',
            'neighbour_edges', 'vertex_value_indices', 'number_of_boundaries',
            'node_index', 'number_of_triangles_per_node',
        ],
        'work arrays': [
            'edge_flux_work', 'neigh_work', 'pressuregrad_work',
            'already_computed_flux', 'work_centroid_values',
            'x_centroid_work', 'y_centroid_work',
            '_grad_workspace_x', '_grad_workspace_y', '_phi_workspace',
            'max_speed',
        ],
        'river wall': [
            'edge_flux_type', 'edge_river_wall_counter',
        ],
        'flags/parallel': [
            'tri_full_flag', 'node_full_flag',
            'boundary_flux_sum',
            'boundary_cells', 'boundary_edges',
        ],
    }

    # Collect actual arrays from domain
    cat_bytes = {}
    cat_detail = {}
    accounted_attrs = set()

    for cat, attrs in _CATEGORIES.items():
        total = 0
        details = []
        for attr in attrs:
            arr = getattr(domain, attr, None)
            if isinstance(arr, np.ndarray):
                b = arr.nbytes
                total += b
                details.append((attr, arr.shape, arr.dtype, b))
                accounted_attrs.add(attr)
        cat_bytes[cat] = total
        cat_detail[cat] = details

    # Catch any numpy arrays not yet in our category list
    other = []
    other_bytes = 0
    for attr in sorted(dir(domain)):
        if attr.startswith('__') or attr in accounted_attrs:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', DeprecationWarning)
                val = getattr(domain, attr)
        except Exception:
            continue
        if isinstance(val, np.ndarray) and val.nbytes > 0:
            other_bytes += val.nbytes
            other.append((attr, val.shape, val.dtype, val.nbytes))

    # ------------------------------------------------------------------ #
    # 2.  Quantity totals (re-use quantity_memory_stats logic)            #
    # ------------------------------------------------------------------ #
    _QTY_ARRAYS = [
        'centroid_values', 'edge_values', '_vertex_values',
        'boundary_values', 'explicit_update', 'semi_implicit_update',
        'centroid_backup_values', '_x_gradient', '_y_gradient', '_phi',
    ]
    qty_bytes = 0
    for qty in domain.quantities.values():
        for attr in _QTY_ARRAYS:
            arr = getattr(qty, attr, None)
            if isinstance(arr, np.ndarray):
                qty_bytes += arr.nbytes

    # ------------------------------------------------------------------ #
    # 3.  Totals and RSS                                                  #
    # ------------------------------------------------------------------ #
    domain_numpy_bytes = sum(cat_bytes.values()) + other_bytes
    total_python_bytes = domain_numpy_bytes + qty_bytes
    rss_mb = _get_rss_mb()
    rss_bytes = rss_mb * 1024**2

    # ------------------------------------------------------------------ #
    # 4.  Format output                                                   #
    # ------------------------------------------------------------------ #
    N = domain.number_of_elements
    w = 10   # column width for MB values

    def _fmt(b):
        return f'{b / 1024**2:.2f}'

    sep = '-' * 55
    lines = [
        f'Domain memory (N={N:,} tri) — MB',
        sep,
        f'  {"Category":<22} {"MB":>{w}}',
        sep,
    ]

    for cat in _CATEGORIES:
        lines.append(f'  {cat:<22} {_fmt(cat_bytes[cat]):>{w}}')

    if other:
        lines.append(f'  {"other domain arrays":<22} {_fmt(other_bytes):>{w}}')

    gap_bytes = rss_bytes - total_python_bytes
    if gap_bytes >= 0:
        gap_label = 'unaccounted (C/libs)'
        gap_str   = _fmt(gap_bytes)
    else:
        # Allocated numpy bytes exceed RSS: the difference is virtual (cold)
        # memory — zero-initialised arrays whose pages have not yet been
        # paged into RAM (OS lazy commit; common before the first evolve call
        # for work arrays such as edge_flux_work, neigh_work, update arrays).
        gap_label = 'cold/virtual (untouch)'
        gap_str   = f'-{_fmt(-gap_bytes)}'

    lines += [
        f'  {"quantities":<22} {_fmt(qty_bytes):>{w}}',
        sep,
        f'  {"total Python numpy":<22} {_fmt(total_python_bytes):>{w}}',
        f'  {"process RSS":<22} {_fmt(rss_bytes):>{w}}',
        f'  {gap_label:<22} {gap_str:>{w}}',
        sep,
    ]

    if other:
        lines.append('  Uncategorised domain arrays:')
        for attr, shape, dtype, b in sorted(other, key=lambda x: -x[3]):
            lines.append(f'    {_fmt(b):>6} MB  shape={str(shape):<20} dtype={dtype}  {attr}')
        lines.append(sep)

    return '\n'.join(lines)


def print_domain_memory_stats(domain):
    """Print a breakdown of memory used by the domain and all its quantities.

    Convenience wrapper around :func:`domain_memory_stats`.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain to inspect.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> anuga.print_domain_memory_stats(domain)
    Domain memory (N=400 tri) — MB
    ...
    """
    print(domain_memory_stats(domain))


def domain_struct_stats(domain):
    """Return a breakdown of the C domain struct and GPU domain struct memory.

    The C domain struct (``_Domain_C_struct``) is a fixed-size block of
    scalars and raw C pointers that mirror the Python domain.  It adds only
    ~1 kB of overhead and holds **no additional array data** — all pointers
    point back into the numpy arrays already counted by
    :func:`domain_memory_stats`.

    The GPU domain struct (``gpu_dom``) extends the C struct with MPI/GPU
    state and boundary sub-structs.  For OpenMP GPU offloading the mapped
    numpy arrays live on the GPU device; ``estimate_required_memory`` gives
    the estimated device-side footprint.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain to inspect.

    Returns
    -------
    str
        A multi-line diagnostic string.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> print(anuga.domain_struct_stats(domain))
    C / GPU struct diagnostics (N=400 tri)
    ...
    """
    N = domain.number_of_elements
    try:
        nb = domain.boundary_length
    except AttributeError:
        nb = 0

    sep = '-' * 55
    lines = [f'C / GPU struct diagnostics (N={N:,} tri)', sep]

    # ------------------------------------------------------------------ #
    # 1.  C domain struct                                                  #
    # ------------------------------------------------------------------ #
    # sizeof(struct domain):
    #   10 int64 scalars + 13 double scalars + ~64 void* pointers
    #   10×8 + 13×8 + 64×8 = 696 bytes (+ alignment padding ≈ 720 B)
    _C_STRUCT_SCALARS_BYTES = (10 + 13) * 8
    _C_STRUCT_PTRS = 64     # pointer fields in struct domain
    _C_STRUCT_SIZE = _C_STRUCT_SCALARS_BYTES + _C_STRUCT_PTRS * 8  # ~696 B

    c_struct = getattr(domain, '_Domain_C_struct', None)
    c_status = 'allocated' if c_struct is not None else 'not allocated'
    lines += [
        '  C domain struct (struct domain)',
        f'    status          : {c_status}',
        f'    sizeof estimate : {_C_STRUCT_SIZE} B  '
        f'({10} int64 scalars, {13} double scalars, {_C_STRUCT_PTRS} pointers)',
        '    note            : pointers reference Python numpy arrays — no extra data',
    ]

    # ------------------------------------------------------------------ #
    # 2.  GPU domain struct                                                #
    # ------------------------------------------------------------------ #
    # sizeof(struct gpu_domain):
    #   struct domain (embedded) + MPI state + ~10 sub-structs + flop counters
    #   Rough estimate: 696 + ~800 ≈ 1.5 kB (struct overhead, not device data)
    _GPU_STRUCT_SIZE_EST = _C_STRUCT_SIZE + 800

    gpu_iface = getattr(domain, 'gpu_interface', None)
    lines.append('')
    lines.append('  GPU domain struct (struct gpu_domain)')

    if gpu_iface is None:
        lines.append('    status          : not initialised (GPU mode not active)')
    else:
        initialized = getattr(gpu_iface, 'initialized', False)
        gpu_dom = getattr(gpu_iface, 'gpu_dom', None)
        lines += [
            f'    status          : {"mapped to device" if initialized else "struct created, not yet mapped"}',
            f'    sizeof estimate : {_GPU_STRUCT_SIZE_EST} B  '
            f'(struct domain + MPI state + boundary sub-structs + flop counters)',
        ]

        # Estimated device memory (arrays mapped to GPU)
        try:
            from anuga.shallow_water.sw_domain_gpu_ext import (
                estimate_required_memory, gpu_available,
            )
            est_bytes = estimate_required_memory(N, nb)
            est_mb = est_bytes / 1024**2
            gpu_avail = gpu_available()
            lines += [
                f'    gpu available   : {gpu_avail}',
                f'    device memory   : {est_mb:.1f} MB  (estimate for N={N:,}, nb={nb:,})',
                '    note            : device arrays are OpenMP-mapped copies of numpy arrays',
            ]
        except ImportError:
            lines.append('    device memory   : (sw_domain_gpu_ext not importable)')

    lines.append(sep)
    return '\n'.join(lines)


def print_domain_struct_stats(domain):
    """Print a breakdown of the C domain struct and GPU domain struct memory.

    Convenience wrapper around :func:`domain_struct_stats`.

    Parameters
    ----------
    domain : anuga.Domain
        The simulation domain to inspect.

    Examples
    --------
    >>> import anuga
    >>> domain = anuga.rectangular_cross_domain(10, 10)
    >>> anuga.print_domain_struct_stats(domain)
    C / GPU struct diagnostics (N=400 tri)
    ...
    """
    print(domain_struct_stats(domain))


