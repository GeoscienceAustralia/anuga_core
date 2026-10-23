""" Classes to read an SWW file.
"""

import os
import numpy
import numpy as num
from anuga.utilities.file_utils import create_filename
from anuga.coordinate_transforms.geo_reference import \
    ensure_geo_reference
from .sts import Write_sts
from anuga.config import minimum_storable_height as default_minimum_storable_height
from anuga.config import institution as default_institution
from anuga.file.netcdf import NetCDFFile
import anuga.utilities.log as log
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.config import max_float
from anuga.config import netcdf_float, netcdf_float32, netcdf_int, netcdf_float64
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.coordinate_transforms.geo_reference import Geo_reference



class DataFileNotOpenError(Exception):
    pass


class SWWFileInUseError(Exception):
    """Another live process is writing this SWW file."""
    pass


class SWWTimeOrderError(Exception):
    """A frame about to be written is earlier than the last frame on file."""
    pass


# ---------------------------------------------------------------------------
# Writer lock
#
# The SWW path is set entirely by domain.set_name() / set_datadir(), and the
# writer opens the file afresh in append mode at every yieldstep, so two runs
# that happen to share a name both append to the same file and their frames
# interleave on the unlimited time dimension. The result opens fine and every
# frame holds plausible data; only the time variable, which is no longer
# monotonic, gives it away (anuga-community/anuga_core#232). The netCDF
# library holds no lock between those opens, so a sidecar '<file>.sww.lock'
# is taken when the file is created and held until the interpreter exits.
# ---------------------------------------------------------------------------

_held_sww_locks = set()


def sww_lock_path(filename):
    """The sidecar lock file for an SWW file."""
    return filename + '.lock'


def _release_all_sww_locks():
    for path in list(_held_sww_locks):
        release_sww_lock(path)


def _pid_is_alive(pid):
    """Whether a process with this id exists. Never signals or touches it."""
    if os.name == 'nt':
        # os.kill(pid, 0) is NOT a probe on Windows: it calls TerminateProcess
        # on a live process and raises a generic OSError on a dead one. Ask
        # the kernel for a query-only handle instead.
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL('kernel32', use_last_error=True)
        PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
        ERROR_ACCESS_DENIED = 5
        STILL_ACTIVE = 259
        handle = k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not handle:
            # No such process (ERROR_INVALID_PARAMETER), or one we may not
            # open, which nonetheless exists.
            return ctypes.get_last_error() == ERROR_ACCESS_DENIED
        try:
            code = wintypes.DWORD()
            if k32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return code.value == STILL_ACTIVE
            return True
        finally:
            k32.CloseHandle(handle)
    try:
        os.kill(pid, 0)       # signal 0: existence check only, on POSIX
    except ProcessLookupError:
        return False
    except PermissionError:
        return True           # exists, owned by someone else
    except OSError:
        return True           # cannot tell: treat as alive
    return True


def _read_sww_lock(path):
    """What a lock file records, as a dict. Missing keys are absent, and an
    unreadable file gives {}; `pid` is an int or absent."""
    fields = {}
    try:
        with open(path) as f:
            for line in f:
                key, sep, value = line.strip().partition('=')
                if sep and value:
                    fields[key] = value
    except OSError:
        return {}
    if 'pid' in fields:
        try:
            fields['pid'] = int(fields['pid'])
        except ValueError:
            del fields['pid']
    return fields


def _describe_sww_lock(fields, host):
    """Who holds a lock, for the error message: the process and host, plus
    the start time and script the lock records. A bare pid is little help
    when the other run is a notebook in another tab."""
    who = 'process %s on %s' % (fields.get('pid', '?'), fields.get('host') or host)
    started, script = fields.get('started'), fields.get('script')
    if started:
        who += ', started %s' % started
    if script:
        who += ', by %s' % script
    return who


def acquire_sww_lock(filename):
    """Take the writer lock for `filename`, or raise SWWFileInUseError.

    A lock left by a process that is no longer running (a crash, a kill)
    is taken over SILENTLY: the file is debris, and the run that left it
    can no longer be writing. Warning about it made an unrelated run fail
    under `warnings.simplefilter('error')`, which is a state a test suite
    reaches easily; the takeover is logged at DEBUG instead, so it is still
    on the record for anyone who turns that level on. A lock held by this very process is
    likewise taken over silently: re-running a model under the same name in
    one session is normal. A lock from another host cannot be checked for
    liveness and is respected; the message says how to clear it.
    """
    import atexit
    import socket
    import sys
    import time as _time
    import warnings

    path = sww_lock_path(filename)
    me = os.getpid()
    host = socket.gethostname()

    if os.path.exists(path):
        fields = _read_sww_lock(path)
        pid = fields.get('pid')
        other_host = fields.get('host')
        same_host = (other_host is None or other_host == host)
        if pid == me:
            pass                                    # our own earlier run
        elif same_host and pid is not None and not _pid_is_alive(pid):
            log.debug('Taking over the SWW writer lock %s left by process '
                      '%d, which is no longer running' % (path, pid))
        else:
            who = _describe_sww_lock(fields, host)
            raise SWWFileInUseError(
                'The SWW file %s is being written by another run (%s, lock '
                'file %s). Two runs appending to one SWW file interleave '
                'their frames and produce a file that is silently wrong. '
                'Give this run its own name with domain.set_name() or its '
                'own directory with domain.set_datadir(); if that run is '
                'known to be dead, delete the lock file.'
                % (filename, who, path))
        try:
            os.remove(path)
        except OSError:
            pass

    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, 'w') as f:
        f.write('pid=%d\nhost=%s\nstarted=%s\nscript=%s\n'
                % (me, host, _time.strftime('%Y-%m-%dT%H:%M:%S'),
                   sys.argv[0] if sys.argv else ''))
    if not _held_sww_locks:
        atexit.register(_release_all_sww_locks)
    _held_sww_locks.add(path)
    return path


def release_sww_lock(path):
    """Remove a lock this process holds. Safe to call twice."""
    _held_sww_locks.discard(path)
    pid = _read_sww_lock(path).get('pid')
    if pid == os.getpid():
        try:
            os.remove(path)
        except OSError:
            pass



class DataMissingValuesError(Exception):
    pass


class NewQuantity(Exception):
    pass


class DataDomainError(Exception):
    pass


class DataTimeError(Exception):
    pass


class Data_format:
    """Generic interface to data formats
    """

    def __init__(self, domain, extension, mode=netcdf_mode_w):
        assert mode[0] in ['r', 'w', 'a'], \
            "Mode %s must be either:\n" % mode + \
            "   'w' (write)\n" + \
            "   'r' (read)\n" + \
            "   'a' (append)"

        # Create filename
        self.filename = create_filename(domain.get_datadir(),
                                        domain.get_name(), extension)

        self.timestep = 0
        self.domain = domain

        # Probably should exclude ghosts in case this is a parallel domain

        self.number_of_nodes = domain.number_of_nodes
        self.number_of_volumes = domain.number_of_triangles
        #self.number_of_volumes = len(domain)

        # FIXME: Should we have a general set_precision function?


class SWW_file(Data_format):
    """Interface to native NetCDF format (.sww) for storing model output

    There are two kinds of data

    1: Constant data: Vertex coordinates and field values. Stored once
    2: Variable data: Conserved quantities. Stored once per timestep.

    All data is assumed to reside at vertex locations.
    """

    def __init__(self, domain,
                 mode=netcdf_mode_w, max_size=200000000000, recursion=False):

        self.precision = netcdf_float32  # Use single precision for quantities
        self.recursion = recursion
        self.mode = mode

        if hasattr(domain, 'max_size'):
            self.max_size = domain.max_size  # File size max is 2Gig
        else:
            self.max_size = max_size

        if hasattr(domain, 'store_centroids'):
            self.store_centroids = domain.store_centroids
        else:
            self.store_centroids = False

        if hasattr(domain, 'institution'):
            self.institution = domain.institution
        else:
            self.institution = default_institution

        if hasattr(domain, 'minimum_storable_height'):
            self.minimum_storable_height = domain.minimum_storable_height
        else:
            self.minimum_storable_height = default_minimum_storable_height

        if hasattr(domain, 'timezone'):
            self.timezone = str(domain.get_timezone())
        else:
            self.timezone = 'UTC'


        # Call parent constructor
        Data_format.__init__(self, domain, 'sww', mode)

        # Get static and dynamic quantities from domain
        static_quantities = []
        dynamic_quantities = []
        static_c_quantities = []
        dynamic_c_quantities = []
        self._overwrite_c_quantities = []  # flag=4: static 1D vars overwritten each yield step

        for q in domain.quantities_to_be_stored:
            flag = domain.quantities_to_be_stored[q]

            msg = 'Quantity %s is requested to be stored ' % q
            msg += 'but it does not exist in domain.quantities'
            assert q in domain.quantities, msg

            assert flag in [1, 2, 3, 4]
            if flag == 1:
                static_quantities.append(q)
                if self.store_centroids:
                    static_c_quantities.append(q+'_c')

            if flag == 2:
                dynamic_quantities.append(q)
                if self.store_centroids:
                    dynamic_c_quantities.append(q+'_c')

            if flag == 3:
                # centroid-only dynamic quantity: no vertex storage
                dynamic_c_quantities.append(q+'_c')

            if flag == 4:
                # centroid-only static variable overwritten at each yield step
                static_c_quantities.append(q+'_c')
                self._overwrite_c_quantities.append(q+'_c')

        # Passive tracers, and therefore suspended sediment classes, which are
        # tracers with settling parameters (see Domain.add_sediment_fraction).
        #
        # Tracers are NOT Quantity objects: they live in one contiguous
        # (n_tracers, N) block so the C kernel can stride them and the device
        # can map them in one go. They are deliberately kept out of
        # domain.quantities -- that dict is walked by the solver, not just by
        # this writer (mesh reordering rewrites vertex/edge/gradient arrays,
        # set_beta pushes a per-quantity beta where tracers share one scalar),
        # so registering a view there would put tracers in the path of
        # machinery that assumes a real Quantity.
        #
        # They are centroid quantities, so they take the same route as a flag-3
        # quantity: one dynamic <name>_c variable, no vertex storage, nothing
        # interpolated that the solver never computed.
        self._tracer_names = list(getattr(domain, '_tracer_names', []) or [])
        if self._tracer_names and getattr(domain, 'store_tracers', True):
            for name in self._tracer_names:
                cname = name + '_c'
                if cname not in dynamic_c_quantities:
                    dynamic_c_quantities.append(cname)
        else:
            self._tracer_names = []

        # Refuse to create a file another live run is writing (#232).
        self._lock = None
        if mode[0] == 'w':
            self._lock = acquire_sww_lock(self.filename)

        # NetCDF file definition
        fid = NetCDFFile(self.filename, mode)
        if mode[0] == 'w':
            description = 'Output from anuga.file.sww ' \
                          'suitable for plotting'

            self.writer = Write_sww(static_quantities,
                                    dynamic_quantities,
                                    static_c_quantities,
                                    dynamic_c_quantities,
                                    tracer_names=self._tracer_names)

            self.writer.store_header(fid,
                                     domain.starttime,
                                     self.number_of_volumes,
                                     self.domain.number_of_nodes,
                                     description=description,
                                     institution=self.institution,
                                     smoothing=domain.smooth,
                                     order=domain.default_order,
                                     sww_precision=self.precision,
                                     timezone=self.timezone)

            # Extra optional information
            if hasattr(domain, 'texture'):
                fid.texture = domain.texture

            if domain.quantities_to_be_monitored is not None:
                fid.createDimension('singleton', 1)
                fid.createDimension('two', 2)

                poly = domain.monitor_polygon
                if poly is not None:
                    N = len(poly)
                    fid.createDimension('polygon_length', N)
                    fid.createVariable('extrema.polygon',
                                       self.precision,
                                       ('polygon_length', 'two'))
                    fid.variables['extrema.polygon'][:] = poly

                interval = domain.monitor_time_interval
                if interval is not None:
                    fid.createVariable('extrema.time_interval',
                                       self.precision,
                                       ('two',))
                    fid.variables['extrema.time_interval'][:] = interval

                for q in domain.quantities_to_be_monitored:
                    fid.createVariable(q + '.extrema', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.min_location', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.max_location', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.min_time', self.precision,
                                       ('singleton',))
                    fid.createVariable(q + '.max_time', self.precision,
                                       ('singleton',))

        fid.close()

    def release(self):
        """Give up the writer lock. Called when the writer moves on to a
        continuation file; otherwise the lock is released at exit."""
        lock = getattr(self, '_lock', None)      # absent on writers unpickled
        if lock is not None:                     # from before the lock existed
            release_sww_lock(lock)
            self._lock = None

    def store_connectivity(self):
        """Store information about nodes, triangles and static quantities

        Writes x,y coordinates of triangles and their connectivity.

        Store also any quantity that has been identified as static.
        """

        # FIXME: Change name to reflect the fact thta this function
        # stores both connectivity (triangulation) and static quantities

        domain = self.domain

        # append to the NetCDF file
        fid = NetCDFFile(self.filename, netcdf_mode_a)

        # Get X, Y from one (any) of the quantities
        Q = list(domain.quantities.values())[0]
        X, Y, _, V = Q.get_vertex_values(xy=True, precision=self.precision)

        # store the connectivity data
        points = num.concatenate(
            (X[:, num.newaxis], Y[:, num.newaxis]), axis=1)
        self.writer.store_triangulation(fid,
                                        points,
                                        V.astype(num.float32),
                                        points_georeference=domain.geo_reference)

        if domain.parallel:
            self.writer.store_parallel_data(fid,
                                            domain.number_of_global_triangles,
                                            domain.number_of_global_nodes,
                                            domain.tri_full_flag,
                                            domain.tri_l2g,
                                            domain.node_l2g)

        # Get names of static quantities
        static_quantities = {}
        static_quantities_centroid = {}

        for name in self.writer.static_quantities:
            Q = domain.quantities[name]
            A, _ = Q.get_vertex_values(xy=False,
                                       precision=self.precision)
            static_quantities[name] = A

        # print domain.quantities
        # print self.writer.static_c_quantities

        for name in self.writer.static_c_quantities:
            Q = domain.quantities[name[:-2]]  # rip off _c from name
            static_quantities_centroid[name] = Q.centroid_values

        # Store static quantities
        self.writer.store_static_quantities(fid, **static_quantities)
        self.writer.store_static_quantities_centroid(
            fid, **static_quantities_centroid)

        fid.close()

    def store_timestep(self):
        """Store time and time dependent quantities
        """

        #import types
        from time import sleep
        from os import stat

        # Get NetCDF
        retries = 0
        file_open = False
        while not file_open and retries < 10:
            try:
                # Open existing file
                fid = NetCDFFile(self.filename, netcdf_mode_a)
            except OSError:
                # This could happen if someone was reading the file.
                # In that case, wait a while and try again
                msg = 'Warning (store_timestep): File %s could not be opened' \
                      % self.filename
                msg += ' - trying step %s again' % self.domain.relative_time
                log.info(msg)
                retries += 1
                sleep(1)
            else:
                file_open = True

        if not file_open:
            msg = 'File %s could not be opened for append' % self.filename
            raise DataFileNotOpenError(msg)

        # Check to see if the file is already too big:
        time = fid.variables['time'][:]

        # Frames must not go backwards. A frame at or before the last one on
        # file is legitimate only when it REWRITES an existing frame: a
        # checkpoint resume, or a second evolve() re-storing its starting
        # frame (Write_sww.store_quantities finds that slot). A time earlier
        # than the last frame that matches no frame at all is what a second
        # process appending to this file produces (#232); it used to be
        # appended with a log warning, leaving a file that opens fine, holds
        # plausible data and has a non-monotonic time variable.
        t_new = self.domain.relative_time
        if len(time) > 0 and t_new < time[-1]:
            tol = 1.0e-9 * max(1.0, abs(t_new))       # as store_quantities
            if not num.any(num.abs(time - t_new) <= tol):
                fid.close()
                raise SWWTimeOrderError(
                    'Refusing to append a frame at t=%g to %s, whose last '
                    'frame is at t=%g and which holds no frame at that time '
                    'to rewrite. Another process is writing this file (see '
                    'the lock file %s), or the model time went backwards.'
                    % (t_new, self.filename, time[-1],
                       sww_lock_path(self.filename)))

        i = len(time) + 1
        file_size = stat(self.filename)[6]
        file_size_increase = file_size//i
        if file_size + file_size_increase > self.max_size * 2**self.recursion:
            # In order to get the file name and start time correct,
            # I change the domain.filename and domain.starttime.
            # This is the only way to do this without changing
            # other modules (I think).

            # Write a filename addon that won't break the anuga viewers
            # (10.sww is bad)
            filename_ext = '_time_%s' % self.domain.relative_time
            filename_ext = filename_ext.replace('.', '_')

            # Remember the old filename, then give domain a
            # name with the extension
            old_domain_filename = self.domain.get_name()
            if not self.recursion:
                self.domain.set_name(old_domain_filename + filename_ext)

            # Temporarily change the domain starttime to the current time
            old_domain_starttime = self.domain.starttime
            self.domain.starttime = self.domain.get_time()

            # Build a new data_structure.
            next_data_structure = SWW_file(self.domain, mode=self.mode,
                                           max_size=self.max_size,
                                           recursion=self.recursion+1)
            if not self.recursion:
                log.info('    file_size = %s' % file_size)
                log.info('    saving file to %s'
                             % next_data_structure.filename)

            # Set up the new data_structure; this file is finished with
            self.domain.writer = next_data_structure
            self.release()

            # Store connectivity and first timestep
            next_data_structure.store_connectivity()
            next_data_structure.store_timestep()
            fid.sync()
            fid.close()

            # Restore the old starttime and filename
            self.domain.starttime = old_domain_starttime
            self.domain.set_name(old_domain_filename)
        else:
            self.recursion = False
            domain = self.domain

            # Get the variables
            time = fid.variables['time'][:]
            i = len(time)

            if 'stage' in self.writer.dynamic_quantities:
                # Select only those values for stage,
                # xmomentum and ymomentum (if stored) where
                # depth exceeds minimum_storable_height
                #
                # In this branch it is assumed that elevation
                # is also available as a quantity

                # Smoothing for the get_vertex_values will be obtained
                # from the smooth setting in domain

                Q = domain.quantities['stage']
                w, _ = Q.get_vertex_values(xy=False)

                Q = domain.quantities['elevation']
                z, _ = Q.get_vertex_values(xy=False)

                storable_indices = num.array(
                    w-z >= self.minimum_storable_height)

                # print numpy.sum(storable_indices), len(z), self.minimum_storable_height, numpy.min(w-z)
            else:
                # Very unlikely branch
                storable_indices = None  # This means take all

            # Now store dynamic quantities
            dynamic_quantities = {}
            dynamic_quantities_centroid = {}

            for name in self.writer.dynamic_quantities:
                #netcdf_array = fid.variables[name]

                Q = domain.quantities[name]
                A, _ = Q.get_vertex_values(xy=False,
                                           precision=self.precision)

                if storable_indices is not None:
                    if name == 'stage':
                        A = num.choose(storable_indices, (z, A))

                    if name in ['xmomentum', 'ymomentum']:
                        # Get xmomentum where depth exceeds
                        # minimum_storable_height

                        # Define a zero vector of same size and type as A
                        # for use with momenta
                        null = num.zeros(num.size(A), A.dtype.char)
                        A = num.choose(storable_indices, (null, A))

                dynamic_quantities[name] = A

            for name in self.writer.dynamic_c_quantities:
                base = name[:-2]
                if base in self._tracer_names:
                    # A tracer: a row of the (n_tracers, N) block, not a
                    # Quantity. The block is REALLOCATED by add_tracer, so it
                    # is looked up on the domain each step rather than cached.
                    s_idx = domain._tracer_names.index(base)
                    dynamic_quantities_centroid[name] = \
                        domain.tracer_centroid_values[s_idx]
                else:
                    Q = domain.quantities[base]
                    dynamic_quantities_centroid[name] = Q.centroid_values

            # Store dynamic quantities
            slice_index = self.writer.store_quantities(fid,
                                                       time=self.domain.relative_time,
                                                       sww_precision=self.precision,
                                                       **dynamic_quantities)

            # Store dynamic quantities
            if self.store_centroids or self.writer.dynamic_c_quantities:
                self.writer.store_quantities_centroid(fid,
                                                      slice_index=slice_index,
                                                      sww_precision=self.precision,
                                                      **dynamic_quantities_centroid)

            # Overwrite centroid-only static quantities (flag=4) each yield step
            for name in self._overwrite_c_quantities:
                Q = domain.quantities[name[:-2]]
                fid.variables[name][:] = Q.centroid_values.astype(self.precision)

            # Update extrema if requested
            domain = self.domain
            if domain.quantities_to_be_monitored is not None:
                for q, info in list(domain.quantities_to_be_monitored.items()):
                    if info['min'] is not None:
                        fid.variables[q + '.extrema'][0] = info['min']
                        fid.variables[q + '.min_location'][:] = \
                            info['min_location']
                        fid.variables[q + '.min_time'][0] = info['min_time']

                    if info['max'] is not None:
                        fid.variables[q + '.extrema'][1] = info['max']
                        fid.variables[q + '.max_location'][:] = \
                            info['max_location']
                        fid.variables[q + '.max_time'][0] = info['max_time']

            # Flush and close
            # fid.sync()
            fid.close()


class Read_sww:

    def __init__(self, source):
        """The source parameter is assumed to be a NetCDF sww file.
        """

        self.source = source

        self.frame_number = 0

        fin = NetCDFFile(self.source, 'r')

        self.time = num.array(fin.variables['time'][:], float)
        self.last_frame_number = self.time.shape[0] - 1

        self.frames = num.arange(self.last_frame_number+1)

        fin.close()

        self.read_mesh()

        self.quantities = {}

        self.read_quantities()

    def read_mesh(self):
        """ Read and store the mesh data contained within this sww file.
        """
        fin = NetCDFFile(self.source, 'r')

        self.vertices = num.array(fin.variables['volumes'][:], int)

        self.x = x = num.array(fin.variables['x'][:], float)
        self.y = y = num.array(fin.variables['y'][:], float)

        assert len(self.x) == len(self.y)

        self.xmin = num.min(x)
        self.xmax = num.max(x)
        self.ymin = num.min(y)
        self.ymax = num.max(y)

        fin.close()

    def read_quantities(self, frame_number=0):
        """
        Read the quantities contained in this file.
        frame_number is the time index to load.
        """
        assert frame_number >= 0 and frame_number <= self.last_frame_number

        self.frame_number = frame_number

        M = len(self.x)//3

        fin = NetCDFFile(self.source, 'r')

        for q in [n for n in list(fin.variables.keys()) if n != 'x' and n != 'y' and n != 'time' and n != 'volumes' and
                  '_range' not in n and '_c' not in n]:
            # print q
            if len(fin.variables[q].shape) == 1:  # Not a time-varying quantity
                self.quantities[q] = num.ravel(
                    num.array(fin.variables[q][:], float)).reshape(M, 3)
            else:  # Time-varying, get the current timestep data
                self.quantities[q] = num.array(
                    fin.variables[q][self.frame_number], float).reshape(M, 3)
        fin.close()
        return self.quantities

    def get_bounds(self):
        """
            Get the bounding rect around the mesh.
        """
        return [self.xmin, self.xmax, self.ymin, self.ymax]

    def get_last_frame_number(self):
        """
            Return the last loaded frame index.
        """
        return self.last_frame_number

    def get_time(self):
        """
            Get time at the current frame num, in secs.
        """
        return self.time[self.frame_number]


class Write_sww(Write_sts):
    """
        A class to write an SWW file.

        It is domain agnostic, and requires all the data to be fed in
        manually.
    """

    def __init__(self,
                 static_quantities,
                 dynamic_quantities,
                 static_c_quantities=None,
                 dynamic_c_quantities=None,
                 tracer_names=None):
        """Initialise Write_sww with two (or 4) list af quantity names:

        static_quantities (e.g. elevation or friction):
            Stored once at the beginning of the simulation in a 1D array
            of length number_of_points
        dynamic_quantities (e.g stage):
            Stored every timestep in a 2D array with
            dimensions number_of_points X number_of_timesteps

        static_c_quantities (e.g. elevation_c or friction_c):
            Stored once at the beginning of the simulation in a 1D array
            of length number_of_triangles
        dynamic_c_quantities (e.g stage_c):
            Stored every timestep in a 2D array with
            dimensions number_of_triangles X number_of_timesteps

        tracer_names (e.g. ['salinity', 'sand']):
            Which of the dynamic_c_quantities are tracers. Recorded on the
            file as a global attribute so a reader can tell `salinity_c`
            from `stage_c` -- both are <name>_c and nothing else
            distinguishes them. Written even when empty, so an old file
            (no attribute at all) is distinguishable from a new one that
            simply has no tracers.

        """
        self.static_quantities = static_quantities
        self.dynamic_quantities = dynamic_quantities
        self.static_c_quantities = static_c_quantities if static_c_quantities is not None else []
        self.dynamic_c_quantities = dynamic_c_quantities if dynamic_c_quantities is not None else []

        self.tracer_names = list(tracer_names) if tracer_names else []

        self.store_centroids = False
        if static_c_quantities or dynamic_c_quantities:
            self.store_centroids = True

    def store_header(self,
                     outfile,
                     times,
                     number_of_volumes,
                     number_of_points,
                     institution='Geosciences Australia',
                     description='Generated by ANUGA',
                     smoothing=True,
                     order=1,
                     sww_precision=netcdf_float32,
                     timezone='UTC',
                     verbose=False):
        """Write an SWW file header.

        Writes the first section of the .sww file.

        outfile - the open file that will be written
        times - A list of the time slice times OR a start time
        Note, if a list is given the info will be made relative.
        number_of_volumes - the number of triangles
        number_of_points - the number of vertices in the mesh
        """

        from anuga import get_revision_number
        from anuga import get_revision_date
        from anuga import get_version

        outfile.institution = institution
        outfile.description = description
        # Which dynamic centroid variables are tracers (space separated, and
        # '' when there are none). Lets a reader pick the tracers out without
        # having to keep a list of every quantity name that is not one.
        outfile.tracer_names = ' '.join(self.tracer_names)

        # For sww compatibility
        if smoothing is True:
            # Smoothing to be depreciated
            outfile.smoothing = 'Yes'
            outfile.vertices_are_stored_uniquely = 'False'
        else:
            # Smoothing to be depreciated
            outfile.smoothing = 'No'
            outfile.vertices_are_stored_uniquely = 'True'
        outfile.order = order

        try:
            revision_number = get_revision_number()
        except Exception:
            # This will be triggered if the system cannot get the
            # revision number.
            revision_number = None
        # Allow None to be stored as a string
        outfile.revision_number = str(revision_number)

        try:
            revision_date = get_revision_date()
        except Exception:
            # This will be triggered if the system cannot get the
            # revision date.
            revision_date = None
        # Allow None to be stored as a string
        outfile.revision_date = str(revision_date)

        try:
            anuga_version = get_version()
        except Exception:
            # This will be triggered if the system cannot get the
            # version.
            version = None
        # Allow None to be stored as a string
        outfile.anuga_version = str(anuga_version)

        # This is being used to seperate one number from a list.
        # what it is actually doing is sorting lists from numeric arrays.
        if isinstance(times, (list, num.ndarray)):
            number_of_times = len(times)
            times = ensure_numeric(times)
            if number_of_times == 0:
                starttime = 0
            else:
                starttime = times[0]
                times = times - starttime  # Store relative times
        else:
            number_of_times = 0
            starttime = times

        outfile.starttime = starttime
        outfile.timezone = timezone

        # dimension definitions
        outfile.createDimension('number_of_volumes', number_of_volumes)
        outfile.createDimension(
            'number_of_triangle_vertices', number_of_points)
        outfile.createDimension('number_of_vertices', 3)
        outfile.createDimension('numbers_in_range', 2)

        if smoothing is True:
            outfile.createDimension('number_of_points', number_of_points)
            # FIXME(Ole): This will cause sww files for parallel domains to
            # have ghost nodes stored (but not used by triangles).
            # To clean this up, we have to change get_vertex_values and
            # friends in quantity.py (but I can't be bothered right now)
        else:
            outfile.createDimension('number_of_points', 3*number_of_volumes)

        outfile.createDimension('number_of_timesteps', number_of_times)

        # variable definitions
        outfile.createVariable('x', sww_precision, ('number_of_points',))
        outfile.createVariable('y', sww_precision, ('number_of_points',))

        outfile.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                       'number_of_vertices'))

        for q in self.static_quantities:

            outfile.createVariable(q, sww_precision,
                                   ('number_of_points',))

            outfile.createVariable(q + Write_sww.RANGE, sww_precision,
                                   ('numbers_in_range',))

            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q+Write_sww.RANGE][0] = max_float  # Min
            outfile.variables[q+Write_sww.RANGE][1] = -max_float  # Max

        for q in self.static_c_quantities:
            outfile.createVariable(q, sww_precision,
                                   ('number_of_volumes',))

        self.write_dynamic_quantities(outfile, times, precis=sww_precision)

        outfile.sync()

    def store_triangulation(self,
                            outfile,
                            points_utm,
                            volumes,
                            zone=None,
                            new_origin=None,
                            points_georeference=None,
                            verbose=False):
        """
        Store triangulation data in the underlying file.

        Stores the points and triangle indices in the sww file

        outfile Open handle to underlying file.

        new_origin georeference that the points can be set to.

        points_georeference The georeference of the points_utm.

        verbose True if this function is to be verbose.

        new_origin - qa georeference that the points can be set to. (Maybe
        do this before calling this function.)

        points_utm - currently a list or array of the points in UTM.
        points_georeference - the georeference of the points_utm

        How about passing new_origin and current_origin.
        If you get both, do a convertion from the old to the new.

        If you only get new_origin, the points are absolute,
        convert to relative

        if you only get the current_origin the points are relative, store
        as relative.

        if you get no georefs create a new georef based on the minimums of
        points_utm.  (Another option would be to default to absolute)

        Yes, and this is done in another part of the code.
        Probably geospatial.

        If you don't supply either geo_refs, then supply a zone. If not
        the default zone will be used.

        precon:
            header has been called.
        """

        number_of_points = len(points_utm)
        volumes = num.array(volumes, num.int32).reshape(-1, 3)

        points_utm = num.array(points_utm)

        # Given the two geo_refs and the points, do the stuff
        # described in the method header
        # if this is needed else where, pull out as a function
        points_georeference = ensure_geo_reference(points_georeference)
        new_origin = ensure_geo_reference(new_origin)
        if new_origin is None and points_georeference is not None:
            points = points_utm
            geo_ref = points_georeference
        else:
            if new_origin is None:
                new_origin = Geo_reference(zone, min(points_utm[:, 0]),
                                           min(points_utm[:, 1]))
            points = new_origin.change_points_geo_ref(points_utm,
                                                      points_georeference)
            geo_ref = new_origin

        # At this stage I need a georef and points
        # the points are relative to the georef
        geo_ref.write_NetCDF(outfile)

        # This will put the geo ref in the middle
        #geo_ref = Geo_reference(refzone,(max(x)+min(x))/2.0,(max(x)+min(y))/2.)

        x = points[:, 0]
        y = points[:, 1]

        #x = x.astype(netcdf_float32)
        #y = y.astype(netcdf_float32)

        if verbose:
            log.info('------------------------------------------------')
            log.info('More Statistics:')
            log.info('  Extent (/lon):')
            log.info('    x in [%f, %f], len(lat) == %d'
                         % (min(x), max(x), len(x)))
            log.info('    y in [%f, %f], len(lon) == %d'
                         % (min(y), max(y), len(y)))
            # log.info('    z in [%f, %f], len(z) == %d'
            #             % (min(elevation), max(elevation), len(elevation)))
            log.info('geo_ref: %s' % str(geo_ref))
            log.info('------------------------------------------------')

        outfile.variables['x'][:] = x  # - geo_ref.get_xllcorner()
        outfile.variables['y'][:] = y  # - geo_ref.get_yllcorner()

        msg = 'Mismatch between shape of volumes array and (number_of_volumes , 3)'
        assert volumes.shape == outfile.variables['volumes'].shape, msg

        outfile.variables['volumes'][:] = volumes

    def write_dynamic_quantities(self, outfile,
                                 times, precis=netcdf_float32, verbose=False):
        """
            Write out given quantities to file.
        """

        for q in self.dynamic_quantities:
            outfile.createVariable(q, precis, ('number_of_timesteps',
                                               'number_of_points'))
            outfile.createVariable(q + Write_sts.RANGE, precis,
                                   ('numbers_in_range',))

            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q+Write_sts.RANGE][0] = max_float  # Min
            outfile.variables[q+Write_sts.RANGE][1] = -max_float  # Max

        for q in self.dynamic_c_quantities:
            outfile.createVariable(q, precis, ('number_of_timesteps',
                                               'number_of_volumes'))

        # Doing sts_precision instead of Float gives cast errors.
        outfile.createVariable('time', netcdf_float, ('number_of_timesteps',))

        if isinstance(times, (list, num.ndarray)):
            outfile.variables['time'][:] = times    # Store time relative

        if verbose:
            log.info('------------------------------------------------')
            log.info('Statistics:')
            log.info('    t in [%f, %f], len(t) == %d'
                         % (num.min(times), num.max(times), len(times.flat)))

    def store_parallel_data(self,
                            outfile,
                            number_of_global_triangles,
                            number_of_global_nodes,
                            tri_full_flag=None,
                            tri_l2g=None,
                            node_l2g=None,
                            sww_precision=netcdf_float32,
                            verbose=False):

        # dimension definitions
        #outfile.createDimension('number_of_volumes', number_of_volumes)
        #outfile.createDimension('number_of_vertices', 3)
        #outfile.createDimension('numbers_in_range', 2)

        # print 'store parallel data'
        outfile.number_of_global_triangles = number_of_global_triangles
        outfile.number_of_global_nodes = number_of_global_nodes

        # variable definitions
        outfile.createVariable('tri_l2g',  netcdf_int, ('number_of_volumes',))
        outfile.createVariable('node_l2g', netcdf_int,
                               ('number_of_triangle_vertices',))
        outfile.createVariable(
            'tri_full_flag', netcdf_int, ('number_of_volumes',))

        # print tri_l2g.shape
        # print tri_l2g
        # print outfile.variables['tri_l2g'].shape

        outfile.variables['tri_l2g'][:] = tri_l2g.astype(num.int32)

        # print node_l2g.shape
        # print node_l2g
        # print outfile.variables['node_l2g'].shape

        outfile.variables['node_l2g'][:] = node_l2g.astype(num.int32)

        # print tri_full_flag.shape
        # print tri_full_flag
        # print outfile.variables['tri_full_flag'].shape

        outfile.variables['tri_full_flag'][:] = tri_full_flag.astype(num.int32)

    def store_static_quantities(self,
                                outfile,
                                sww_precision=num.float32,
                                verbose=False,
                                **quant):
        """
        Write the static quantity info.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        The argument sww_precision allows for storing as either
        * single precision (default): num.float32
        * double precision: num.float64 or float

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be quantities from Domain such as friction
        #
        # Arrays not listed in static_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity,
        # e.g. stage_range, xmomentum_range and ymomentum_range
        for q in self.static_quantities:
            if q not in quant:
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity(msg)
            else:
                q_values = ensure_numeric(quant[q])

                x = q_values.astype(sww_precision)
                outfile.variables[q][:] = x

                # This populates the _range values
                outfile.variables[q + Write_sww.RANGE][0] = num.min(x)
                outfile.variables[q + Write_sww.RANGE][1] = num.max(x)

        # FIXME: Hack for backwards compatibility with old viewer
        # if 'elevation' in self.static_quantities:
        #    outfile.variables['z'][:] = outfile.variables['elevation'][:]

    def store_static_quantities_centroid(self,
                                         outfile,
                                         sww_precision=num.float32,
                                         verbose=False,
                                         **quant):
        """
        Write the static centroid quantity info.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        The argument sww_precision allows for storing as either
        * single precision (default): num.float32
        * double precision: num.float64 or float

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be quantities from Domain such as friction
        #
        # Arrays not listed in static_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity,
        # e.g. stage_range, xmomentum_range and ymomentum_range

        # print outfile.variables.keys()
        # print self.static_c_quantities

        for q in self.static_c_quantities:
            if q not in quant:
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity(msg)
            else:
                q_values = ensure_numeric(quant[q])

                x = q_values.astype(sww_precision)
                outfile.variables[q][:] = x

    def store_quantities(self,
                         outfile,
                         sww_precision=num.float32,
                         slice_index=None,
                         time=None,
                         verbose=False,
                         **quant):
        """
        Write the quantity info at each timestep.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        if the time array is already been built, use the slice_index
        to specify the index.

        Otherwise, use time to increase the time dimension

        Maybe make this general, but the viewer assumes these quantities,
        so maybe we don't want it general - unless the viewer is general

        The argument sww_precision allows for storing as either
        * single precision (default): num.float32
        * double precision: num.float64 or float

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        if time is not None:
            file_time = outfile.variables['time']
            slice_index = len(file_time)
            # check if time already saved as in check pointing
            if slice_index > 0:
                if time <= file_time[slice_index-1]:
                    # This time is being (re)written rather than appended
                    # (e.g. resuming from a checkpoint, or a near-duplicate
                    # yield time): find the existing slot to overwrite.
                    # The tolerance must scale with the time magnitude — a
                    # fixed 1e-14 is *below* the float64 ULP for any t beyond
                    # a few tens of seconds, so the match would spuriously
                    # fail at large t. If nothing matches, append rather than
                    # crash with an IndexError on an empty result.
                    tol = 1.0e-9 * max(1.0, abs(time))
                    check = numpy.where(numpy.abs(file_time[:]-time) <= tol)[0]
                    if len(check) > 0:
                        slice_index = int(check[0])
                    else:
                        log.warning(
                            'store_quantities: time %.17g <= last stored time '
                            '%.17g but no slot matches within %g; appending'
                            % (time, file_time[slice_index-1], tol))
                        slice_index = len(file_time)
            file_time[slice_index] = time
        else:
            # Has to be cast in case it was numpy.int
            slice_index = int(slice_index)

        # Write the named dynamic quantities
        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be the conserved quantities from Domain
        # (Typically stage,  xmomentum, ymomentum).
        #
        # Arrays not listed in dynamic_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity,
        # e.g. stage_range, xmomentum_range and ymomentum_range
        for q in self.dynamic_quantities:
            if q not in quant:
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity(msg)
            else:
                q_values = ensure_numeric(quant[q])

                q_retyped = q_values.astype(sww_precision)
                outfile.variables[q][slice_index] = q_retyped

                # This updates the _range values
                q_range = outfile.variables[q + Write_sww.RANGE][:]
                q_values_min = num.min(q_values)
                if q_values_min < q_range[0]:
                    outfile.variables[q + Write_sww.RANGE][0] = q_values_min
                q_values_max = num.max(q_values)
                if q_values_max > q_range[1]:
                    outfile.variables[q + Write_sww.RANGE][1] = q_values_max

        return slice_index

    def store_quantities_centroid(self,
                                  outfile,
                                  sww_precision=num.float32,
                                  slice_index=None,
                                  verbose=False,
                                  **quant):
        """
        Write the quantity centroid info at each timestep.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        if the time array is already been built, use the slice_index
        to specify the index.

        Otherwise, use time to increase the time dimension

        Maybe make this general, but the viewer assumes these quantities,
        so maybe we don't want it general - unless the viewer is general

        The argument sww_precision allows for storing as either
        * single precision (default): num.float32
        * double precision: num.float64 or float

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        assert slice_index is not None, 'slice_index should be set in store_quantities'

        # Write the named dynamic quantities
        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be the conserved quantities from Domain
        # (Typically stage,  xmomentum, ymomentum).
        #
        # Arrays not listed in dynamic_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity,
        # e.g. stage_range, xmomentum_range and ymomentum_range

        # print 50*"="
        # print quant
        # print self.dynamic_c_quantities

        for q in self.dynamic_c_quantities:
            if q not in quant:
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity(msg)
            else:
                q_values = ensure_numeric(quant[q])

                q_retyped = q_values.astype(sww_precision)
                outfile.variables[q][slice_index] = q_retyped

    def verbose_quantities(self, outfile):
        log.info('------------------------------------------------')
        log.info('More Statistics:')
        for q in self.dynamic_quantities:
            log.info('  %s in [%f, %f]'
                         % (q, outfile.variables[q+Write_sww.RANGE][0],
                            outfile.variables[q+Write_sww.RANGE][1]))
        log.info('------------------------------------------------')


def extent_sww(file_name):
    """Read in an sww file, then get its extents

    Input:
    file_name - the sww file

    Output:
    A list: [min(x),max(x),min(y),max(y),min(stage.flat),max(stage.flat)]
    """

    # Get NetCDF
    fid = NetCDFFile(file_name, netcdf_mode_r)

    # Get the variables
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    stage = fid.variables['stage'][:]

    fid.close()

    return [min(x), max(x), min(y), max(y), num.min(stage), num.max(stage)]


def Xload_sww_as_domain(filename, boundary=None, t=None,
                       fail_if_NaN=True, NaN_filler=0,
                       verbose=False, very_verbose=False):
    """
    Load an sww file into a domain.


    DEPRECATE THIS. It is not robust and doesn't for instance provide a proper
    boundary map (as this information isn't currently stored in sww files)


    Usage: domain = load_sww_as_domain('file.sww',
                        t=time (default = last time in file))

    Boundary is not recommended if domain.smooth is not selected, as it
    uses unique coordinates, but not unique boundaries. This means that
    the boundary file will not be compatable with the coordinates, and will
    give a different final boundary, or crash.
    """

    from anuga.shallow_water.shallow_water_domain import Domain

    # initialise NaN.
    NaN = 9.969209968386869e+036

    if verbose:
        log.info('Reading from %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)    # Open existing file for read
    time = fid.variables['time'][:]       # Timesteps
    if t is None:
        t = time[-1]
    time_interp = get_time_interp(time, t)

    # Get the variables as numeric arrays
    x = fid.variables['x'][:]                   # x-coordinates of vertices
    y = fid.variables['y'][:]                   # y-coordinates of vertices
    # elevation = fid.variables['elevation']      # Elevation
    # stage = fid.variables['stage']              # Water level
    # xmomentum = fid.variables['xmomentum']      # Momentum in the x-direction
    # ymomentum = fid.variables['ymomentum']      # Momentum in the y-direction

    starttime = float(fid.starttime)
    #starttime = fid.starttime[0]
    volumes = fid.variables['volumes'][:]       # Connectivity
    coordinates = num.transpose(num.asarray([x.tolist(), y.tolist()]))
    # FIXME (Ole): Something like this might be better:
    #                 concatenate((x, y), axis=1)
    # or              concatenate((x[:,num.newaxis], x[:,num.newaxis]), axis=1)

    dynamic_quantities = []
    interpolated_quantities = {}
    static_quantities = []

    # get geo_reference
    try:                             # sww files don't have to have a geo_ref
        geo_reference = Geo_reference(NetCDFObject=fid)
    except Exception:  # AttributeError, e:
        geo_reference = None

    if verbose:
        log.info('    getting quantities')

    for quantity in list(fid.variables.keys()):
        dimensions = fid.variables[quantity].dimensions
        if 'number_of_timesteps' in dimensions:
            dynamic_quantities.append(quantity)
            interpolated_quantities[quantity] = \
                interpolated_quantity(fid.variables[quantity][:], time_interp)
        else:
            static_quantities.append(quantity)

    # print static_quantities
    # print dynamic_quantities

    try:
        dynamic_quantities.remove('stage_c')
        dynamic_quantities.remove('xmomentum_c')
        dynamic_quantities.remove('ymomentum_c')
        dynamic_quantities.remove('elevation_c')
        dynamic_quantities.remove('friction_c')
    except ValueError:
        pass

    try:
        static_quantities.remove('elevation_c')
        static_quantities.remove('friction_c')
    except ValueError:
        pass

    static_quantities.remove('x')
    static_quantities.remove('y')
    # other_quantities.remove('z')
    static_quantities.remove('volumes')
    try:
        static_quantities.remove('stage_range')
        static_quantities.remove('xmomentum_range')
        static_quantities.remove('ymomentum_range')
        static_quantities.remove('elevation_range')
        static_quantities.remove('friction_range')
    except ValueError:
        pass

    dynamic_quantities.remove('time')

    if verbose:
        log.info('    building domain')

    #    From domain.Domain:
    #    domain = Domain(coordinates, volumes,\
    #                    conserved_quantities = conserved_quantities,\
    #                    other_quantities = other_quantities,zone=zone,\
    #                    xllcorner=xllcorner, yllcorner=yllcorner)

    # From shallow_water.Domain:
    coordinates = coordinates.tolist()
    volumes = volumes.tolist()
    # FIXME:should this be in mesh? (peter row)
    if fid.smoothing == 'Yes':
        unique = False
    else:
        unique = True

    if unique:
        coordinates, volumes, boundary = weed(coordinates, volumes, boundary)

    try:
        domain = Domain(coordinates, volumes, boundary,
                        starttime=(float(starttime) + float(t)))
    except AssertionError as e:
        fid.close()
        msg = 'Domain could not be created: %s. ' \
              'Perhaps use "fail_if_NaN=False and NaN_filler = ..."' % e
        raise DataDomainError(msg)

    if boundary is not None:
        domain.boundary = boundary

    domain.geo_reference = geo_reference

    for quantity in static_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except AttributeError:
            pass                       # quantity has no missing_value number
        X = fid.variables[quantity][:]
        if very_verbose:
            log.info('       %s' % str(quantity))
            log.info('        NaN = %s' % str(NaN))
            log.info('        max(X)')
            log.info('       %s' % str(max(X)))
            log.info('        max(X)==NaN')
            log.info('       %s' % str(max(X) == NaN))
            log.info('')
        if max(X) == NaN or min(X) == NaN:
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry' % quantity
                raise DataMissingValuesError(msg)
            else:
                data = (X != NaN)
                X = (X*data) + (data == 0)*NaN_filler
        if unique:
            X = num.resize(X, (len(X)//3, 3))
        domain.set_quantity(quantity, X)
    #
    for quantity in dynamic_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except AttributeError:
            pass                       # quantity has no missing_value number
        X = interpolated_quantities[quantity]
        if very_verbose:
            log.info('       %s' % str(quantity))
            log.info('        NaN = %s' % str(NaN))
            log.info('        max(X)')
            log.info('       %s' % str(max(X)))
            log.info('        max(X)==NaN')
            log.info('       %s' % str(max(X) == NaN))
            log.info('')
        if max(X) == NaN or min(X) == NaN:
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry' % quantity
                raise DataMissingValuesError(msg)
            else:
                data = (X != NaN)
                X = (X*data) + (data == 0)*NaN_filler
        if unique:
            X = num.resize(X, (X.shape[0]//3, 3))
        domain.set_quantity(quantity, X)

    fid.close()

    return domain


def get_mesh_and_quantities_from_file(filename,
                                      quantities=None,
                                      verbose=False):
    """Get and rebuild mesh structure and associated quantities from sww file

    Input:
        filename - Name os sww file
        quantities - Names of quantities to load

    Output:
        mesh - instance of class Interpolate
               (including mesh and interpolation functionality)
        quantities - arrays with quantity values at each mesh node or
                    each triangle vertex. (depending on whethr stored as smooth or not)
        time - vector of stored timesteps

    This function is used by e.g.:
        get_interpolated_quantities_at_polyline_midpoints
    """

    # FIXME (Ole): Maybe refactor filefunction using this more fundamental code.

    import types
    from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

    if verbose:
        log.info('Reading from %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)    # Open existing file for read
    time = fid.variables['time'][:]    # Time vector
    #time += fid.starttime[0]
    time += fid.starttime

    # Get the variables as numeric arrays
    x = fid.variables['x'][:]                   # x-coordinates of nodes
    y = fid.variables['y'][:]                   # y-coordinates of nodes

    elevation = fid.variables['elevation'][:]   # Elevation
    stage = fid.variables['stage'][:]           # Water level
    xmomentum = fid.variables['xmomentum'][:]   # Momentum in the x-direction
    ymomentum = fid.variables['ymomentum'][:]   # Momentum in the y-direction

    # Mesh (nodes (Mx2), triangles (Nx3))
    nodes = num.concatenate((x[:, num.newaxis], y[:, num.newaxis]), axis=1)
    triangles = fid.variables['volumes'][:]

    # Get geo_reference
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except Exception:  # AttributeError, e:
        # Sww files don't have to have a geo_ref
        geo_reference = None

    if verbose:
        log.info('    building mesh from sww file %s' % filename)

    boundary = None

    # FIXME (Peter Row): Should this be in mesh?
    if fid.smoothing != 'Yes':
        nodes = nodes.tolist()
        triangles = triangles.tolist()
        nodes, triangles, boundary = weed(nodes, triangles, boundary)

    try:
        mesh = Mesh(nodes, triangles, boundary, geo_reference=geo_reference)
    except AssertionError as e:
        fid.close()
        msg = 'Domain could not be created: %s. "' % e
        raise DataDomainError(msg)

    def gather(quantity):

        shape = quantity.shape

        def my_num_add_at(a, indices, b):
            """
            Use the numpy add.at operation if it is available, (numpy version >1.8)
            otherwise just use a quick and dirty implementation via a python loop
            """

            try:
                num.add.at(a, indices, b)
            except AttributeError:
                n_ids = len(indices)
                b_array = num.zeros_like(indices, dtype=float)
                b_array[:] = b

                for n in range(n_ids):
                    a[indices[n]] = a[indices[n]] + b_array[n]

        if len(shape) == 2:
            # time array
            if shape[1] == len(mesh.nodes):
                return quantity
            # need to calculate unique vertex values
            n_time = shape[0]
            n_nodes = len(mesh.nodes)

            mesh_ids = mesh.triangles.ravel()
            temp_uv = num.zeros((n_time, n_nodes))

            count_uv = num.zeros(len(mesh.nodes))
            my_num_add_at(count_uv, mesh_ids, 1)

            for i in range(n_time):
                my_num_add_at(temp_uv[i, :], mesh_ids, quantity[i, :])
                temp_uv[i, :] = temp_uv[i, :]/count_uv

        elif len(shape) == 1:
            # non time array
            if shape[0] == len(mesh.nodes):
                return quantity

            mesh_ids = mesh.triangles.ravel()
            temp_uv = num.zeros(len(mesh.nodes))

            count_uv = num.zeros(len(mesh.nodes))
            my_num_add_at(count_uv, mesh_ids, 1)
            my_num_add_at(temp_uv, mesh_ids, quantity)

        else:
            raise Exception

        return temp_uv

    quantities = {}

    if fid.smoothing != 'Yes':
        quantities['elevation'] = elevation[:]
        quantities['stage'] = stage[:]
        quantities['xmomentum'] = xmomentum[:]
        quantities['ymomentum'] = ymomentum[:]
    else:
        quantities['elevation'] = gather(elevation)
        quantities['stage'] = gather(stage)
        quantities['xmomentum'] = gather(xmomentum)
        quantities['ymomentum'] = gather(ymomentum)

    fid.close()

    return mesh, quantities, time


def get_time_interp(time, t=None):
    """Finds the ratio and index for time interpolation.
        time is an array of time steps
        t is the sample time.
        returns a tuple containing index into time, and ratio
    """
    if t is None:
        t = time[-1]
        index = -1
        ratio = 0.
    else:
        T = time
        tau = t
        index = 0
        msg = 'Time interval derived from file %s [%s:%s]' \
              % ('FIXMEfilename', T[0], T[-1])
        msg += ' does not match model time: %s' % tau
        if tau < time[0]:
            raise DataTimeError(msg)
        if tau > time[-1]:
            raise DataTimeError(msg)
        while tau > time[index]:
            index += 1
        while tau < time[index]:
            index -= 1
        if tau == time[index]:
            # Protect against case where tau == time[-1] (last time)
            # - also works in general when tau == time[i]
            ratio = 0
        else:
            # t is now between index and index+1
            ratio = (tau - time[index])/(time[index+1] - time[index])

    return (index, ratio)


def interpolated_quantity(saved_quantity, time_interp):
    """Interpolate a quantity with respect to time.

    saved_quantity  the quantity to interpolate
    time_interp     (index, ratio)

    Returns a vector of interpolated values.
    """

    index, ratio = time_interp

    Q = saved_quantity

    if ratio > 0:
        q = (1-ratio)*Q[index] + ratio*Q[index+1]
    else:
        q = Q[index]

    # Return vector of interpolated values
    return q


def weed(coordinates, volumes, boundary=None):
    """ Excise all duplicate points.
    """
    if isinstance(coordinates, num.ndarray):
        coordinates = coordinates.tolist()
    if isinstance(volumes, num.ndarray):
        volumes = volumes.tolist()

    unique = False
    point_dict = {}
    same_point = {}
    for i in range(len(coordinates)):
        point = tuple(coordinates[i])
        if point in point_dict:
            unique = True
            same_point[i] = point
            # to change all point i references to point j
        else:
            point_dict[point] = i
            same_point[i] = point

    coordinates = []
    i = 0
    for point in list(point_dict.keys()):
        point = tuple(point)
        coordinates.append(list(point))
        point_dict[point] = i
        i += 1

    for volume in volumes:
        for i in range(len(volume)):
            index = volume[i]
            if index > -1:
                volume[i] = point_dict[same_point[index]]

    new_boundary = {}
    if boundary is not None:
        for segment in list(boundary.keys()):
            point0 = point_dict[same_point[segment[0]]]
            point1 = point_dict[same_point[segment[1]]]
            label = boundary[segment]
            # FIXME should the bounday attributes be concaterated
            # ('exterior, pond') or replaced ('pond')(peter row)

            if (point0, point1) in new_boundary:
                new_boundary[(point0, point1)] = new_boundary[(point0, point1)]

            elif (point1, point0) in new_boundary:
                new_boundary[(point1, point0)] = new_boundary[(point1, point0)]
            else:
                new_boundary[(point0, point1)] = label

        boundary = new_boundary

    return coordinates, volumes, boundary


def sww_files_are_equal(filename1, filename2):
    """Read and compare numerical values of two sww files: filename1 and filename2

    If they are identical (up to a tolerance) the return value is True
    If anything substantial is different, the return value is False.
    """

    import anuga.utilities.plot_utils as util

    if not (filename1.endswith('.sww') and filename2.endswith('.sww')):
        msg = f'Filenames {filename1} and {filename2} must both end with .sww'
        raise Exception(msg)


    domain1_v = util.get_output(filename1)
    domain1_c = util.get_centroids(domain1_v)

    domain2_v = util.get_output(filename2)
    domain2_c = util.get_centroids(domain2_v)

    if not num.allclose(domain1_c.stage, domain2_c.stage):
        return False

    if not num.allclose(domain1_c.xmom, domain2_c.xmom):
        return False

    if not num.allclose(domain1_c.ymom, domain2_c.ymom):
        return False

    if not num.allclose(domain1_c.xvel, domain2_c.xvel):
        return False

    if not num.allclose(domain1_c.yvel, domain2_c.yvel):
        return False

    if not num.allclose(domain1_v.x, domain2_v.x):
        return False

    if not num.allclose(domain1_v.y, domain2_v.y):
        return False

    # Otherwise, they are deemed to be identical
    return True

