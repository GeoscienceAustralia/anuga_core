"""Read an SWW file's centroid quantities one timestep at a time.

SWW_plotter used to read stage, xmomentum and ymomentum for every timestep
and then derive depth, velocity, speed and speed x depth for every timestep
too: about eight arrays the size of the whole run, held in memory for as
long as the plotter lives.  A 66k-triangle, 2000-step run needs several GB,
and anuga_sww_gui builds one plotter per animation worker.

FrameArray stands in for those (time, triangle) arrays.  It reads what each
access needs from the file -- one timestep for ``a[frame]``, a few
triangles' timeseries for ``a[:, triangles]`` -- and streams over the
timesteps for ``min``, ``max`` and ``max(axis=0)``.  Anything else (numpy
functions, arithmetic, other methods) gets the whole array, so existing
code keeps working, at the old memory cost.
"""

import collections
import os

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin


class SWWSource:
    """An SWW file's variables, read a timestep or a few triangles at a time.

    The file stays open between reads, and is reopened after ``close()``
    or in a process forked from the one that opened it.
    """

    # Variable timesteps kept, so the quantities derived from the same
    # frame (depth, speed, ...) do not each read it again.
    CACHED_READS = 16

    def __init__(self, filename):
        self.filename = filename
        self._fid = None
        self._pid = None
        self._cache = collections.OrderedDict()

    def _variable(self, name):
        if self._fid is None or self._pid != os.getpid():
            from anuga.file.netcdf import NetCDFFile
            self._fid = NetCDFFile(self.filename)
            self._pid = os.getpid()
            self._cache.clear()
        return self._fid.variables[name]

    def close(self):
        if self._fid is not None and self._pid == os.getpid():
            self._fid.close()
        self._fid = None
        self._cache.clear()

    def __getstate__(self):
        # An open netCDF handle cannot be pickled; the copy opens its own.
        return {'filename': self.filename}

    def __setstate__(self, state):
        self.__init__(state['filename'])

    def is_animated(self, name):
        return self._variable(name).ndim == 2

    def frame(self, name, index):
        """The variable at one timestep (all of it, if it is static)."""
        v = self._variable(name)
        key = (name, index if v.ndim == 2 else None)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        data = v[index, :] if v.ndim == 2 else v[:]
        self._cache[key] = data
        if len(self._cache) > self.CACHED_READS:
            self._cache.popitem(last=False)
        return data

    def columns(self, name, triangles):
        """The variable over time for some triangles: shape (time,) + the
        shape of triangles, or just that shape if the variable is static."""
        v = self._variable(name)
        if isinstance(triangles, slice):
            return v[:, triangles] if v.ndim == 2 else v[triangles]
        triangles = np.asarray(triangles)
        if triangles.dtype == bool:
            triangles = np.flatnonzero(triangles)
        flat = triangles.ravel() % v.shape[-1]
        # netCDF wants increasing indices without repeats.
        unique, inverse = np.unique(flat, return_inverse=True)
        if v.ndim == 1:
            return v[unique][inverse].reshape(triangles.shape)
        data = v[:, unique][:, inverse]
        return data.reshape((v.shape[0],) + triangles.shape)


def _stack(frames):
    if any(isinstance(f, np.ma.MaskedArray) for f in frames):
        return np.ma.stack(frames)
    return np.stack(frames)


class FrameArray(NDArrayOperatorsMixin):
    """A (time, triangle) array computed from an SWW file as it is used.

    compute(get) returns the quantity from get(name), the named file
    variable for the requested selection: one timestep, or some triangles
    over time.  It must work elementwise, so the same function serves both,
    and involve at least one variable that changes with time.
    """

    # Computed timesteps kept, for the frame drawn now and then again.
    CACHED_FRAMES = 2

    def __init__(self, source, compute, ntime, ntriangles):
        self._source = source
        self._compute = compute
        self.shape = (ntime, ntriangles)
        self.ndim = 2
        self._frames = collections.OrderedDict()
        self._reductions = {}

    def __len__(self):
        return self.shape[0]

    @property
    def size(self):
        return self.shape[0] * self.shape[1]

    @property
    def dtype(self):
        return np.asarray(self._frame(0)).dtype

    def _frame(self, index):
        index = range(self.shape[0])[index]          # negative, and bounds
        if index in self._frames:
            self._frames.move_to_end(index)
            return self._frames[index]
        data = self._compute(lambda name: self._source.frame(name, index))
        self._frames[index] = data
        if len(self._frames) > self.CACHED_FRAMES:
            self._frames.popitem(last=False)
        return data

    def _columns(self, triangles):
        return self._compute(lambda name: self._source.columns(name, triangles))

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if len(key) > 2:
            raise IndexError('too many indices for a (time, triangle) array')
        when = key[0]
        where = key[1] if len(key) == 2 else slice(None)
        if isinstance(when, (int, np.integer)):
            data = self._frame(int(when))
            return data if _is_everything(where) else data[where]
        if _is_everything(when) and not _is_everything(where):
            return self._columns(where)
        steps = np.arange(self.shape[0])[when]
        return _stack([self._frame(int(i)) for i in steps])[:, where]

    def __iter__(self):
        for index in range(self.shape[0]):
            yield self._frame(index)

    def __array__(self, dtype=None, copy=None):
        data = np.asarray(_stack(list(self)))
        return data if dtype is None else data.astype(dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        inputs = [np.asarray(i) if isinstance(i, FrameArray) else i for i in inputs]
        return getattr(ufunc, method)(*inputs, **kwargs)

    def _reduce(self, name, combine, axis):
        if axis is not None and axis < 0:
            axis += 2
        if axis == 1:
            return np.array([getattr(f, name)() for f in self])
        key = (name, axis)
        if key not in self._reductions:
            result = None
            for frame in self:
                result = frame.copy() if result is None else combine(result, frame)
            self._reductions[key] = result if axis == 0 else getattr(result, name)()
        result = self._reductions[key]
        return result.copy() if axis == 0 else result

    def max(self, axis=None, out=None, keepdims=False, **kwargs):
        """The largest value, or with axis=0 each triangle's largest over
        time, found a timestep at a time."""
        if out is not None or keepdims or kwargs:
            return np.asarray(self).max(axis=axis, out=out, keepdims=keepdims, **kwargs)
        return self._reduce('max', np.maximum, axis)

    def min(self, axis=None, out=None, keepdims=False, **kwargs):
        """The smallest value, or with axis=0 each triangle's smallest over
        time, found a timestep at a time."""
        if out is not None or keepdims or kwargs:
            return np.asarray(self).min(axis=axis, out=out, keepdims=keepdims, **kwargs)
        return self._reduce('min', np.minimum, axis)

    def __getattr__(self, name):
        # Other ndarray methods (mean, copy, astype, ...) on the whole array.
        if name.startswith('_'):
            raise AttributeError(name)
        return getattr(np.asarray(self), name)

    def __repr__(self):
        return 'FrameArray(shape=%s, file=%r)' % (self.shape, self._source.filename)


def _is_everything(index):
    return isinstance(index, slice) and index == slice(None)
