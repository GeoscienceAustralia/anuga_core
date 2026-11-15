"""
Rate operators driven by raster or ARR 2016 rainfall data.

Classes
-------
Raster_rate_operator
    Rate_operator backed by a Raster_time_slice_data instance.
    Centroid rates are pre-interpolated at initialisation; each
    yield-step only does an array lookup, keeping the GPU path free
    of file I/O.

ARR_rate_operator
    Rate_operator driven by an ARR 2016 temporal pattern and a
    per-centroid total event depth array.  Like Raster_rate_operator
    the per-centroid rates for every pattern time step are precomputed
    at initialisation.
"""

import warnings

import numpy as np
from anuga.operators.rate_operators import Rate_operator


class _Precomputed_rate_operator(Rate_operator):
    """Rate_operator switching between precomputed centroid rate arrays.

    Outside the time span of the data the rate is ``default_rate``, kept as a
    centroid array so the operator stays on the ``centroid_array`` GPU path.
    """

    def _set_slice(self, idx):
        """Apply slice ``idx``, or ``default_rate`` when ``idx`` is None."""
        if idx is not None:
            if idx != self._last_slice_idx:
                self.set_rate(self._rates_cache[idx])
                self._last_slice_idx = idx
            return

        default = self.default_rate
        value = 0.0 if default is None else float(default(self.domain.get_time()))
        if self._last_slice_idx != -2 or value != self._last_default_value:
            self.set_rate(np.full(self.domain.number_of_elements, value))
            self._last_slice_idx = -2
            self._last_default_value = value


class Raster_rate_operator(_Precomputed_rate_operator):
    """Rate_operator backed by a :class:`~anuga.rain.Raster_time_slice_data`.

    The raster data is interpolated to all domain centroids once during
    :meth:`__init__`; after that each call only performs an array index
    lookup and, when the time slice changes, a :meth:`~Rate_operator.set_rate`
    call that keeps the ``centroid_array`` rate type required by the GPU path.

    Parameters
    ----------
    domain : anuga Domain
    raster_data : Raster_time_slice_data
        Object whose ``times`` (s), ``time_step`` (s), ``x``, ``y``, and
        ``data_slices`` attributes describe the gridded rain data.
        Data slices are assumed to be in *metres per time_step* (SI).
        Slice *i* covers ``times[i]`` to ``times[i] + time_step``.
    time_offset : float, optional
        Offset added to ``domain.get_time()`` before looking up the raster
        time array.  Use this when raster times are absolute epoch seconds
        and the domain starts at 0.
    factor : float or callable, optional
        Multiplicative factor applied to all rates (passed to Rate_operator).
    polygon, region, indices, center, radius : optional
        Spatial restriction arguments passed to Rate_operator.
    default_rate : float or callable(t), optional
        Rate (m/s) outside the time span of the data.  Default 0.
    description, label, logging, verbose : optional

    Attributes
    ----------
    time_step : float
        Time (s) between changes of rate, ``raster_data.time_step``. Use it, or a
        whole fraction of it, as the ``yieldstep`` of
        :meth:`~anuga.Domain.evolve` to deliver the exact depth.
    """

    def __init__(self, domain, raster_data,
                 time_offset=0.0,
                 factor=1.0,
                 region=None, indices=None, polygon=None,
                 center=None, radius=None,
                 default_rate=0.0,
                 description=None, label=None,
                 logging=False, verbose=False):

        N = domain.number_of_elements
        super().__init__(domain,
                         rate=np.zeros(N),
                         factor=factor,
                         region=region,
                         indices=indices,
                         polygon=polygon,
                         center=center,
                         radius=radius,
                         default_rate=default_rate,
                         description=description,
                         label=label,
                         logging=logging,
                         verbose=verbose)

        self.time_step      = raster_data.time_step
        self._raster        = raster_data
        self._time_offset   = time_offset
        self._last_slice_idx = -1

        # Pre-extract rates (m/s) at all centroids for every time slice.
        # extract_data_at_locations returns shape (n_slices, n_centroids).
        centroids = domain.get_centroid_coordinates(absolute=True)
        depths = raster_data.extract_data_at_locations(centroids)  # (n_slices, N)
        self._rates_cache = depths / raster_data.time_step  # m/s

    def __call__(self):
        times = self._raster.times
        t     = self.domain.get_time() + self._time_offset
        if times[0] <= t < times[-1] + self._raster.time_step:
            self._set_slice(int(np.searchsorted(times, t, side='right')) - 1)
        else:
            self._set_slice(None)

        super().__call__()


class ARR_rate_operator(_Precomputed_rate_operator):
    """Rate_operator driven by ARR 2016 spatial depth and temporal pattern.

    The per-centroid rates for every pattern time step are precomputed
    at :meth:`__init__`.  During the run only an array lookup occurs each
    yield step, which keeps the GPU execution path free of file I/O.

    Parameters
    ----------
    domain : anuga Domain
    depth_array : array-like, shape (n_elements,)
        Total event depth (mm) at each domain centroid.  Typically obtained
        from :meth:`~anuga.rain.Arr_grd.get_rain_at_points` called with the
        centroid longitudes and latitudes.  Values on the cells the operator
        acts on must be finite and non-negative.
    pattern : Single_pattern
        Temporal rainfall pattern.  ``pattern.Tplot`` (minutes) gives the
        time axis and ``pattern.Rplot`` (mm) gives the incremental depth
        delivered in each ``Tstep``-minute interval.
    factor : float or callable, optional
    polygon, region, indices, center, radius : optional
        Spatial restriction arguments passed to Rate_operator.
    default_rate : float or callable(t), optional
        Rate (m/s) before the pattern starts and after it ends.  Default 0.
    description, label, logging, verbose : optional
    ifd_duration : int, optional
        Duration (minutes) of the IFD grid behind ``depth_array``, for example
        ``grd.dur``.  A warning is given if it differs from ``pattern.Ev_dur``.

    Attributes
    ----------
    time_step : float
        Time (s) between changes of rate, ``pattern.Tstep`` in seconds. Use it, or
        a whole fraction of it, as the ``yieldstep`` of
        :meth:`~anuga.Domain.evolve` to deliver the exact depth.

    Notes
    -----
    The rate at centroid *c* during pattern step *i* is

    .. math::

        r_{c,i} = \\frac{D_c \\cdot R_i}{E_{\\text{dep}} \\cdot \\Delta t \\cdot 1000}
        \\quad (\\text{m s}^{-1})

    where :math:`D_c` is ``depth_array[c]`` (mm), :math:`R_i` is
    ``pattern.Rplot[i]`` (mm), :math:`E_{\\text{dep}}` is ``pattern.Ev_dep``
    (mm), and :math:`\\Delta t` is ``pattern.Tstep`` in seconds.

    ``pattern.Rplot`` is itself scaled by ``Ev_dep``, so ``Ev_dep`` cancels:
    the depth comes from ``depth_array`` alone and any positive ``Ev_dep``
    gives the same rates.
    """

    def __init__(self, domain, depth_array, pattern,
                 factor=1.0,
                 region=None, indices=None, polygon=None,
                 center=None, radius=None,
                 default_rate=0.0,
                 description=None, label=None,
                 logging=False, verbose=False,
                 ifd_duration=None):

        N = domain.number_of_elements
        super().__init__(domain,
                         rate=np.zeros(N),
                         factor=factor,
                         region=region,
                         indices=indices,
                         polygon=polygon,
                         center=center,
                         radius=radius,
                         default_rate=default_rate,
                         description=description,
                         label=label,
                         logging=logging,
                         verbose=verbose)

        depth_array = np.asarray(depth_array, dtype=float)
        assert depth_array.shape == (N,), (
            f'depth_array shape {depth_array.shape} must be ({N},)'
        )

        used = depth_array if self.indices is None else depth_array[self.indices]
        bad = ~np.isfinite(used) | (used < 0.0)
        if bad.any():
            raise ValueError(
                f'depth_array has {int(bad.sum())} NaN, infinite or negative '
                f'values on the cells this operator acts on (NODATA in the IFD '
                f'grid?)'
            )
        # Cells outside the operator may hold NODATA: keep it out of the rates.
        depth_array = np.where(np.isfinite(depth_array), depth_array, 0.0)

        if ifd_duration is not None and int(ifd_duration) != int(pattern.Ev_dur):
            warnings.warn(
                f'ARR_rate_operator: the IFD grid is for {int(ifd_duration)} min '
                f'but the temporal pattern is for {int(pattern.Ev_dur)} min.',
                stacklevel=2)

        self.time_step       = pattern.Tstep * 60.0
        self._pattern        = pattern
        self._last_slice_idx = -1

        # Pattern time axis converted to seconds (Tplot is in minutes)
        self._times_sec = pattern.Tplot * 60.0

        Ev_dep    = pattern.Ev_dep          # mm
        Tstep_sec = pattern.Tstep * 60.0    # s

        # Pre-compute rate (m/s) for every pattern step at every centroid.
        # rate[i, c] = depth_c * (Rplot[i] / Ev_dep) / Tstep_sec / 1000
        n_steps = len(pattern.Rplot)
        self._rates_cache = np.zeros((n_steps, N), dtype=float)
        if Ev_dep > 0.0 and Tstep_sec > 0.0:
            scale = depth_array / (Ev_dep * Tstep_sec * 1000.0)
            for i, r in enumerate(pattern.Rplot):
                self._rates_cache[i] = scale * r

    def __call__(self):
        t = self.domain.get_time()
        # Rplot[i] is the depth for the interval ending at _times_sec[i].
        # Search against the interval end-times (index 1 onward) so that
        # t=0 correctly selects Rplot[1] (the first non-zero step).
        if self._times_sec[0] <= t < self._times_sec[-1]:
            self._set_slice(
                int(np.searchsorted(self._times_sec[1:], t, side='right')) + 1)
        else:
            self._set_slice(None)

        super().__call__()
