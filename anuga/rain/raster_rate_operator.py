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

import numpy as np
from anuga.operators.rate_operators import Rate_operator


class Raster_rate_operator(Rate_operator):
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
    time_offset : float, optional
        Offset added to ``domain.get_time()`` before looking up the raster
        time array.  Use this when raster times are absolute epoch seconds
        and the domain starts at 0.
    factor : float or callable, optional
        Multiplicative factor applied to all rates (passed to Rate_operator).
    polygon, region, indices, center, radius : optional
        Spatial restriction arguments passed to Rate_operator.
    default_rate : float, optional
    description, label, logging, verbose : optional
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

        self._raster        = raster_data
        self._time_offset   = time_offset
        self._last_slice_idx = -1

        # Pre-extract rates (m/s) at all centroids for every time slice.
        # extract_data_at_locations returns shape (n_slices, n_centroids).
        centroids = domain.centroid_coordinates
        depths = raster_data.extract_data_at_locations(centroids)  # (n_slices, N)
        self._rates_cache = depths / raster_data.time_step  # m/s

    def __call__(self):
        times = self._raster.times
        t     = self.domain.get_time() + self._time_offset
        idx   = int(np.searchsorted(times, t, side='right')) - 1
        idx   = max(0, min(idx, len(times) - 1))

        if idx != self._last_slice_idx:
            self.set_rate(self._rates_cache[idx])
            self._last_slice_idx = idx

        super().__call__()


class ARR_rate_operator(Rate_operator):
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
        centroid longitudes and latitudes.
    pattern : Single_pattern
        Temporal rainfall pattern.  ``pattern.Tplot`` (minutes) gives the
        time axis and ``pattern.Rplot`` (mm) gives the incremental depth
        delivered in each ``Tstep``-minute interval.
    factor : float or callable, optional
    polygon, region, indices, center, radius : optional
        Spatial restriction arguments passed to Rate_operator.
    default_rate : float, optional
    description, label, logging, verbose : optional

    Notes
    -----
    The rate at centroid *c* during pattern step *i* is

    .. math::

        r_{c,i} = \\frac{D_c \\cdot R_i}{E_{\\text{dep}} \\cdot \\Delta t \\cdot 1000}
        \\quad (\\text{m s}^{-1})

    where :math:`D_c` is ``depth_array[c]`` (mm), :math:`R_i` is
    ``pattern.Rplot[i]`` (mm), :math:`E_{\\text{dep}}` is ``pattern.Ev_dep``
    (mm), and :math:`\\Delta t` is ``pattern.Tstep`` in seconds.
    """

    def __init__(self, domain, depth_array, pattern,
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

        depth_array = np.asarray(depth_array, dtype=float)
        assert depth_array.shape == (N,), (
            f'depth_array shape {depth_array.shape} must be ({N},)'
        )

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
        # t=0 correctly selects Rplot[1] (the first non-zero step), and
        # the last increment is applied right through to finaltime.
        idx = int(np.searchsorted(self._times_sec[1:], t, side='right')) + 1
        idx = min(idx, len(self._times_sec) - 1)

        if idx != self._last_slice_idx:
            self.set_rate(self._rates_cache[idx])
            self._last_slice_idx = idx

        super().__call__()
