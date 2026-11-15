"""
Rate operators driven by raster or ARR 2016 rainfall data.

Classes
-------
Hyetograph_rate_operator
    Rain that is uniform in space and changes in steps through time.

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

    losses = None
    cumulative_loss = 0.0
    _loss_slice_idx = -1      # last slice given to the loss model
    _loss_discharge = 0.0     # m^3/s lost on this operator's cells, this slice

    def _gross_rate(self, idx):
        """Rain rate (m/s) at every centroid in slice ``idx``."""
        return self._rates_cache[idx]

    def _slice_duration(self, idx):
        """Length (s) of slice ``idx``."""
        return self.time_step

    def _set_losses(self, losses):
        """Attach a rainfall loss model (see :mod:`anuga.rain.losses`)."""
        self.losses = losses
        if losses is not None:
            losses.bind(self.domain)

    def _excess_rate(self, idx):
        """Rate (m/s) of slice ``idx`` once the losses are taken out.

        The loss model keeps the state of the storm, so it is given every
        slice up to ``idx`` once and in order, including any the run stepped
        over.
        """
        if idx <= self._loss_slice_idx:
            raise RuntimeError(
                'rainfall losses follow the storm forwards in time; the '
                'operator cannot go back to an earlier rain time slice')
        for i in range(self._loss_slice_idx + 1, idx + 1):
            gross = self._gross_rate(i)
            dt = self._slice_duration(i)
            rate = self.losses.excess(gross*dt, dt)/dt
        self._loss_slice_idx = idx

        lost = gross - rate
        if self.indices is not None:
            lost = lost[self.indices]
        self._loss_discharge = float(
            (lost*self.areas)[self.full_indices].sum())
        return rate

    def _set_slice(self, idx):
        """Apply slice ``idx``, or ``default_rate`` when ``idx`` is None."""
        if idx is not None:
            if idx != self._last_slice_idx:
                if self.losses is None:
                    self.set_rate(self._gross_rate(idx))
                else:
                    self.set_rate(self._excess_rate(idx))
                self._last_slice_idx = idx
            return

        self._loss_discharge = 0.0

        default = self.default_rate
        value = 0.0 if default is None else float(default(self.domain.get_time()))
        if self._last_slice_idx != -2 or value != self._last_default_value:
            self.set_rate(np.full(self.domain.number_of_elements, value))
            self._last_slice_idx = -2
            self._last_default_value = value

    def _count_loss(self):
        """Add this step's lost rain to ``cumulative_loss`` (m^3)."""
        if self.losses is not None and self._loss_discharge != 0.0:
            self.cumulative_loss += (self.get_factor(self.domain.get_time())
                                     * self.domain.get_timestep()
                                     * self._loss_discharge)


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
    losses : Rainfall_loss, optional
        A rainfall loss model, such as :class:`~anuga.rain.SCS_curve_number`
        or :class:`~anuga.rain.Initial_continuing_loss`.  The operator then
        applies the excess rain only.  ``factor`` scales the excess, and
        ``default_rate`` is applied with no losses.

    Attributes
    ----------
    time_step : float
        Time (s) between changes of rate, ``raster_data.time_step``. Use it, or a
        whole fraction of it, as the ``yieldstep`` of
        :meth:`~anuga.Domain.evolve` to deliver the exact depth.
    cumulative_loss : float
        Volume (m^3) of rain lost so far when ``losses`` is given.
        ``cumulative_influx`` is the water added, which is the excess.
    """

    def __init__(self, domain, raster_data,
                 time_offset=0.0,
                 factor=1.0,
                 region=None, indices=None, polygon=None,
                 center=None, radius=None,
                 default_rate=0.0,
                 description=None, label=None,
                 logging=False, verbose=False,
                 losses=None):

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
        self._set_losses(losses)

    def __call__(self):
        times = self._raster.times
        t     = self.domain.get_time() + self._time_offset
        if times[0] <= t < times[-1] + self._raster.time_step:
            self._set_slice(int(np.searchsorted(times, t, side='right')) - 1)
        else:
            self._set_slice(None)

        super().__call__()
        self._count_loss()


class Hyetograph_rate_operator(_Precomputed_rate_operator):
    """Rain that is the same everywhere and changes in steps through time.

    This is the usual rain-on-mesh input: a hyetograph from a gauge or a
    design storm, applied over the whole domain or over a region, with
    losses that can vary from cell to cell.

    Parameters
    ----------
    domain : anuga Domain
    times : array-like, length n + 1
        Start of each interval and the end of the last (s), increasing.
        Interval *i* runs from ``times[i]`` to ``times[i+1]``.
    depths : array-like, length n, optional
        Rain (mm) falling in each interval.
    intensities : array-like, length n, optional
        Rain intensity (mm/h) in each interval. Give this or ``depths``.
    losses : Rainfall_loss, optional
        A rainfall loss model, such as :class:`~anuga.rain.SCS_curve_number`
        or :class:`~anuga.rain.Initial_continuing_loss`.  The operator then
        applies the excess rain only.  ``factor`` scales the excess, and
        ``default_rate`` is applied with no losses.
    factor : float or callable, optional
    polygon, region, indices, center, radius : optional
        Spatial restriction arguments passed to Rate_operator.
    default_rate : float or callable(t), optional
        Rate (m/s) before ``times[0]`` and from ``times[-1]`` on.  Default 0.
    description, label, logging, verbose : optional

    Attributes
    ----------
    time_step : float or None
        Time (s) between changes of rate when the intervals are all the same
        length, otherwise None. End each ``yieldstep`` of
        :meth:`~anuga.Domain.evolve` on a change to deliver the exact depth.
    cumulative_loss : float
        Volume (m^3) of rain lost so far when ``losses`` is given.
        ``cumulative_influx`` is the water added, which is the excess.
    """

    def __init__(self, domain, times, depths=None, intensities=None,
                 losses=None,
                 factor=1.0,
                 region=None, indices=None, polygon=None,
                 center=None, radius=None,
                 default_rate=0.0,
                 description=None, label=None,
                 logging=False, verbose=False):

        times = np.asarray(times, dtype=float)
        if (depths is None) == (intensities is None):
            raise ValueError('give either depths (mm) or intensities (mm/h)')
        if times.ndim != 1 or len(times) < 2 or not (np.diff(times) > 0.0).all():
            raise ValueError('times must be increasing, with at least two values')
        durations = np.diff(times)
        if depths is not None:
            rates = np.asarray(depths, dtype=float) * 1.0e-3 / durations
        else:
            rates = np.asarray(intensities, dtype=float) * 1.0e-3 / 3600.0
        if rates.shape != durations.shape:
            raise ValueError(
                f'{len(times)} times mark {len(durations)} intervals, but '
                f'{rates.size} rain values were given')
        if not np.isfinite(rates).all() or (rates < 0.0).any():
            raise ValueError('rain values must be finite and not negative')

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

        self._times = times
        self._durations = durations
        self._rates = rates                     # m/s, one per interval
        self.time_step = (float(durations[0])
                          if np.allclose(durations, durations[0]) else None)
        self._last_slice_idx = -1
        self._set_losses(losses)

    def _gross_rate(self, idx):
        return np.full(self.domain.number_of_elements, self._rates[idx])

    def _slice_duration(self, idx):
        return self._durations[idx]

    def __call__(self):
        t = self.domain.get_time()
        if self._times[0] <= t < self._times[-1]:
            self._set_slice(int(np.searchsorted(self._times, t, side='right')) - 1)
        else:
            self._set_slice(None)

        super().__call__()
        self._count_loss()


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
    losses : Rainfall_loss, optional
        A rainfall loss model, such as :class:`~anuga.rain.SCS_curve_number`
        or :class:`~anuga.rain.Initial_continuing_loss`.  The operator then
        applies the excess rain only.  ``factor`` scales the excess, and
        ``default_rate`` is applied with no losses.

    Attributes
    ----------
    time_step : float
        Time (s) between changes of rate, ``pattern.Tstep`` in seconds. Use it, or
        a whole fraction of it, as the ``yieldstep`` of
        :meth:`~anuga.Domain.evolve` to deliver the exact depth.
    cumulative_loss : float
        Volume (m^3) of rain lost so far when ``losses`` is given.
        ``cumulative_influx`` is the water added, which is the excess.

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
                 ifd_duration=None,
                 losses=None):

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
        self._set_losses(losses)

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
        self._count_loss()
