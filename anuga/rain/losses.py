"""
Rainfall loss models for rain-on-mesh.

A loss model turns the rain falling on each cell into the part that runs
off (the *excess*). It sees the rain alone, never the flow: water running on
to a cell from upslope does not infiltrate.

Classes
-------
SCS_curve_number
    The NRCS (SCS) curve number method, as used in the US.
Initial_continuing_loss
    Initial loss followed by a continuing loss rate, as used with Australian
    Rainfall and Runoff (ARR).
Initial_proportional_loss
    Initial loss followed by the loss of a fixed fraction of the rain: ARR's
    IL/PL, and a fixed runoff coefficient or percentage runoff.

Pass one to a rain operator::

    losses = anuga.rain.SCS_curve_number(curve_number=75.0)
    anuga.Raster_rate_operator(domain, rain_data, losses=losses)

The operator hands the model each rain time slice in turn, and applies the
excess it gets back. Rain is taken as uniform within a slice.
"""

import numpy as np

MM = 1.0e-3            # metres in a millimetre
MM_PER_HOUR = MM / 3600.0


def per_cell(value, domain, name):
    """``value`` as one finite float per cell of ``domain``.

    ``value`` may be a number, an array with one entry per cell, a Quantity
    (its centroid values are used) or a function ``f(x, y)`` of the centroid
    coordinates, which are relative to the domain's origin as they are for
    ``domain.set_quantity``.
    """
    N = domain.number_of_elements
    if hasattr(value, 'centroid_values'):
        out = np.array(value.centroid_values, dtype=float)
    elif callable(value):
        xy = domain.get_centroid_coordinates()
        out = np.asarray(value(xy[:, 0], xy[:, 1]), dtype=float)
        if out.ndim == 0:
            out = np.full(N, float(out))
    elif np.ndim(value) == 0:
        out = np.full(N, float(value))
    else:
        out = np.array(value, dtype=float)

    if out.shape != (N,):
        raise ValueError(f'{name} has shape {out.shape}; it must be a number '
                         f'or have one value per cell ({N},)')
    if not np.isfinite(out).all():
        raise ValueError(f'{name} has NaN or infinite values')
    return out


class Rainfall_loss:
    """Base class of the rainfall loss models.

    A subclass sets its per-cell parameters in :meth:`_setup` and gives the
    excess of one slice of rain in :meth:`_excess`.

    Attributes
    ----------
    cumulative_rain : ndarray
        Rain (m) given to the model so far, per cell.
    cumulative_excess : ndarray
        Excess rain (m) returned so far, per cell.
    """

    def __init__(self):
        self.domain = None
        self.cumulative_rain = None
        self.cumulative_excess = None

    def bind(self, domain):
        """Attach the model to ``domain`` and start with no rain fallen.

        Called by the rain operator. A model keeps the state of one storm,
        so it can serve one operator only.
        """
        if self.domain is not None:
            raise RuntimeError(
                f'this {type(self).__name__} is already used by a rain '
                f'operator; make a separate one for each operator')
        self.domain = domain
        N = domain.number_of_elements
        self.cumulative_rain = np.zeros(N)
        self.cumulative_excess = np.zeros(N)
        self._setup(domain)

    @property
    def cumulative_loss(self):
        """Rain (m) lost so far, per cell."""
        return self.cumulative_rain - self.cumulative_excess

    def excess(self, depth, duration):
        """Excess rain (m) per cell from a slice of rain.

        Parameters
        ----------
        depth : ndarray
            Rain (m) falling on each cell in the slice.
        duration : float
            Length of the slice (s).

        The slices of a storm must be given once each and in order: the call
        moves the model's state on to the end of the slice.
        """
        if self.domain is None:
            raise RuntimeError('bind the loss model to a domain first')
        depth = np.maximum(np.asarray(depth, dtype=float), 0.0)
        out = self._excess(depth, float(duration))
        self.cumulative_rain += depth
        self.cumulative_excess += out
        return out

    def _setup(self, domain):
        raise NotImplementedError

    def _excess(self, depth, duration):
        raise NotImplementedError


class SCS_curve_number(Rainfall_loss):
    """NRCS (SCS) curve number losses.

    The runoff depth :math:`Q` from a cumulative rain depth :math:`P` is

    .. math::

        Q = \\frac{(P - I_a)^2}{P - I_a + S} \\quad (P > I_a),
        \\qquad Q = 0 \\quad (P \\le I_a)

    with the potential retention :math:`S = 25400/CN - 254` (mm) and the
    initial abstraction :math:`I_a = \\lambda S`. The excess of a slice of
    rain is the change in :math:`Q` across it.

    Parameters
    ----------
    curve_number : number, array, Quantity or f(x, y)
        The curve number CN, above 0 and at most 100, from the NRCS tables
        for the land use, hydrologic soil group and antecedent moisture.
        CN = 100 loses nothing.
    initial_abstraction_ratio : float, optional
        :math:`\\lambda`. The traditional value 0.2 is the default; 0.05 is
        used with curve numbers that were derived for it.
    """

    def __init__(self, curve_number, initial_abstraction_ratio=0.2):
        super().__init__()
        if initial_abstraction_ratio < 0.0:
            raise ValueError('initial_abstraction_ratio must not be negative')
        self.curve_number = curve_number
        self.initial_abstraction_ratio = float(initial_abstraction_ratio)

    def _setup(self, domain):
        cn = per_cell(self.curve_number, domain, 'curve_number')
        if (cn <= 0.0).any() or (cn > 100.0).any():
            raise ValueError('curve_number must be above 0 and at most 100')
        self.retention = (25400.0/cn - 254.0) * MM                 # S (m)
        self.initial_abstraction = self.initial_abstraction_ratio * self.retention

    def runoff(self, rain):
        """Cumulative runoff depth (m) for a cumulative rain depth (m)."""
        over = np.maximum(rain - self.initial_abstraction, 0.0)
        denominator = over + self.retention
        return np.divide(over*over, denominator,
                         out=np.zeros_like(over), where=denominator > 0.0)

    def _excess(self, depth, duration):
        before = self.cumulative_rain
        return self.runoff(before + depth) - self.runoff(before)


class Initial_continuing_loss(Rainfall_loss):
    """Initial loss followed by a continuing loss rate (ARR IL/CL).

    Rain first fills the initial loss, and none of it runs off. After that
    the continuing loss is taken from the rain at a steady rate, and the
    excess is what is left, never below zero.

    Within a slice the rain is uniform, so a slice that fills the initial
    loss part of the way through loses the continuing rate over the rest of
    it only.

    Parameters
    ----------
    initial : number, array, Quantity or f(x, y)
        Initial loss (mm).
    continuing : number, array, Quantity or f(x, y)
        Continuing loss (mm/h).

    Attributes
    ----------
    initial_remaining : ndarray
        Initial loss (m) not yet filled, per cell.
    """

    def __init__(self, initial, continuing):
        super().__init__()
        self.initial = initial
        self.continuing = continuing

    def _setup(self, domain):
        initial = per_cell(self.initial, domain, 'initial')
        continuing = per_cell(self.continuing, domain, 'continuing')
        if (initial < 0.0).any() or (continuing < 0.0).any():
            raise ValueError('initial and continuing losses must not be negative')
        self.initial_remaining = initial * MM            # m
        self.continuing_rate = continuing * MM_PER_HOUR  # m/s

    def _excess(self, depth, duration):
        to_initial = np.minimum(depth, self.initial_remaining)
        self.initial_remaining = self.initial_remaining - to_initial
        after = depth - to_initial
        # The share of the slice left once the initial loss is full
        share = np.divide(after, depth, out=np.zeros_like(depth), where=depth > 0.0)
        return np.maximum(after - self.continuing_rate*duration*share, 0.0)


class Initial_proportional_loss(Rainfall_loss):
    """Initial loss followed by the loss of a fixed fraction of the rain.

    Rain first fills the initial loss, and none of it runs off. After that
    the fraction ``proportional`` of the rain is lost and the rest runs off.
    This is ARR's initial loss / proportional loss model. With no initial
    loss it is a runoff coefficient (or fixed percentage runoff) of
    ``1 - proportional``.

    Parameters
    ----------
    proportional : number, array, Quantity or f(x, y)
        Fraction of the rain that is lost, from 0 to 1.
    initial : number, array, Quantity or f(x, y), optional
        Initial loss (mm). Default 0.

    Attributes
    ----------
    initial_remaining : ndarray
        Initial loss (m) not yet filled, per cell.
    """

    def __init__(self, proportional, initial=0.0):
        super().__init__()
        self.proportional = proportional
        self.initial = initial

    def _setup(self, domain):
        fraction = per_cell(self.proportional, domain, 'proportional')
        initial = per_cell(self.initial, domain, 'initial')
        if (fraction < 0.0).any() or (fraction > 1.0).any():
            raise ValueError('proportional must be from 0 to 1')
        if (initial < 0.0).any():
            raise ValueError('initial loss must not be negative')
        self.runoff_fraction = 1.0 - fraction
        self.initial_remaining = initial * MM            # m

    def _excess(self, depth, duration):
        to_initial = np.minimum(depth, self.initial_remaining)
        self.initial_remaining = self.initial_remaining - to_initial
        return (depth - to_initial)*self.runoff_fraction
