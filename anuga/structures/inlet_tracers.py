"""Tracer mass carried by the inlet and structure operators.

Inlet operators: see InletTracers below. Structures (culverts, weirs,
bridges, internal boundaries): see StructureTracers. Rate operators: see
RateTracers at the end.


An inlet operator adds or removes water over a region of the mesh by changing
the stage there. Without this module it never touched the tracers, so water
added at an inlet came in clean (diluting whatever was there) and water taken
out at an outlet left its tracer behind, which then concentrated without bound
in the outlet cells as they drained.

The rules, applied after the water has been moved:

* **Inflow** (volume >= 0). The added water carries the inflow concentration
  given for each tracer (0 when none is given, which is the old behaviour). The
  inlet fills its lowest cells first to a level stage, so each cell gains water
  ``dh_i >= 0`` and gains tracer ``c_in * dh_i``.

* **Extraction** (volume < 0). The operator treats the inlet as one level pool
  (it sets a uniform depth), so the water leaves at the pool's mean
  concentration ``C = sum(m_i A_i) / sum(h_i A_i)`` and the cells left behind
  hold ``C * h_new``. The tracer removed is exactly ``C`` times the water
  removed, and draining the inlet dry removes all of it.

Both rules conserve tracer: the change in the inlet's tracer mass is the
tracer carried by the water that crossed it, no more and no less.
"""

import numpy as num


class InletTracers:
    """Tracer bookkeeping for one inlet operator.

    Parameters
    ----------
    domain : Domain
    concentrations : dict, optional
        Tracer name -> inflow concentration, a float or a callable of time.
        Tracers not named come in at zero concentration.
    label : str
        The operator's label, for messages.
    """

    def __init__(self, domain, concentrations=None, label=''):
        self.domain = domain
        self.concentrations = dict(concentrations or {})
        self.label = label
        # Tracer volume gained / lost by this operator's cells, per tracer name,
        # summed step by step by sign. In parallel a rank can gain tracer from
        # an extraction (water moves between ranks inside the pool), so across
        # ranks it is the net, total_in - total_out, that is the transfer.
        self.total_in = {}
        self.total_out = {}
        self._checked_ns = -1

    def active(self):
        return getattr(self.domain, 'number_of_tracers', 0) > 0

    def _check_names(self):
        ns = self.domain.number_of_tracers
        if ns == self._checked_ns:
            return
        names = self.domain.get_tracer_names()
        unknown = sorted(set(self.concentrations) - set(names))
        if unknown:
            raise ValueError(
                'inlet operator %r: tracer_concentrations names %s, which the '
                'domain does not have; its tracers are %s'
                % (self.label, unknown, names))
        self._checked_ns = ns

    def inflow_concentrations(self, t, timestep):
        """Inflow concentration of every tracer over the step, (ns,).

        A callable is averaged over the step's two ends, as the discharge is.
        """
        self._check_names()
        names = self.domain.get_tracer_names()
        c = num.zeros(len(names))
        for i, name in enumerate(names):
            value = self.concentrations.get(name, 0.0)
            if callable(value):
                value = 0.5 * (float(num.ravel(value(t))[0])
                               + float(num.ravel(value(t + timestep))[0]))
            c[i] = float(value)
        return c

    def capture(self, indices):
        """State before the water moves: depths (n,) and tracer mass (ns, n)."""
        d = self.domain
        h = num.maximum(d.quantities['stage'].centroid_values[indices]
                        - d.quantities['elevation'].centroid_values[indices], 0.0)
        m = d.tracer_conserved_values[:, indices].copy()
        return h, m

    def apply(self, indices, areas, h_old, m_old, volume, t, timestep,
              global_sum=None):
        """Update the tracer mass in the inlet cells after the water has moved.

        `global_sum(x)` sums an array over every process holding part of the
        inlet; given for the parallel inlet, where the pool spans ranks.
        Returns the change in tracer volume, (ns,).
        """
        d = self.domain
        self._check_names()
        h_new = num.maximum(d.quantities['stage'].centroid_values[indices]
                            - d.quantities['elevation'].centroid_values[indices], 0.0)
        if volume >= 0.0:
            c_in = self.inflow_concentrations(t, timestep)
            dh = num.maximum(h_new - h_old, 0.0)
            m_new = m_old + c_in[:, None] * dh[None, :]
        else:
            M = (m_old * areas[None, :]).sum(axis=1)
            V = float((h_old * areas).sum())
            if global_sum is not None:
                M = num.asarray(global_sum(M), dtype=float)
                V = float(global_sum(num.array([V]))[0])
            C = M / V if V > 0.0 else num.zeros_like(M)
            m_new = C[:, None] * h_new[None, :]

        d.tracer_conserved_values[:, indices] = m_new
        hmin = d.minimum_allowed_height
        with num.errstate(divide='ignore', invalid='ignore'):
            inv_h = num.where(h_new > hmin, 1.0 / h_new, 0.0)
        d.tracer_centroid_values[:, indices] = m_new * inv_h[None, :]

        dmass = ((m_new - m_old) * areas[None, :]).sum(axis=1)
        self.record(dmass)
        return dmass

    def record(self, dmass):
        for name, dm in zip(self.domain.get_tracer_names(), dmass):
            if dm >= 0.0:
                self.total_in[name] = self.total_in.get(name, 0.0) + float(dm)
            else:
                self.total_out[name] = self.total_out.get(name, 0.0) - float(dm)


class StructureTracers:
    """Tracer carried by a structure's water transfer (culvert, weir, bridge,
    internal boundary).

    A structure draws its inflow region down and fills its outflow region
    (Inlet.set_average_depth levels the surface: on a transfer a drawn-down
    cell only loses water and a filled cell only gains it). The rule, per
    cell, after the water has moved:

    * a cell that LOST water keeps its concentration, so it gives up
      ``c_i * dh_i`` of tracer;
    * a cell that GAINED water receives the total tracer given up, divided by
      the total water gained, in proportion to its own gain.

    The tracer lost and gained balance exactly. Being per cell, an idle
    structure (no transfer) leaves the tracer untouched: nothing is mixed
    across an inlet region that no water crossed.
    """

    def __init__(self, domain):
        self.domain = domain
        self.total_moved = {}   # tracer volume moved through, per tracer name

    def active(self):
        return getattr(self.domain, 'number_of_tracers', 0) > 0

    def capture(self, inlets):
        """State of each inlet region before the water moves (None if absent
        on this rank)."""
        d = self.domain
        stage = d.quantities['stage'].centroid_values
        bed = d.quantities['elevation'].centroid_values
        state = []
        for inlet in inlets:
            if inlet is None or len(inlet.triangle_indices) == 0:
                state.append(None)
                continue
            idx = num.asarray(inlet.triangle_indices, dtype=int)
            h = num.maximum(stage[idx] - bed[idx], 0.0)
            state.append((idx, d.areas[idx].copy(), h,
                          d.tracer_conserved_values[:, idx].copy()))
        return state

    def apply(self, state, global_sum=None):
        """Move the tracer with the water. `global_sum(x)` sums an array over
        every rank of the structure (parallel); returns the tracer moved."""
        d = self.domain
        ns = d.number_of_tracers
        stage = d.quantities['stage'].centroid_values
        bed = d.quantities['elevation'].centroid_values
        T = num.zeros(ns)
        G = 0.0
        work = []
        for s in state:
            if s is None:
                work.append(None)
                continue
            idx, areas, h_old, m_old = s
            h_new = num.maximum(stage[idx] - bed[idx], 0.0)
            dh = h_new - h_old
            lose = dh < 0.0
            gain = dh > 0.0
            m_new = m_old.copy()
            with num.errstate(divide='ignore', invalid='ignore'):
                keep = num.where(h_old > 0.0, h_new / h_old, 0.0)
            m_new[:, lose] = m_old[:, lose] * keep[None, lose]
            T += ((m_old - m_new) * areas[None, :]).sum(axis=1)
            G += float((dh[gain] * areas[gain]).sum())
            work.append((idx, h_new, dh, gain, m_new))
        if global_sum is not None:
            T = num.asarray(global_sum(T), dtype=float)
            G = float(global_sum(num.array([G]))[0])
        c_t = T / G if G > 0.0 else num.zeros(ns)

        hmin = d.minimum_allowed_height
        for w in work:
            if w is None:
                continue
            idx, h_new, dh, gain, m_new = w
            m_new[:, gain] += c_t[:, None] * dh[None, gain]
            d.tracer_conserved_values[:, idx] = m_new
            with num.errstate(divide='ignore', invalid='ignore'):
                inv_h = num.where(h_new > hmin, 1.0 / h_new, 0.0)
            d.tracer_centroid_values[:, idx] = m_new * inv_h[None, :]

        moved = c_t * G if G > 0.0 else num.zeros(ns)
        for name, v in zip(d.get_tracer_names(), moved):
            self.total_moved[name] = self.total_moved.get(name, 0.0) + float(v)
        return moved


class RateTracers:
    """Tracer carried by a Rate_operator's water, cell by cell.

    * Water ADDED (rain, a distributed inflow) carries the concentration given
      in ``tracer_concentrations`` (0 when none is given: clean rain, the old
      behaviour).
    * Water REMOVED (a negative rate) either takes each tracer with it at the
      cell's concentration (``'carry'``: pumping, drains, abstraction, or a
      dissolved tracer lost to infiltration) or leaves it behind (``'retain'``:
      salt under evaporation, or sediment the bed filters out). ``'retain'``
      is the old behaviour and the default.

    Exact per cell. Nothing is done when no concentration is given and every
    tracer retains, so an ordinary rainfall operator costs nothing extra.
    """

    MODES = ('retain', 'carry')

    def __init__(self, domain, concentrations=None, extraction='retain', label=''):
        self.domain = domain
        self.concentrations = dict(concentrations or {})
        if isinstance(extraction, str):
            if extraction not in self.MODES:
                raise ValueError('tracer_extraction must be %r or a dict of them, got %r'
                                 % (self.MODES, extraction))
            self.default_mode = extraction
            self.modes = {}
        else:
            self.default_mode = 'retain'
            self.modes = dict(extraction or {})
            bad = {k: v for k, v in self.modes.items() if v not in self.MODES}
            if bad:
                raise ValueError('tracer_extraction values must be %r, got %r'
                                 % (self.MODES, bad))
        self.label = label
        self._checked_ns = -1

    def _check_names(self):
        ns = self.domain.number_of_tracers
        if ns == self._checked_ns:
            return
        names = self.domain.get_tracer_names()
        unknown = sorted((set(self.concentrations) | set(self.modes)) - set(names))
        if unknown:
            raise ValueError('rate operator %r: tracer names %s, which the domain '
                             'does not have; its tracers are %s'
                             % (self.label, unknown, names))
        self._checked_ns = ns

    def needed(self):
        """True when this operator has anything to do with the tracers."""
        if getattr(self.domain, 'number_of_tracers', 0) == 0:
            return False
        if any(callable(v) or v != 0.0 for v in self.concentrations.values()):
            return True
        return self.default_mode == 'carry' or 'carry' in self.modes.values()

    def settings(self, t, timestep):
        """(c_in, carry) arrays over the domain's tracers, c_in averaged over
        the step as a callable Q is."""
        self._check_names()
        names = self.domain.get_tracer_names()
        c_in = num.zeros(len(names))
        carry = num.zeros(len(names), dtype=num.intc)
        for i, name in enumerate(names):
            v = self.concentrations.get(name, 0.0)
            if callable(v):
                v = 0.5 * (float(num.ravel(v(t))[0]) + float(num.ravel(v(t + timestep))[0]))
            c_in[i] = float(v)
            carry[i] = 1 if self.modes.get(name, self.default_mode) == 'carry' else 0
        return c_in, carry

    def capture(self, indices):
        d = self.domain
        stage = d.quantities['stage'].centroid_values
        bed = d.quantities['elevation'].centroid_values
        if indices is None:
            return num.maximum(stage - bed, 0.0)
        return num.maximum(stage[indices] - bed[indices], 0.0)

    def apply(self, indices, h_old, t, timestep):
        d = self.domain
        c_in, carry = self.settings(t, timestep)
        stage = d.quantities['stage'].centroid_values
        bed = d.quantities['elevation'].centroid_values
        sel = slice(None) if indices is None else indices
        h_new = num.maximum(stage[sel] - bed[sel], 0.0)
        dh = h_new - h_old
        m = d.tracer_conserved_values[:, sel]
        gain = dh > 0.0
        if gain.any() and num.any(c_in != 0.0):
            m[:, gain] += c_in[:, None] * dh[None, gain]
        lose = dh < 0.0
        if lose.any() and carry.any():
            with num.errstate(divide='ignore', invalid='ignore'):
                keep = num.where(h_old > 0.0, h_new / h_old, 0.0)
            rows = num.nonzero(carry)[0]
            m[num.ix_(rows, num.nonzero(lose)[0])] *= keep[None, lose]
        d.tracer_conserved_values[:, sel] = m
        hmin = d.minimum_allowed_height
        with num.errstate(divide='ignore', invalid='ignore'):
            inv_h = num.where(h_new > hmin, 1.0 / h_new, 0.0)
        d.tracer_centroid_values[:, sel] = m * inv_h[None, :]
