.. _mathematical_background:

.. currentmodule:: anuga

Mathematical Background
=======================

This page outlines the mathematics underpinning ANUGA.  It is a
web-searchable version of Chapter 16 of the
`ANUGA User Manual <https://github.com/anuga-community/anuga_user_manual>`_.
ANUGA has been applied to hydrodynamic modelling of coastal inundation
[Nielsen2005]_ and many other shallow water problems.

.. contents:: Contents
   :local:
   :depth: 2


Shallow Water Wave Equations
-----------------------------

The shallow water wave equations are a system of differential conservation
equations describing the flow of a thin layer of fluid over terrain:

.. math::

   \frac{\partial \mathbf{U}}{\partial t}
   + \frac{\partial \mathbf{E}}{\partial x}
   + \frac{\partial \mathbf{G}}{\partial y}
   = \mathbf{S}

The vector of **conserved quantities** is

.. math::

   \mathbf{U} = \begin{bmatrix} h \\ uh \\ vh \end{bmatrix}

where :math:`h` is water depth, :math:`uh` is momentum in the
:math:`x`-direction, and :math:`vh` is momentum in the :math:`y`-direction.
Other quantities entering the system are bed elevation :math:`z` and stage
(absolute water level) :math:`w`, related by :math:`w = z + h` at all times.

The **flux vectors** in the :math:`x` and :math:`y` directions are

.. math::

   \mathbf{E} = \begin{bmatrix} uh \\ u^2 h + \tfrac{1}{2}gh^2 \\ uvh \end{bmatrix},
   \qquad
   \mathbf{G} = \begin{bmatrix} vh \\ uvh \\ v^2 h + \tfrac{1}{2}gh^2 \end{bmatrix}

and the **source term** (gravity and friction) is

.. math::

   \mathbf{S} = \begin{bmatrix} 0 \\ -gh(z_x + S_{fx}) \\ -gh(z_y + S_{fy}) \end{bmatrix}

where :math:`S_f` is the bed friction modelled by Manning's resistance law:

.. math::

   S_{fx} = \frac{u \eta^2 \sqrt{u^2 + v^2}}{h^{4/3}},
   \qquad
   S_{fy} = \frac{v \eta^2 \sqrt{u^2 + v^2}}{h^{4/3}}

with :math:`\eta` the Manning resistance coefficient.  The model does not
include kinematic viscosity or dispersion (though :class:`Kinematic_viscosity_operator`
can be added as an explicit operator).


.. _tracer_transport_equation:

Passive Tracer Transport
------------------------

A **passive tracer** is a depth-averaged concentration :math:`c` carried by the
flow: a salinity, a dye, a pollutant. It is passive in the sense that it is
transported by the water but does not act back on it, so the equations above are
unchanged by its presence.

What is conserved is not the concentration but the **mass per unit area**

.. math::

   m = h\,c

which obeys a conservation law with the same structure as the mass equation, and
with the same water flux:

.. math::

   \frac{\partial m}{\partial t}
   + \frac{\partial (u m)}{\partial x}
   + \frac{\partial (v m)}{\partial y}
   = 0

For :math:`N_s` tracers this is one such equation per tracer, so the state vector
carried by the solver is

.. math::

   \mathbf{U} = \begin{bmatrix} h & uh & vh & m_1 & \dots & m_{N_s}\end{bmatrix}^{T}

Two consequences follow from conserving :math:`m` rather than :math:`c`, and both
are visible in the interface.

*Concentration is derived, not stored.* :math:`c = m/h` is recomputed from the
conserved variable each substep, exactly as height is derived from stage. Where
the cell is dry the concentration is taken to be zero rather than the last wet
value, since :math:`m/h` is meaningless there.

*The tracer flux uses the same edge flux as the water.* Writing
:math:`F` for the water mass flux already computed across an edge, the tracer
flux across that edge is :math:`F\,c^{\ast}`, where :math:`c^{\ast}` is the
concentration on the upwind side:

.. math::

   c^{\ast} =
   \begin{cases}
     c_{\text{interior}}, & F \le 0 \quad \text{(outflow)} \\
     c_{\text{boundary}}, & F > 0 \quad \text{(inflow)}
   \end{cases}

Using the one water flux for both cells sharing an edge is what makes tracer mass
conservation structural: whatever leaves one cell enters the other, independently
of the two cells' sizes. It is also why a boundary concentration is only
meaningful on inflow -- on outflow the interior value is carried out, and
prescribing anything there would over-determine the advection.

A source term may be added on the right-hand side; suspended sediment is exactly
that case, with erosion and deposition supplying it. See
:ref:`sediment_governing_equations`.


Finite Volume Method
---------------------

ANUGA uses a finite-volume method to solve the shallow water wave equations
[ZR1999]_.  The study area is represented by a mesh of triangular cells;
the conserved quantities :math:`h`, :math:`uh`, :math:`vh` are stored at the
centroid of each cell.  The size of the triangles may vary to provide greater
resolution in regions of particular interest.

**Cell update equation**

Integrating the conservation equations over triangular cell :math:`T_i` and
applying the divergence theorem gives the rate equation for the cell-averaged
conserved quantities:

.. math::

   \frac{d\mathbf{U}_i}{dt}
   + \frac{1}{A_i} \sum_{j \in \mathcal{N}_i} \mathbf{H}_{ij}\, l_{ij}
   = \mathbf{S}_i

where

- :math:`\mathbf{U}_i` — conserved quantities averaged over cell :math:`i`
- :math:`A_i` — area of cell :math:`i`
- :math:`\mathcal{N}_i` — set of indices :math:`j` of cells neighbouring cell :math:`i`
- :math:`l_{ij}` — length of edge :math:`e_{ij}` between cells :math:`i` and :math:`j`
- :math:`\mathbf{H}_{ij}\,l_{ij}` — outward normal flux across edge :math:`e_{ij}`
- :math:`\mathbf{S}_i` — average source term over cell :math:`i`

**Reconstruction**

A second-order piecewise-linear reconstruction of the conserved quantities
is built within each cell from the cell-average values of that cell and its
neighbours.  The slope is limited (see `Slope Limiting`_ below) to suppress
spurious oscillations.  The reconstruction is allowed to be discontinuous
across edges.  The values on either side of edge :math:`e_{ij}` are denoted

.. math::

   \mathbf{U}^i_{ij} = \lim_{\mathbf{x} \to m_{ij}} \mathbf{U}_i(\mathbf{x}),
   \qquad
   \mathbf{U}^j_{ij} = \lim_{\mathbf{x} \to m_{ij}} \mathbf{U}_j(\mathbf{x})

where :math:`m_{ij}` is the midpoint of edge :math:`e_{ij}`.

**Numerical flux — central-upwind scheme**

The normal flux :math:`\mathbf{H}(\mathbf{U}) = \mathbf{E}(\mathbf{U})\,n_1 + \mathbf{G}(\mathbf{U})\,n_2`
is approximated using the central-upwind scheme of Kurganov, Noelle & Petrova
[KurNP2001]_:

.. math::

   \mathbf{H}_{ij}
   = \frac{a^+_{ij}\,\mathbf{H}(\mathbf{U}^i_{ij})
           - a^-_{ij}\,\mathbf{H}(\mathbf{U}^j_{ij})}
          {a^+_{ij} - a^-_{ij}}
   + \frac{a^+_{ij}\,a^-_{ij}}{a^+_{ij} - a^-_{ij}}
     \bigl[\mathbf{U}^j_{ij} - \mathbf{U}^i_{ij}\bigr]

with one-sided wave speed bounds

.. math::

   a^+_{ij} &= \max\!\left\{
       u^i_{ij} + \sqrt{g h^i_{ij}},\;
       u^j_{ij} + \sqrt{g h^j_{ij}},\;
       0
   \right\} \\[4pt]
   a^-_{ij} &= \min\!\left\{
       u^i_{ij} - \sqrt{g h^i_{ij}},\;
       u^j_{ij} - \sqrt{g h^j_{ij}},\;
       0
   \right\}

**Time stepping**

An explicit Euler scheme with adaptive timestep is used:

.. math::

   \mathbf{U}^{n+1}_i
   = \mathbf{U}^n_i
   - \frac{\Delta t}{A_i} \sum_{j \in \mathcal{N}_i} \mathbf{H}_{ij}\,l_{ij}
   + \Delta t\,\mathbf{S}_i

Stability requires the timestep to satisfy the Courant–Friedrichs–Lewy (CFL)
condition:

.. math::

   \Delta t \le \min_{i}\; \min_{j \in \mathcal{N}_i}
   \left( \frac{r_i}{a_{ij}},\; \frac{r_j}{a_{ij}} \right)

where :math:`a_{ij} = \max\{a^+_{ij},\, -a^-_{ij}\}` and :math:`r_i` is the
inradius (inscribed circle radius) of cell :math:`i`.


Flux Limiting
--------------

Near wet/dry boundaries the water depth :math:`h` becomes very small, making
the velocity recovery :math:`u = uh/h` numerically unreliable and potentially
producing unphysical speeds.

ANUGA replaces the raw velocity calculation with the limited approximation

.. math::

   \hat{u} = \frac{uh \cdot h}{h^2 + h_0},
   \qquad
   \hat{v} = \frac{vh \cdot h}{h^2 + h_0}

where :math:`h_0` is a regularisation parameter.  Taking limits:

.. math::

   \lim_{h \to 0} \hat{u} = 0,
   \qquad
   \lim_{h \to \infty} \hat{u} = \frac{uh}{h} = u

so the limiter smoothly recovers zero velocity at dryness and the true
velocity in deep water.

ANUGA exposes a global minimum-depth parameter :math:`H_0` (typically
:math:`10^{-3}` m, set in ``anuga/config.py`` as ``minimum_allowed_height``).
Setting

.. math::

   h_0 = H_0^2

provides a good balance between accuracy and stability.  At depth
:math:`h = N H_0` the predicted speed is scaled by

.. math::

   \frac{1}{1 + 1/N^2}

which converges quadratically to the true value with :math:`N`.  The limiter
is applied only for depths less than :math:`10\,H_0` (roughly 1 cm for the
default :math:`H_0`), where it affects the computed velocity by less than 1 %.


Slope Limiting
---------------

A multidimensional slope-limiting technique (MinMod limiter) achieves
second-order spatial accuracy and prevents spurious oscillations.  Near the
bed, the limiter must additionally ensure that no negative depths occur.

Let :math:`w`, :math:`z`, :math:`h` be stage, bed elevation, and depth at the
centroid, and let :math:`w_i`, :math:`z_i`, :math:`h_i` be the corresponding
vertex values.  Define the minimum vertex depth as
:math:`h_{\min} = \min_i h_i`.

Two reconstructed stages are computed:

- :math:`\tilde{w}_i` — stage from a gradient limiter applied to stage only
  (suitable for deep water where the bed slope can be ignored)
- :math:`\bar{w}_i` — stage from a gradient limiter applied to depth,
  respecting the bed slope (suitable for shallow water near wet/dry fronts)

The **balanced stage** is the linear combination

.. math::

   w_i = \alpha\,\tilde{w}_i + (1-\alpha)\,\bar{w}_i, \qquad \alpha \in [0,1]

In deep water (:math:`h_{\min} \ge \varepsilon`) we set :math:`\alpha = 1`
(stage-limited reconstruction only).

In shallow water (:math:`h_{\min} < \varepsilon`) we choose the largest
:math:`\alpha` that still guarantees non-negative vertex depths:

.. math::

   \alpha = \min_i \frac{\bar{h}_i - \varepsilon}{\bar{h}_i - \tilde{h}_i}

where :math:`\bar{h}_i = \bar{w}_i - z_i` and
:math:`\tilde{h}_i = \tilde{w}_i - z_i`.  Should :math:`\alpha` fall outside
:math:`[0, 1]` it is clamped to that interval.


.. seealso::

   `ANUGA User Manual — Chapter 16: Mathematical Background
   <https://github.com/anuga-community/anuga_user_manual>`_
   is the original source for this page and includes the parallel
   implementation theory, scalability results, and the full TikZ mesh diagrams.

References
-----------

.. [ZR1999] Zoppou, C. & Roberts, S. (1999). *Catastrophic collapse of water
   supply reservoirs in urban areas*. Journal of Hydraulic Engineering, 125(7),
   686–695.

.. [KurNP2001] Kurganov, A., Noelle, S. & Petrova, G. (2001). *Semidiscrete
   central-upwind schemes for hyperbolic conservation laws and Hamilton–Jacobi
   equations*. SIAM Journal on Scientific Computing, 23(3), 707–740.

.. [Nielsen2005] Nielsen, O., Roberts, S., Gray, D., McPherson, A. & Hitchman,
   A. (2005). *Hydrodynamic modelling of coastal inundation*. MODSIM 2005
   International Congress on Modelling and Simulation, 518–523.
