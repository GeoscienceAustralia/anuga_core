Vegetation drag
===============

Emergent or submerged vegetation slows the flow far more than a bare bed of
the same material, and on a delta or a floodplain the stem field is what sets
where the water goes. ANUGA can carry that drag alongside Manning friction:

.. code-block:: python

   domain.set_vegetation_drag(density=m, diameter=D, height=hv)

``density`` is the stem density :math:`m` in stems per square metre,
``diameter`` the stem diameter :math:`D` in metres and ``height`` the stem
height :math:`h_v` in metres. Each is given as :meth:`~anuga.Domain.set_quantity` takes a
value: a constant, an array with one value per cell centroid, or a function
of ``x`` and ``y``. They become the quantities ``veg_density``,
``veg_diameter`` and ``veg_height``. A cell with zero density or zero height
has no vegetation and is left to Manning friction alone, so the usual
arrangement is a Manning ``n`` on the open water and channels, and the
stem field, with ``n`` = 0, on the vegetated classes.

The formulation is that of Baptist et al. (2007): the vegetated Chezy
coefficient

.. math::

   C_v = \left( C_b^{-2} + \frac{C_D\, m\, D\, \min(h, h_v)}{2 g} \right)^{-1/2}
         + \frac{\sqrt{g}}{\kappa} \ln \frac{\max(h, h_v)}{h_v}

gives the friction slope :math:`g\, q\, |q| / (C_v^2 h^2)`, applied
semi-implicitly like Manning friction and on top of it. Stems taller than
the flow are emergent and the first term is the whole story; once the flow
overtops them the logarithmic term adds the faster layer above the canopy.
The stem drag coefficient :math:`C_D` (1.68) and the bed Chezy coefficient
:math:`C_b` (65) inside the stems can be changed with the ``Cd`` and
``bed_chezy`` arguments. The values above are those of the Delta-X Wax Lake
Delta model of Wright and Passalacqua, where the wetland classes carry
120 to 200 stems per square metre of 1 to 1.5 cm diameter.

The drag runs in the friction kernel of both compute modes, so it offloads
with the rest of the step on a GPU build. Call :meth:`~anuga.Domain.set_vegetation_drag`
before the first ``evolve``: the three fields are copied to the device once.
On a distributed run call it after ``distribute()``, on every rank, as for
any :meth:`~anuga.Domain.set_quantity`. :meth:`set_vegetation_drag(formulation='off') <anuga.Domain.set_vegetation_drag>` removes it.

Vegetation and sediment
-----------------------

With sediment transport on, the vegetation also decides what bed shear the
sediment kernel sees in a vegetated cell -- for entrainment, the Rouse
near-bed profile and bedload alike. The stems take most of the drag, so the
bed feels far less than the total resistance. Three choices, through
``sediment_shear``:

.. list-table::
   :header-rows: 1
   :widths: 14 46 40

   * - ``sediment_shear``
     - bed shear :math:`\tau_b/\rho`
     - when
   * - ``'bed'`` (default)
     - :math:`g\,u_v^2 / C_b^2`, with the canopy velocity
       :math:`u_v = U\,C_{v,r}/C_v` (:math:`u_v = U` for emergent stems)
     - the bed's share of the resistance, as Baptist et al. (2007) split it
   * - ``'total'``
     - :math:`g\,U^2 / C_v^2`, the whole vegetated resistance
     - to reproduce models that apply the total friction to the bed, such as
       the Delta-X Wax Lake Delta sediment model
   * - ``'ignore'``
     - the cell's own sediment friction closure
     - with :math:`n = 0` on the vegetated classes this is no shear at all

:math:`C_{v,r} = (C_b^{-2} + C_D m D \min(h, h_v)/(2g))^{-1/2}` is the
canopy part of the vegetated Chezy coefficient :math:`C_v`. Cells without
stems keep their own closure whatever the setting.

.. code-block:: python

   domain.set_vegetation_drag(density=m, diameter=D, height=hv,
                              sediment_shear='total')

.. seealso::

   :ref:`sediment_physics` section 8 for the Kean and Smith formulation
   that anugaSed used, which is not yet available here.
