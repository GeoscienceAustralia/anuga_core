.. currentmodule:: anuga

Setting up Initial Conditions
=============================

The domain class incorporates a number of important quantities:

 - ``stage`` — absolute water-surface elevation (m)
 - ``elevation`` — bed elevation / bathymetry (m)
 - ``xmomentum`` — depth-integrated momentum in the x-direction (m²/s)
 - ``ymomentum`` — depth-integrated momentum in the y-direction (m²/s)
 - ``friction`` — Manning roughness coefficient *n* (dimensionless)

Water **depth** is the derived quantity ``stage - elevation`` (a cell is dry
where these are equal). Elevation and stage are measured from the same vertical
datum, so ``stage`` is *not* the water depth.

These variables are stored in the domain as quantities.
The quantities are stored as a dictionary with the key being the name of the 
quantity and the value being the quantity itself. They all have default values of 0.0.

The setting of the initial conditions is done by setting the values of these quantities.
The values can be set by using the :meth:`set_quantity <Domain.set_quantity>` 
method of the domain object.

For instance, to set the elevation to a function of x and y, and the stage to a constant value,
the following code can be used:

.. code-block:: python

    domain.set_quantity('elevation', function = lambda x,y : x/10)
    domain.set_quantity('stage', expression = "elevation + 0.2" )


The `set_quantity` method can also be used to set the initial conditions 
for the xmomentum, ymomentum and friction, indeed any quantity that is 
stored in the domain.

.. seealso::

   `ANUGA User Manual — Chapter 8: Initial Conditions and set_quantity
   <https://github.com/anuga-community/anuga_user_manual>`_
   covers ``set_quantity`` in depth, including raster file inputs, spatial
   averaging, expressions involving other quantities, and fitting point clouds
   onto the mesh.

.. seealso::

   :ref:`api_domain` — full API of ``Domain.set_quantity``.