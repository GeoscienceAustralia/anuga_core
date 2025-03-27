.. currentmodule:: anuga

Setting up Initial Conditions
=============================

The domain class incorporates a number of important quantities

 - stage
 - elevation (or bed)
 - xmomentum (xmom)
 - ymomentum (ymom)
 - friction

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

Reference
---------

.. automethod:: Domain.set_quantity