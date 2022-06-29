

Introduction
~~~~~~~~~~~~

After starting Python, import the :code:`anuga` module with

>>> import anuga

To save repetition, in the documentation we assume that 
ANUGA has been imported this way.

If importing ANUGA fails, it means that Python cannot find the installed
module. Check your installation and your PYTHONPATH.

The following :class:`Domain` class is available:

:class:`Domain`
   This class initializes a domain object.

Initialize a ANUGA Model with

>>> domain = anuga.Domain()

Once a :class:`Domain` is initialized, there are several options available to 
setup the domain (initial conditions, boundary conditions, operators) and run the model (evolve). 

.. autosummary::
   :toctree: generated 
   
   anuga
 
