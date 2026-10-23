.. _quickstart:

.. currentmodule:: anuga

Quick start
===========

.. code-block:: bash

   # 1. install (conda-forge, recommended)
   conda create -n anuga -c conda-forge anuga
   conda activate anuga

.. code-block:: python

   # 2. run a minimal model
   import anuga
   domain = anuga.rectangular_cross_domain(10, 5)
   domain.set_quantity('elevation', lambda x, y: -x / 10)
   domain.set_quantity('stage', 0.0)
   domain.set_boundary({b: anuga.Reflective_boundary(domain)
                        for b in domain.get_boundary_tags()})
   for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
       print(domain.timestepping_statistics())

Then view the ``.sww`` output (see :doc:`Visualisation <visualisation/index>`).
The sections below build this up step by step; new to ANUGA? Start with
:doc:`Examples <examples/index>`.
