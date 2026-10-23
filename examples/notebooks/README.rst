Example Notebooks
=================

Interactive Jupyter notebooks demonstrating ANUGA usage from beginner to advanced level.

.. list-table::
   :header-rows: 1
   :widths: 40 45 15

   * - Notebook
     - Description
     - Level
   * - ``notebook_simple_example.ipynb``
     - Wave runup on a linearly sloping bed — core ANUGA workflow.
     - Beginner
   * - ``notebook_create_domain_from_regions.ipynb``
     - Unstructured mesh from polygon regions with tagged boundary segments.
     - Intermediate
   * - ``notebook_create_domain_with_riverwalls.ipynb``
     - Levee/riverwall structures via ``breaklines`` and ``create_riverwalls``.
     - Intermediate
   * - ``notebook_flooding_example.ipynb``
     - Real-world urban flood model of Merewether, Newcastle (NSW).
     - Advanced
   * - ``notebook_tsunami_benchmark.ipynb``
     - Okushiri Island tsunami runup (IWTS Benchmark 2).
     - Advanced
   * - ``notebook_sediment_transport.ipynb``
     - Suspended sediment: scour in a sloping channel, and settling in still
       water checked against its reference solution.
     - Intermediate

Data files required by the advanced notebooks are in ``anuga_core/examples/data/``.
The notebooks locate this directory automatically via ``anuga.__file__``.

These notebooks are also rendered in the online documentation at
https://anuga-core.readthedocs.io/en/latest/examples/index.html
