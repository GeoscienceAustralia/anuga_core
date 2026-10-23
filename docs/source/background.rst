Background
==========

Modelling the effects on the built environment of natural hazards such
as riverine flooding, storm surges and tsunami is critical for
understanding their economic and social impact on our urban
communities.  Geoscience Australia and the Australian National
University have developed a hydrodynamic inundation modelling tool
called ANUGA to help simulate the impact of these hazards.

The core of ANUGA is the fluid dynamics object, called :code:`anuga.Domain`,
which is based on a finite-volume method for solving the Shallow Water
Wave Equation.  The study area is represented by a mesh of triangular
cells.  By solving the governing equation within each cell, water
depth and horizontal momentum are tracked over time.

A major capability of ANUGA is that it can model the process of
wetting and drying as water enters and leaves an area.  This means
that it is suitable for simulating water flow onto a beach or dry land
and around structures such as buildings.  ANUGA is also capable
of modelling hydraulic jumps due to the ability of the finite-volume
method to accommodate discontinuities in the solution and the bed.

To set up a particular scenario the user specifies the geometry
(bathymetry and topography), the initial water level (stage),
boundary conditions such as tide, and any operators that may
drive the system such as rainfall, abstraction of water, erosion or culverts.
See :doc:`setup_anuga_script/operators` for the operators and structures
available in ANUGA.

The procedure :code:`anuga.create_domain_from_regions` lets the user set up the
geometry of the problem from polygon regions and identify boundary segments and
regions using symbolic tags.  These tags may then be used to set the
actual boundary conditions and attributes for different regions
(e.g. the Manning friction coefficient) for each simulation.

Most ANUGA components are written in the object-oriented programming
language Python.  Software written in Python can be produced quickly
and can be readily adapted to changing requirements throughout its
lifetime.  Computationally intensive components are written for
efficiency in :code:`C` routines working directly with Python :code:`numpy`
structures.

ANUGA results are stored as ``.sww`` (NetCDF) files.  For fast, interactive
viewing of large ``.sww`` files the recommended tool is the OpenSceneGraph-based
`ANUGA Viewer <https://anuga-viewer.readthedocs.io/en/latest/>`_.  Results can
also be explored with the built-in ``anuga_sww_gui`` viewer, the ``SWW_plotter``
and ``Domain_plotter`` Python plotters, or by exporting to GIS such as QGIS.  See
:doc:`visualisation/index` for the visualisation options.

See :doc:`mathematical_background` for a full derivation of the governing
equations, the finite volume discretisation, and the flux and slope limiting
schemes used by ANUGA.

.. _acknowledgements:

Acknowledgements
----------------

ANUGA was created by Geoscience Australia and the Mathematical Sciences
Institute at the Australian National University, and is developed and
maintained by a community of volunteers.

The GPU and multicore solver introduced in ANUGA 4.0 (the unified compute
mode and its OpenMP target offloading, see :ref:`compute_modes` and
:ref:`use_gpu_offloading`) was developed in collaboration with the
`Centre for Development of Advanced Computing (C-DAC) <https://www.cdac.in/>`_,
India, and the `National Computational Infrastructure (NCI) <https://nci.org.au/>`_,
Australia. The C-DAC collaboration is part of a catchment flood prediction
project aimed at two-day flood forecasts, and grew out of GPU hackathons held
from 2023; NCI's Gadi supercomputer provided the GPU systems on which the work
was developed, tested and benchmarked. We thank Samir Shaikh, Shweta Das and Rutvik
Gulhane (C-DAC) and Jorge Luis Gálvez Vallejo (NCI) for their
contributions, and both organisations for their support.
