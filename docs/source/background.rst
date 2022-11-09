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
boundary conditions such as tide, and any operators  that may
drive the system such as rainfall, abstraction of water,  erosion, culverts
See section :doc:`operators` for details of operators available in ANUGA.

The built-in mesh generator, called :code:`graphical_mesh_generator`, or 
the procedure :code:`anuga.create_domain_from_regions`
allows the user to set up the geometry
of the problem and to identify boundary segments and
regions using symbolic tags.  These tags may then be used to set the
actual boundary conditions and attributes for different regions
(e.g. the Manning friction coefficient) for each simulation.

Most ANUGA components are written in the object-oriented programming
language Python.  Software written in Python can be produced quickly
and can be readily adapted to changing requirements throughout its
lifetime.  Computationally intensive components are written for
efficiency in :code:`C` routines working directly with Python :code:`numpy`
structures.

The visualisation tool developed for ANUGA is based on
OpenSceneGraph, an Open Source Software (OSS) component allowing high
level interaction with sophisticated graphics primitives.
See \cite{nielsen2005} for more background on ANUGA.