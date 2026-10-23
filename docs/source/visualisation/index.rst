.. currentmodule:: anuga


Visualisation
=============

ANUGA provides several tools for visualising simulation results. During a
simulation run, the ``domain`` object can generate plots of the current state via
:ref:`Domain_plotter <use_domain_plotter>`. After the simulation, the ``.sww``
output is most quickly viewed — especially for large files — with the interactive
:ref:`ANUGA Viewer <use_anuga_viewer>`. Results can also be opened in the
``anuga_sww_gui`` GUI, explored programmatically via ``anuga.SWW_plotter``,
opened in ParaView or a web browser, or exported to GIS such as QGIS.


Choosing a viewer
-----------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Tool
     - Best for
     - Install
   * - :ref:`ANUGA Viewer <use_anuga_viewer>`
     - Fast, interactive 3D animation of large ``.sww`` files on a
       desktop; map tiles; recording camera paths for movies.
     - Separate application: installers from
       https://github.com/anuga-community/anuga-viewer/releases
   * - :ref:`anuga_sww_gui <sww_gui>`
     - 2D frames and animations (GIF/MP4) of any quantity, maximum
       envelopes, timeseries at picked points, basemaps.
     - ``pip install anuga[gui]``
   * - :ref:`SWW_plotter <use_sww_plotter>`
     - Plots from your own scripts or notebooks.
     - Part of ``anuga``
   * - :ref:`ParaView <use_paraview>`
     - Analysis in 3D: contours, slices, calculators, Python scripting,
       batch rendering, remote data through ``pvserver``.
     - ParaView, plus ``pip install anuga-paraview``
   * - :ref:`Web viewer <use_paraview>`
     - Viewing results on a remote or headless machine (HPC, cloud) in a
       browser, without copying the file.
     - ``pip install "anuga-paraview[web]"``
   * - :ref:`QGIS <use_qgis>`
     - Combining results with other GIS layers; mesh calculator;
       contours for maps.
     - QGIS
   * - :ref:`Domain_plotter <use_domain_plotter>`
     - Plots of the current state while a simulation runs.
     - Part of ``anuga``

.. only:: html

.. toctree::
   :maxdepth: 1

   use_anuga_viewer
   use_sww_gui
   use_domain_plotter
   use_sww_plotter
   use_paraview
   use_qgis

.. only:: html




