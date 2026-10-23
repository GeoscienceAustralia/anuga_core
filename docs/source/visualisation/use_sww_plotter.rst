.. _use_sww_plotter:

.. currentmodule:: anuga

SWW Plotter
===========

``SWW_plotter`` provides a Python interface for visualising ANUGA simulation
output files (``.sww`` format).  It loads the full time series from the SWW
file into memory so any frame can be plotted or saved without re-running the
simulation.

For plotting **during** a running simulation, use
:ref:`use_domain_plotter` instead.

Features
--------

- Plot and animate results from ANUGA SWW files: depth, stage, speed.
- Configurable colourmap, transparency, and optional OpenStreetMap basemap.
- Compute water volumes and cross-section flows from stored results.


Basic usage
-----------

.. code-block:: python

    # Run this code in a Jupyter notebook environment
    import anuga
    import matplotlib.pyplot as plt

    # Allow inline jshtml animations in Jupyter notebooks
    %matplotlib inline
    from matplotlib import rc
    rc('animation', html='jshtml')

    # Create SWW plotter object by opening an SWW file
    splotter = anuga.SWW_plotter('domain.sww')

    # Find min and max depth for plotting
    vmin = splotter.depth.min()
    vmax = splotter.depth.max()

    # Plot depth at last frame
    splotter.plot_depth_frame(-1, vmin=vmin, vmax=vmax)

    # Plot speed at second frame, collect fig, ax and change title
    fig, ax = splotter.plot_speed_frame(1)
    ax.set_title('Speed at second frame')

    # Plot mesh
    splotter.plot_mesh()

    # Save depth frames for all time steps
    for frame in range(len(splotter.time)):
        print(f'Processing frame {frame+1} of {len(splotter.time)}')
        splotter.save_depth_frame(frame, vmin=vmin, vmax=vmax)

    # Assemble saved frames into an animation
    anim = splotter.make_depth_animation()

    # Save animation to file
    anim.save('depth_animation.mp4', writer='ffmpeg', fps=5)

    # Display inline in Jupyter
    anim

    # Calculate water volume over entire domain
    volume = splotter.water_volume()
    print('Water volume at each time step:', volume)

    # Calculate flow through a polyline
    polyline = [[759195.0, 5912922.0],
                [759250.0, 5912892.0]]
    flow = splotter.get_flow_through_cross_section(polyline)
    print('Flow through cross section at each time step:', flow)


Visualisation options
---------------------

All ``plot_*`` and ``save_*`` frame methods accept these optional keyword
arguments in addition to ``frame``, ``figsize``, ``dpi``, ``vmin``, and
``vmax``:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Argument
     - Default
     - Description
   * - ``cmap``
     - ``'viridis'``
     - Any `matplotlib colormap
       <https://matplotlib.org/stable/gallery/color/colormap_reference.html>`_
       name, e.g. ``'plasma'``, ``'Blues'``, ``'jet'``.
   * - ``basemap``
     - ``False``
     - If ``True``, overlay a tile basemap behind the quantity plot.
       Requires the ``contextily`` package and the SWW file to contain a
       valid EPSG code.  Dry cells are left transparent so the map shows
       through.
   * - ``basemap_provider``
     - ``'OpenStreetMap.Mapnik'``
     - Dot-notation provider string.  Any key from
       ``anuga.utilities.animate.BASEMAP_PROVIDERS`` is accepted, or any
       ``contextily`` provider path.  Useful alternatives:
       ``'Esri.WorldImagery'`` (satellite),
       ``'Esri.WorldShadedRelief'`` (hillshade),
       ``'Esri.WorldTopoMap'`` (topographic),
       ``'OpenTopoMap'`` (OSM-based topo),
       ``'CartoDB.Positron'`` (light minimal).
   * - ``alpha``
     - ``1.0``
     - Opacity of the wet-area colour layer (0 = fully transparent,
       1 = fully opaque).  Values around ``0.6`` work well when
       ``basemap=True`` so the underlying map remains visible.


Example: custom colormap
------------------------

.. code-block:: python

    splotter = anuga.SWW_plotter('domain.sww')
    vmax = splotter.depth.max()

    for frame in range(len(splotter.time)):
        splotter.save_depth_frame(frame, vmin=0.0, vmax=vmax, cmap='plasma')

    anim = splotter.make_depth_animation()


Example: OSM basemap
--------------------

When the SWW file was produced from a georeferenced domain (i.e. it stores
an EPSG code), an OpenStreetMap basemap can be composited automatically.
Install ``contextily`` first::

    conda install contextily
    # or
    pip install contextily

.. code-block:: python

    splotter = anuga.SWW_plotter('domain.sww')
    vmax = splotter.depth.max()

    for frame in range(len(splotter.time)):
        splotter.save_depth_frame(frame, vmin=0.0, vmax=vmax,
                                  basemap=True, alpha=0.6, cmap='Blues')

    anim = splotter.make_depth_animation()

If ``contextily`` is not installed, or the SWW file has no EPSG code, a
warning is issued and the basemap is skipped — the frames are still
generated without it.


See Also
--------

.. autosummary::
   :toctree:

   SWW_plotter
