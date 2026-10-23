.. _use_paraview:

ParaView and the Web Viewer
===========================

The `anuga-paraview <https://github.com/anuga-community/anuga_paraview>`_
package opens ``.sww`` files in `ParaView <https://www.paraview.org>`_, and
in a web browser.  It reads the files with ``netCDF4`` and does not need
ANUGA, so it can be installed on a machine (or in a ParaView) without it.

.. contents:: Contents
   :local:
   :depth: 1

.. note::

   ``anuga-paraview`` is early work: the reader and the web viewer are
   prototypes.  Please report problems at
   https://github.com/anuga-community/anuga_paraview/issues.


ParaView
--------

ParaView is a general-purpose scientific visualisation application.  With
the ANUGA reader, a ``.sww`` file opens like any other data set, so
ParaView's filters (contours, slices, clips, calculators, ...), its Python
scripting and its batch rendering can all be used on ANUGA results.

Installation
~~~~~~~~~~~~

The reader needs ``netCDF4`` in ParaView's Python, so the conda-forge build
of ParaView is the easiest::

    conda create -n paraview -c conda-forge paraview netcdf4
    conda activate paraview
    pip install anuga-paraview

It is tested with ParaView 6.2.

Usage
~~~~~

::

    anuga-paraview results.sww

starts ParaView with the reader loaded and a ready-made view: the bed
coloured by elevation and the water by depth, with very shallow water
faded out.  Press Play on the animation toolbar.  Without a file name,
``.sww`` files open with File > Open.

The reader has two outputs, ``Bed`` (z = elevation) and ``Water``
(z = stage, with dry triangles left out).  Both carry every quantity in the
file, plus derived ``depth``, ``speed``, ``momentum`` and a ``velocity``
vector.  Its properties set the vertical exaggeration, switch to
georeferenced coordinates, and add the maxima over all timesteps
(``max_depth``, ``max_speed``, ...).

The plugin also adds an **ANUGA Timeseries** filter
(Filters > ANUGA > ANUGA Timeseries), which plots depth and speed against
time at a point.  It reads only the values at that point from the file, so
it is quick even for long runs.

To load the reader into a ParaView started some other way, use
Tools > Manage Plugins > Load New... and pick ``sww_reader.py`` in the
directory printed by ``anuga-paraview --plugin-path``.  For data on a
remote machine, ParaView's client/server mode (``pvserver``) can be used,
with the plugin loaded on both ends.

See the `anuga-paraview README
<https://github.com/anuga-community/anuga_paraview#readme>`_ for the
details.


Web viewer
----------

``anuga-web-viewer`` shows a ``.sww`` file in a web browser.  The server
reads the file and renders the images, and the browser receives only the
images, so a large result on a remote machine (an HPC or cloud node, say)
can be viewed without copying it.  It does not need ParaView, a display or
a GPU.

Installation
~~~~~~~~~~~~

::

    pip install "anuga-paraview[web]"

Usage
~~~~~

On the machine with the data::

    anuga-web-viewer results.sww --port 8080 --server

then open http://localhost:8080.  For a remote machine, forward the port
from your workstation and open the same address::

    ssh -L 8080:localhost:8080 remote-host

The viewer has a time slider and Play; the side panel sets the quantity
the water is coloured by (depth, speed, stage or momentum), the colour map
and scale maximum, the depth below which water is hidden, and the vertical
scale.  Clicking on the surface plots a timeseries at that point.  It
serves one user at a time.
