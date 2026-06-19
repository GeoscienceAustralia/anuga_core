
.. currentmodule:: anuga

Logging
=======

By default ANUGA prints progress messages to the terminal (stdout).
Calling :func:`set_logfile` activates file logging: from that point
every ``print()`` call and every ``log.*()`` call is written to both the
terminal **and** the named log file simultaneously.

Quick start
-----------

::

    import anuga

    anuga.set_logfile('my_run.log')

    # All subsequent output goes to terminal AND my_run.log
    print('Setting up domain ...')

    domain = anuga.rectangular_cross_domain(10, 5)
    ...

    for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
        domain.print_timestepping_statistics()

After this call the file ``my_run.log`` contains a complete record of the
run including timestep statistics, warnings, and any other printed output.


Output destinations
-------------------

ANUGA's output is split into three categories:

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Call
     - Terminal
     - Log file
     - When to use
   * - ``print(msg)``
     - ✓
     - ✓
     - User script milestones and results
   * - ``log.info(msg)``
     - ✓
     - ✓
     - Significant simulation events
   * - ``log.warning(msg)``
     - ✓
     - ✓
     - Non-fatal issues that need attention
   * - ``log.verbose(msg)``
     - ✗
     - ✓
     - Internal ANUGA chatter (mesh stats, solver steps)
   * - ``log.debug(msg)``
     - ✗
     - ✓
     - Developer diagnostics

``log.verbose()`` and ``log.debug()`` are below the default console
threshold (``INFO``) so they are silently dropped when no log file is
active.


Capturing verbose third-party output
-------------------------------------

Some ANUGA functions (mesh construction, domain distribution) produce
detailed ``print()`` output when called with ``verbose=True``.  Use the
:func:`file_only` context manager to redirect that output to the log
file without showing it on the terminal::

    import anuga
    import anuga.utilities.log as log

    anuga.set_logfile('my_run.log')

    with log.file_only():
        domain = anuga.create_domain_from_regions(
            bounding_polygon,
            boundary_tags=tags,
            maximum_triangle_area=res,
            verbose=True,   # full output goes to file, not screen
        )

Outside the ``with`` block normal tee behaviour resumes immediately.


Showing verbose output on the terminal
---------------------------------------

When debugging it can be useful to see all output on the terminal as
well.  Pass ``verbose_to_screen=True`` to :func:`set_logfile`::

    anuga.set_logfile('debug.log', verbose_to_screen=True)

This lowers the console threshold from ``INFO`` to ``DEBUG`` so
``log.verbose()`` and ``log.debug()`` messages also appear on screen.


Using log levels directly
--------------------------

The ``anuga.utilities.log`` module is available as ``anuga.log`` and
exposes the standard Python logging levels::

    import anuga.utilities.log as log

    log.info('Mesh built successfully')
    log.warning('Elevation data outside domain extent — using zero')
    log.verbose('Triangle count: %d' % n)   # file only
    log.debug('CG solver iteration %d, residual %.2e' % (k, r))

Levels in order of decreasing severity:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Level
     - Constant
   * - ``CRITICAL``
     - ``log.CRITICAL`` (50) — fatal errors
   * - ``ERROR``
     - ``log.ERROR``    (40) — recoverable errors
   * - ``WARNING``
     - ``log.WARNING``  (30) — non-fatal issues
   * - ``INFO``
     - ``log.INFO``     (20) — normal progress (default console threshold)
   * - ``DEBUG`` / ``VERBOSE``
     - ``log.DEBUG``    (10) — detail (default file threshold)


API reference
-------------

.. autofunction:: anuga.set_logfile

.. autoclass:: anuga.utilities.log.TeeStream
   :members: write, flush, close

.. autofunction:: anuga.utilities.log.file_only

.. autofunction:: anuga.utilities.log.verbose
.. autofunction:: anuga.utilities.log.info
.. autofunction:: anuga.utilities.log.warning
.. autofunction:: anuga.utilities.log.debug
.. autofunction:: anuga.utilities.log.critical
