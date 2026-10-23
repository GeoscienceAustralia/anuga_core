.. _checkpointing:

.. currentmodule:: anuga

Checkpointing
=============

Long simulations can be interrupted unexpectedly — a power cut, a job-queue
time limit on an HPC cluster, or a network failure.  ANUGA's checkpointing
mechanism guards against losing all progress by periodically saving the full
domain state to disk.  If the run is interrupted, it can be restarted from
the most recent checkpoint rather than from time zero.

How it works
------------

At regular intervals (controlled by wall-clock time or by the number of
yield steps) ANUGA serialises the entire ``Domain`` object — quantities,
boundary conditions, operators, and simulation time — to a `pickle
<https://docs.python.org/3/library/pickle.html>`_ file.  On restart, that
pickle file is loaded and the evolve loop continues from where it stopped.

.. note::

   The `dill <https://pypi.org/project/dill/>`_ package is used when
   available (it can serialise a wider range of Python objects than the
   standard ``pickle`` module).  Install it with ``pip install dill`` or
   ``conda install dill`` for best results, especially when operators or
   user-defined functions are attached to the domain.

Checkpoint files are named::

   <checkpoint_dir>/<domain_name>_<simulation_time>.pickle

For parallel runs each rank writes its own file::

   <checkpoint_dir>/<domain_name>_P<nproc>_<rank>_<simulation_time>.pickle


Enabling checkpointing
-----------------------

Call ``domain.set_checkpointing()`` after the domain is fully configured
but before the evolve loop:

.. code-block:: python

   domain.set_checkpointing(
       checkpoint_dir='CHECKPOINTS',   # directory to store pickle files
       checkpoint_time=900,            # save every 900 s of wall time (15 min)
   )

.. list-table:: ``set_checkpointing`` parameters
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``checkpoint``
     - ``True``
     - Set to ``False`` to disable checkpointing without removing the call.
   * - ``checkpoint_dir``
     - ``'CHECKPOINTS'``
     - Directory where pickle files are written.  Created automatically if
       it does not exist.
   * - ``checkpoint_step``
     - ``10``
     - Save a checkpoint every *N* yield steps.  Ignored if
       ``checkpoint_time`` is set.
   * - ``checkpoint_time``
     - ``None``
     - Save a checkpoint every *N* seconds of **wall-clock** time.
       Overrides ``checkpoint_step`` when set.  Recommended for long runs
       where the yield step interval is unpredictable.

Use ``checkpoint_time`` rather than ``checkpoint_step`` for production runs
so that checkpoints are written at predictable wall-clock intervals
regardless of how long each yield step takes.


Restarting from a checkpoint
------------------------------

Use the try / except pattern shown below.  On a fresh run the ``try`` block
fails (no checkpoint files exist) and falls through to the normal domain
setup.  On a restart it succeeds and the evolve loop picks up from the last
saved simulation time.

.. code-block:: python

   import anuga
   from anuga import load_checkpoint_file

   DOMAIN_NAME    = 'my_simulation'
   CHECKPOINT_DIR = 'CHECKPOINTS'

   try:
       # Attempt to load the most recent checkpoint
       domain = load_checkpoint_file(domain_name=DOMAIN_NAME,
                                     checkpoint_dir=CHECKPOINT_DIR)
       print(f'Restarting from checkpoint at t = {domain.get_time():.1f} s')

   except Exception:
       # No checkpoint found — build the domain from scratch
       domain = anuga.rectangular_cross_domain(100, 100, len1=10.0, len2=10.0)
       domain.set_name(DOMAIN_NAME)

       domain.set_quantity('elevation', lambda x, y: -1.0 - 0.1 * x)
       domain.set_quantity('stage',     expression='elevation + 0.5')
       domain.set_quantity('friction',  0.03)

       Br = anuga.Reflective_boundary(domain)
       domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

       # Enable checkpointing — save every 15 minutes of wall time
       domain.set_checkpointing(checkpoint_dir=CHECKPOINT_DIR,
                                checkpoint_time=900)

   # Evolve — same code whether starting fresh or restarting
   for t in domain.evolve(yieldstep=60.0, finaltime=86400.0):
       domain.print_timestepping_statistics()

.. list-table:: ``load_checkpoint_file`` parameters
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``domain_name``
     - ``'domain'``
     - Base name of the domain (must match ``domain.set_name()``).
   * - ``checkpoint_dir``
     - ``'.'``
     - Directory to search for pickle files.
   * - ``time``
     - ``None``
     - Load the checkpoint at this specific simulation time.  By default
       the most recent checkpoint is used.


Parallel simulations
---------------------

Checkpointing works with MPI parallel runs.  Each rank reads its own
checkpoint file automatically — the ``_P<nproc>_<rank>`` suffix is appended
internally by ``load_checkpoint_file``.  The try / except pattern is the
same as in the serial case, but the domain setup in the ``except`` branch
must follow the normal parallel structure (``distribute`` on rank 0 etc.):

.. code-block:: python

   import anuga
   from anuga import distribute, myid, numprocs, finalize, barrier
   from anuga import load_checkpoint_file

   DOMAIN_NAME    = 'my_parallel_simulation'
   CHECKPOINT_DIR = 'CHECKPOINTS'

   try:
       domain = load_checkpoint_file(domain_name=DOMAIN_NAME,
                                     checkpoint_dir=CHECKPOINT_DIR)
   except Exception:
       if myid == 0:
           domain = anuga.rectangular_cross_domain(200, 200,
                                                   len1=10.0, len2=10.0)
           domain.set_name(DOMAIN_NAME)
           domain.set_quantity('elevation', lambda x, y: -1.0 - 0.1 * x)
           domain.set_quantity('stage', expression='elevation + 0.5')
       else:
           domain = None

       domain = distribute(domain)

       Br = anuga.Reflective_boundary(domain)
       domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

       domain.set_checkpointing(checkpoint_dir=CHECKPOINT_DIR,
                                checkpoint_time=900)

   barrier()

   for t in domain.evolve(yieldstep=60.0, finaltime=86400.0):
       if myid == 0:
           domain.print_timestepping_statistics()

   domain.sww_merge(delete_old=False)
   finalize()


Practical tips
--------------

**Choose checkpoint interval wisely**
   For a run expected to take several hours, saving every 15–30 minutes
   (``checkpoint_time=900`` to ``1800``) strikes a good balance between
   protection and disk overhead.  Checkpointing too frequently can slow
   the simulation noticeably.

**Checkpoint files accumulate**
   ANUGA does not delete old checkpoint files automatically.  Each
   checkpoint is a separate pickle file.  Clean up stale files manually
   once a run has completed successfully, keeping the ``CHECKPOINTS``
   directory from growing unbounded across multiple runs.

**SWW output is unaffected**
   Checkpointing saves domain state for restart purposes only.  The SWW
   output file continues to be written normally throughout the run, so
   results up to the interruption point are preserved in the SWW file
   even without a restart.

**Keep domain name consistent**
   ``load_checkpoint_file`` matches files by ``domain_name``.  If you
   change the name between runs the checkpoint will not be found.

A working example including MPI checkpointing is in
``examples/checkpointing/runCheckpoint.py``.

.. seealso::

   `ANUGA User Manual — Chapter 13: Checkpointing
   <https://github.com/anuga-community/anuga_user_manual>`_
   discusses the checkpointing mechanism in more detail, including how to
   manage checkpoint files across long cluster runs.
