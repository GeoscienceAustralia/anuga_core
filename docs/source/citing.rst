.. _citing:

Citing ANUGA
============

If you use ANUGA in published work, please cite it.

Software
--------

Cite the ANUGA software using the metadata in the repository's ``CITATION.cff``
file — GitHub's **"Cite this repository"** button (top right of the
`repository page <https://github.com/anuga-community/anuga_core>`_) generates an
up-to-date APA or BibTeX entry from it. A minimal BibTeX entry:

.. code-block:: bibtex

   @software{anuga,
     title     = {{ANUGA}: hydrodynamic modelling of floods, storm surges and tsunamis},
     author    = {Roberts, Stephen and Davies, Gareth and Nielsen, Ole and others},
     publisher = {Australian National University and Geoscience Australia},
     url       = {https://github.com/anuga-community/anuga_core}
   }

From version 4.1.0, each release is also archived on Zenodo with a DOI of its
own, linked from the release on the
`GitHub releases page <https://github.com/anuga-community/anuga_core/releases>`_.
Cite that DOI to identify the exact version you used.

ANUGA was created by Geoscience Australia and the Mathematical Sciences Institute
at the Australian National University, and is now maintained by a community of
volunteers.

Key publications
----------------

Two papers describe the method ANUGA implements, and are the usual citations for
the model itself. The first develops the finite-volume solution of the
two-dimensional shallow water equations on unstructured triangular meshes that
ANUGA is built on; the second describes its application to coastal inundation:

* Zoppou, C. & Roberts, S. (1999). Catastrophic collapse of water supply
  reservoirs in urban areas. *Journal of Hydraulic Engineering*, 125(7),
  686–695. https://doi.org/10.1061/(ASCE)0733-9429(1999)125:7(686)

* Nielsen, O., Roberts, S., Gray, D., McPherson, A. & Hitchman, A. (2005).
  Hydrodynamic modelling of coastal inundation. In A. Zerger & R. M. Argent
  (Eds.), *MODSIM 2005 International Congress on Modelling and Simulation*
  (pp. 518–523). Modelling and Simulation Society of Australia and New Zealand.
  https://www.mssanz.org.au/modsim05/papers/nielsen.pdf

.. code-block:: bibtex

   @article{zoppou1999catastrophic,
     author  = {Zoppou, Christopher and Roberts, Stephen},
     title   = {Catastrophic Collapse of Water Supply Reservoirs in Urban Areas},
     journal = {Journal of Hydraulic Engineering},
     volume  = {125},
     number  = {7},
     pages   = {686--695},
     year    = {1999},
     doi     = {10.1061/(ASCE)0733-9429(1999)125:7(686)}
   }

   @inproceedings{nielsen2005hydrodynamic,
     author    = {Nielsen, Ole and Roberts, Stephen and Gray, Duncan and
                  McPherson, Andrew and Hitchman, Adrian},
     title     = {Hydrodynamic Modelling of Coastal Inundation},
     booktitle = {MODSIM 2005 International Congress on Modelling and Simulation},
     editor    = {Zerger, Andre and Argent, Robert M.},
     pages     = {518--523},
     year      = {2005},
     publisher = {Modelling and Simulation Society of Australia and New Zealand},
     url       = {https://www.mssanz.org.au/modsim05/papers/nielsen.pdf}
   }

Both are also cited from :doc:`mathematical_background`, where the numerical
method is described.

User manual
-----------

The ANUGA user manual describes the numerical method, the finite-volume
discretisation, and the validation suite, and has a DOI:

   https://dx.doi.org/10.13140/RG.2.2.17267.81446

Validation reports
------------------

The validation reports of each final release are archived on Zenodo as
versions of the *ANUGA validation reports* record. Each version has its own
DOI and is linked from the corresponding release on the
`GitHub releases page <https://github.com/anuga-community/anuga_core/releases>`_.
When a result depends on ANUGA reproducing one of the validation cases, cite
the version that matches the ANUGA release you used. See
:ref:`anuga-validation` for what the reports contain.

License
-------

ANUGA is free and open-source software, released under the
`Apache License, Version 2.0 <https://www.apache.org/licenses/LICENSE-2.0>`_.
You are free to use, modify, and redistribute it under the terms of that
licence; see ``LICENSE.txt`` in the repository for the full text.
