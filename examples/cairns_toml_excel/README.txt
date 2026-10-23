NOTE: SIMPLE DATA PREPARATION IS REQUIRED TO RUN THIS EXAMPLE -- SEE BELOW

NOTE: The code in example (in the current directory, and under 'setup') can
optionally be under a BSD3 licence (see the LICENCE file). These 
licences do not apply to the 'xlrd' library, which has its own BSD style licence
conditions, see './xlrd/licences.py'. 

This example illustrates an excel based interface to ANUGA, using the cairns example.

The interface is generic, and allows access to a limited subset of ANUGA's
functionality. It has been used to run both flood and tsunami type models.

In some instances it is easier for new users than writing a python script --
however this necessarily comes at the expense of some flexibility, and we do
not ever expect to support all of ANUGA's functionality in this interface.

More advanced users may however find the code structure provides a good basis
for writing their own scripts. Note that most details occur in scripts in 'setup',
while the run_model.py script is the main driver routine.

## How to run ##

The DEM raster for this example now lives in the shared data directory
'../data/cairns/' (cairns.asc and friends), so a single copy is reused by both
this legacy Excel example and the TOML example in '../run_toml/cairns/'. The
cairns_example.toml here already points at '../data/cairns/cairns.asc'. The
Excel workbook (cairns_excel.xlsx) still refers to 'cairns_initialcond/' for the
elevation; to run the Excel path, copy (or symlink) '../data/cairns/cairns.asc'
and 'cairns.prj' into a 'cairns_initialcond/' folder here first.

Then, it should be able to be run with:
> python run_model.py cairns_excel.xls
or (in parallel, with 6 cores in this example):
> mpirun -np 6 python run_model.py cairns_excel.xls

See the xls file for explanation of the configuration.

## TOML configuration (alternative to Excel) ##

A plain-text TOML configuration file is provided as an alternative to the
Excel interface:
> python run_model.py cairns_example.toml
> mpirun -np 6 python run_model.py cairns_example.toml

TOML files are human-readable, version-control friendly (clean diffs), and
require no external dependencies on Python 3.11+.  For older Python versions
install the 'tomli' package (pip install tomli).

See cairns_example.toml for a fully commented example of every available
setting.

## Post Processing ##

There are also various post-processing scripts (see them for details):
    flow_through_cross_sections.py
    gauge_export.py
    ipython_velocity_vector_plot.py
    make_anugaviewer_movie.py
    points_export.py
    raster_export.py

## NOTE ##
This example is for illustrative purposes, to show how to set up a model with
the excel interface. It is not a realistic case (obviously!). Also 'design'
decisions about the mesh resolution and structure, placement of boundary
conditions, friction, elevation data quality, etc have not been given high
scrutiny or quality control.  In a 'real' study I would probably move the
lateral boundaries further away from the region of interest, do convergence
testing to check the influence of mesh size, potentially use a more carefully
designed mesh, etc.  All those things could be done using this excel interface.

