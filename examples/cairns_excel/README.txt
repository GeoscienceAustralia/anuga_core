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

To run it, you will first need to copy the files 'cairns.asc' and 'cairns.prj'
from '../cairns' to a folder named 'cairns_initialcond' in the current
directory. Note that 'cairns.asc' might be in the .zip file in '../cairns' and
require extraction first. We do not automate this step because the scripts are
generic (and we don't provide the data directly here, in order to reduce the
repository size).

Then, it should be able to be run with:
> python run_model.py cairns_excel.xls
or (in parallel, with 6 cores in this example):
> mpirun -np 6 python run_model.py cairns_excel.xls

See the xls file for explanation of the configuration.

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

