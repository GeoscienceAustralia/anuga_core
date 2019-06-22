"""
Towradgi Creek 17 August 1998 Storm Event Calibration
By Petar Milevski, some revisions by Gareth Davies
"""

import anuga
from os.path import join, isdir

from setup_boundaries import setup_boundaries
from setup_domain     import setup_domain
from setup_rainfall   import setup_rainfall
from setup_structures import setup_structures


if anuga.myid == 0 and not isdir('DEM_bridges'):
    msg = """
################################################################################
#
# Could not the find data directories
#
# You can download these directories using the data_download.py script.
# This will download over 120 MB of data!
#
################################################################################
"""
    raise Exception(msg)

#args = anuga.get_args()
#alg = args.alg
#verbose = args.verbose



towradgi = anuga.Simulation(setup_domain=setup_domain,
                      setup_boundaries=setup_boundaries,
                      setup_rainfall=setup_rainfall,
                      setup_structures=setup_structures)


print('checkpoint_time', towradgi.checkpoint_time)

towradgi.run()
