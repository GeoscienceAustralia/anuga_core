"""
Towradgi Creek 17 August 1998 Storm Event Calibration
By Petar Milevski, some revisions by Gareth Davies
"""

from simulation import Simulation
import anuga
from os.path import join

from setup_boundaries import setup_boundaries
from setup_domain import setup_domain
from setup_rainfall import setup_rainfall
from setup_structures import setup_structures



#args = anuga.get_args()
#alg = args.alg
#verbose = args.verbose



towradgi = Simulation(setup_domain=setup_domain, 
                      setup_boundaries=setup_boundaries,
                      setup_rainfall=setup_rainfall,
                      setup_structures=setup_structures)

towradgi.run()

