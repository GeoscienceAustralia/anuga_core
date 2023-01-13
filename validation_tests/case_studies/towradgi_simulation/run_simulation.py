"""
Generic simulation run script to load problem defined in 
the functions 

setup_domain
setup_boundaries
setup_rainfall
setup_structures

and arguments and parameters defined in the script project.py
"""

import anuga

from setup_boundaries import setup_boundaries
from setup_domain     import setup_domain
from setup_rainfall   import setup_rainfall
from setup_structures import setup_structures

sim = anuga.Simulation(setup_domain=setup_domain,
                            setup_boundaries=setup_boundaries,
                            setup_rainfall=setup_rainfall,
                            setup_structures=setup_structures)


print('checkpoint_time', sim.checkpoint_time)

sim.run()
