#--------------------------------
# import modules
#--------------------------------

import anuga
from anuga.validation_utilities import produce_report

args = anuga.get_args()

produce_report('radial_dam_break.py', args=args)


