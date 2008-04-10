#!/usr/bin/env python
from sww_visualiser import SWWVisualiser
from height_quantity import HeightQuantity

vis = SWWVisualiser(source='../../anuga_viewer/tests/cylinders.sww', recording=False, recordPattern='%02d.png')
vis.add_feature(HeightQuantity('elevation'))
vis.add_feature(HeightQuantity('stage', dynamic=True, #colour=(lambda q:q['stage'], 0.0, 10.0), offset=-0.01))
                               colour=(0.0, 0.0, 0.8)))

import cProfile
cProfile.run('vis.run()')
