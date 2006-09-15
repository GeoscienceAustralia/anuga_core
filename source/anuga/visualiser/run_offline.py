from offline import OfflineVisualiser

o = OfflineVisualiser("../../swollen_viewer/tests/cylinders.sww")
#o = OfflineVisualiser("../../../../swollen_viewer/tests/karratha_100m.sww")
o.render_quantity_height("elevation", dynamic=False)
o.render_quantity_height("stage", dynamic=True)
#o.colour_height_quantity('stage', (0.0, 0.0, 0.8))
o.colour_height_quantity('stage', (lambda q:q['stage'], 0, 10))
o.precache_height_quantities()
o.run()
