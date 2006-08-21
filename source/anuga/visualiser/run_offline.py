from offline import OfflineVisualiser

o = OfflineVisualiser("../../../../swollen_viewer/tests/cylinders.sww")
o.render_quantity_height("elevation", dynamic=False)
o.render_quantity_height("stage", dynamic=True)
o.run()
