""" VTK sww Visualiser for ANUGA

python vtk_viewer swwfile.sww
"""

import sys, os
from anuga.visualiser import OfflineVisualiser

def get_filename():
    if len(sys.argv) > 1:
        filename = sys.argv[1]

        root, ext = os.path.splitext(filename)
            
        if ext != '.sww':
            print('WARNING: I only view sww files.' %filename)
      
        return filename



if __name__ == '__main__':

    
    filename = get_filename()
    if filename is not None:
        # The argument to OfflineVisualiser is the path to a sww file
        o = OfflineVisualiser(filename)

        # Specify the height-based-quantities to render.
        # Remember to set dynamic=True for time-varying quantities
        o.render_quantity_height("elevation", dynamic=False)
        o.render_quantity_height("stage", dynamic=True)

        # Colour the stage:
        # Either with an RGB value as a 3-tuple of Floats,
        # o.colour_height_quantity('stage', (0.0, 0.0, 0.8))
        # Or with a function of the quantities at that point, such as
        # the stage height:
        # 0 and 10 are the minimum and maximum values of the stage.
        o.colour_height_quantity('stage', (lambda q: q['stage'], 1.0, 5.0))
        # Or with the magnitude of the momentum at that point:
        # Needs the sqrt function from numeric. Again, 0 and 10
        # define the colour range.
        from numpy import sqrt
        #o.colour_height_quantity('stage',
        #                          (lambda q:sqrt((q['xmomentum'] ** 2) +
        #                                         (q['ymomentum'] ** 2)),
        #                                          0, 10))
        o.colour_height_quantity('stage',
                                 (lambda q: q['xmomentum'] /
                                  (q['stage'] - q['elevation'])),
                                 0, 5)


        # Draw some axes on the visualiser so we can see how big the wave is


        # Draw some axes on the visualiser so we can see how big the wave is
        o.render_axes()

        # Precaching the height-based quantities reduces the time taken
        # to draw each frame, but increases the time taken when the
        # visualiser starts.
        o.precache_height_quantities()

        # Start the visualiser (in its own thread).
        o.start()
        # Wait for the visualiser to terminate before shutdown
        o.join()
        
