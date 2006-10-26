#!/usr/bin/env python
##########
# Demonstration of the VTK realtime Visualiser.
# Based on run_sw_rectangle.py
# Jack Kelly and Stephen Roberts
# October 2006
##########

# Import the offline visualiser
from anuga.visualiser import RealtimeVisualiser
from vtk import vtkCubeAxesActor2D

import time
from Numeric import array
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular


class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.75, y0=0.0, y1=1.0, h=5.0, h0=0.0):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.h  = h
        self.h0 = h0

    def __call__(self, x, y):
        return self.h0 + self.h*((x>self.x0)&(x<self.x1)&(y>self.y0)&(y<self.y1))

M = 30
points, vertices, boundary = rectangular(M, M, len1 = 1.0, len2 = 1.0)

yieldstep = 0.002
finaltime = 0.8
rect = [0.0, 0.0, 1.0, 1.0]

domain = Domain(points, vertices, boundary)

# Turn on the visualisation. The argument to the realtime visualier
# is a domain object.
v = RealtimeVisualiser(domain)

# Specify the height-based-quantities to render.
# Remember to set dynamic=True for time-varying quantities
v.render_quantity_height("elevation", dynamic=False)
v.render_quantity_height("stage", dynamic=True)

# Colour the stage:
# Either with an RGB value as a 3-tuple of Floats,
#v.colour_height_quantity('stage', (0.0, 0.0, 0.8))
# Or with a function of the quantities at that point, such as the stage height:
# 0 and 1 are the minimum and maximum values of the stage.
v.colour_height_quantity('stage', (lambda q:q['stage'], 0, 1))
# Or with the magnitude of the momentum at that point:
# Needs the sqrt function from Numeric. Again, 0 and 10 define the colour range.
#v.colour_height_quantity('stage', (lambda q:sqrt((q['xmomentum'] ** 2) +
#                                                 (q['ymomentum'] ** 2)), 0, 10))

# Draw some axes on the visualiser so we can see how big the wave is
v.render_axes()

# Increase the number of labels on the axes
v.alter_axes(vtkCubeAxesActor2D.SetNumberOfLabels, (5,))

# Draw a yellow polygon at height 2
v.overlay_polygon([(0, 0), (0, 0.1), (0.1, 0)], 2, colour=(1.0, 1.0, 0.0))

# Start the visualiser (in its own thread).
v.start()

#-----------------------------------------------------------------
# Boundaries and Initial conditions
#-----------------------------------------------------------------
R = Reflective_boundary(domain)
domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))

#-----------------------------------------------------------------
# Evolve
#-----------------------------------------------------------------
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    v.update()
    domain.write_time()
# Unhook the visualiser from the evolve loop.
# It won't shutdown cleanly unless you do this.
v.evolveFinished()

print 'That took %.2f seconds' %(time.time()-t0)

# Wait for the visualiser to be closed
v.join()
