from builtins import object
from types import FloatType
from vtk import vtkActor
class Feature(object):
    def __init__(self, colour=(0.5, 0.5, 0.5), opacity=1.0, dynamic=False):
        '''
        Parameters:
        colour: (float, float, float) - apply a single colour to the feature.
        opacity: float - 1.0 is opaque, 0.0 is invisible
        dynamic: boolean - this quantity changes with time
        '''
        self.actor = vtkActor()
        self.colour = colour
        self.drawn = False
        self.dynamic = dynamic
        self.opacity = opacity
        self.inRenderer = False
        self.visualiser = None
        
    def button(self, tk_component):
        '''
        Construct and return a Tkinter button that allows editing of
        the feature's parameters.
        '''
        raise NotImplementedError('Subclasses must override Feature::button!')

    def draw(self, renderer):
        '''
        Draw this object into the renderer, updating it if necessary.
        '''
        self.drawn = True
        if not self.inRenderer:
            self.inRenderer = True
            if type(self.colour[0]) is FloatType:
                self.actor.GetProperty().SetColor(self.colour)
            renderer.AddActor(self.actor)
        self.actor.GetProperty().SetOpacity(self.opacity)

    def redraw(self, renderer):
        '''
        Force a redraw of this feature.
        '''
        self.drawn = False
        self.draw(renderer)
