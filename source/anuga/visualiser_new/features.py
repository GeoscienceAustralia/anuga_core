from enthought.traits.api import Callable, Float, HasTraits, Instance, Range, Trait, true, Tuple
from vtk import vtkActor

class Feature(HasTraits):
    colour = Trait((0.5, 0.5, 0.5),
                   Tuple(Callable,
                         Float,
                         Float),
                   Tuple(Range(0.0, 1.0),
                         Range(0.0, 1.0),
                         Range(0.0, 1.0)))
    opacity = Range(0.0, 1.0, 1.0)
    actor = Instance(vtkActor, ())

class HeightFeature(Feature):
    z_scale = Range(0.0, None)
    offset = Float
    dynamic = true

class PolygonFeature(Feature):
    vertices = List(Tuple(Float, Float, Float))
