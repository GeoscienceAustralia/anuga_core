
from feature import Feature
from tkinter import Button
from vtk import vtkFloatArray, vtkPoints, vtkPolyData, vtkPolyDataMapper
from time import time
class HeightQuantity(Feature):
    '''
    A height quantity, such as stage, elevation, ...
    '''
    def __init__(self, quantityName, zScale=1.0,
                 offset=0.0, **kwargs):
        '''
        Parameters:
        quantityName: string - name of a quantity
        zScale: float - multiply point z-values by this
        offset: float - add this to point z-values
        '''
        Feature.__init__(self, **kwargs)

        self.quantityName = quantityName
        self.zScale = zScale
        self.offset = offset

    def button(self, tk_component):
        return Button(tk_component,
                      text='H: ' + self.quantityName)

    def draw(self, renderer):
        ### FIXME this should be made faster (C++ Module? How to deal with C++ linkage problems?)
        ### Sticking the vtkPoints objects in a cache would help somewhat but not on the first view.
        ### - Jack
        if not self.drawn:
            vtk_points = vtkPoints()
            points = self.visualiser.getQuantityPoints(self.quantityName, dynamic=self.dynamic)
            nPoints = len(points)
            vtk_points.SetNumberOfPoints(nPoints)
            setPoint = vtkPoints.SetPoint
            for i in range(nPoints):
                z = points[i] * self.zScale + self.offset
                setPoint(vtk_points, i, self.visualiser.xPoints[i], self.visualiser.yPoints[i], z)

            polyData = vtkPolyData()
            polyData.SetPoints(vtk_points)
            polyData.SetPolys(self.visualiser.vtk_cells)
            mapper = vtkPolyDataMapper()
            mapper.SetInput(polyData)
            setValue = vtkFloatArray.SetValue
            if hasattr(self.colour[0], '__call__'):
                scalars = self.colour[0](self.visualiser.getQuantityDict())
                nScalars = len(scalars)
                vtk_scalars = vtkFloatArray()
                vtk_scalars.SetNumberOfValues(nScalars)
                for i in range(nScalars):
                    setValue(vtk_scalars, i, scalars[i])
                polyData.GetPointData().SetScalars(vtk_scalars)
                mapper.SetScalarRange(self.colour[1:3])
            mapper.Update()
            self.actor.SetMapper(mapper)
        Feature.draw(self, renderer)
