from builtins import range
from Numeric import array, Float
from anuga.file.netcdf import NetCDFFile
from visualiser import Visualiser
from vtk import vtkCellArray

class SWWVisualiser(Visualiser):
    '''
    A Visualiser that views SWW files. The source for this
    visualiser is a string containing the sww file name.
    '''

    def __init__(self, *args, **kwargs):
        Visualiser.__init__(self, *args, **kwargs)
        fin = NetCDFFile(self.vis_source, 'r')
        self.xPoints = array(fin.variables['x'], Float)
        self.yPoints = array(fin.variables['y'], Float)
        self.quantityCache = {}
        fin.close()

    def setupGrid(self):
        fin = NetCDFFile(self.vis_source, 'r')
        nTri = fin.variables['volumes'].shape[0]
        insertNextCell = vtkCellArray.InsertNextCell
        insertCellPoint = vtkCellArray.InsertCellPoint
        
        for v in range(nTri):
            insertNextCell(self.vtk_cells, 3)
            for i in range(3):
                insertCellPoint(self.vtk_cells, fin.variables['volumes'][v][i])

        fin.close()

    def getMaxFrameNumber(self):
        fin = NetCDFFile(self.vis_source, 'r')
        rv = fin.variables['time'].shape[0] - 1
        fin.close()
        return rv

    def getQuantityPoints(self, quantityName, dynamic=False):
        try:
            if dynamic:
                q = self.quantityCache[quantityName][self.vis_frame]
            else:
                q = self.quantityCache[quantityName]
        except KeyError:
            fin = NetCDFFile(self.vis_source, 'r')
            if dynamic:
                if quantityName not in self.quantityCache:
                    self.quantityCache[quantityName] = {}
                q = array(fin.variables[quantityName][self.vis_frame], Float)
                self.quantityCache[quantityName][self.vis_frame] = q
            else:
                q = array(fin.variables[quantityName], Float)
                self.quantityCache[quantityName] = q
            fin.close()
        return q

    def getQuantityDict(self):
        quantities = {}
        fin = NetCDFFile(self.vis_source, 'r')
        names = [ k for k in list(fin.variables.keys())
                  if k != 'x'
                  and k != 'y'
                  and k != 'z'
                  and k != 'time'
                  and k != 'volumes' ]
        argss = [ (name, len(fin.variables[name].shape) != 1) for name in names ]
        fin.close()
        for args in argss:
            quantities[name] = self.getQuantityPoints(*args)
        return quantities
