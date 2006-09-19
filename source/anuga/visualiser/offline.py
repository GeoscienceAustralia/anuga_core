from Numeric import array, Float, ravel, zeros
from Scientific.IO.NetCDF import NetCDFFile
from Tkinter import Button, E, W
from visualiser import Visualiser
from vtk import vtkCellArray, vtkPoints, vtkPolyData

class OfflineVisualiser(Visualiser):
    """A VTK-powered offline visualiser which runs in its own thread.
    In addition to the functions provided by the standard visualiser,
    the following additional functions are provided:

    precache_height_quantities() - Precache all the vtkpoints
    structures for any dynamic height based quantities to render.
    """
    def __init__(self, source):
        """The source parameter is assumed to be a NetCDF sww file.
        """
        Visualiser.__init__(self, source)

        self.frameNumber = 0
        fin = NetCDFFile(self.source, 'r')
        self.maxFrameNumber = fin.variables['time'].shape[0] - 1
        fin.close()

        self.vtk_heightQuantityCache = []
        for i in range(self.maxFrameNumber + 1): # maxFrameNumber is zero indexed.
            self.vtk_heightQuantityCache.append({})

        self.paused = False
        
    def setup_grid(self):
        fin = NetCDFFile(self.source, 'r')
        self.vtk_cells = vtkCellArray()
        N_tri = fin.variables['volumes'].shape[0]
        for v in range(N_tri):
            self.vtk_cells.InsertNextCell(3)
            for i in range(3):
                self.vtk_cells.InsertCellPoint(fin.variables['volumes'][v][i])
        fin.close()
        
    def update_height_quantity(self, quantityName, dynamic=True):
        polydata = self.vtk_polyData[quantityName] = vtkPolyData()
        if dynamic is True:
            if not self.vtk_heightQuantityCache[self.frameNumber].has_key(quantityName):
                self.vtk_heightQuantityCache[self.frameNumber][quantityName]\
                    = self.read_height_quantity(quantityName, True, self.frameNumber);
            polydata.SetPoints(self.vtk_heightQuantityCache[self.frameNumber][quantityName])
        else:
            polydata.SetPoints(self.read_height_quantity(quantityName, False))
        polydata.SetPolys(self.vtk_cells)
            
    def read_height_quantity(self, quantityName, dynamic=True, frameNumber=0):
        """Read in a height based quantity from the NetCDF source file
        and return a vtkPoints object. frameNumber is ignored if
        dynamic is false."""
        fin = NetCDFFile(self.source, 'r')
        points = vtkPoints()
        if dynamic is True:
            N_vert = fin.variables[quantityName].shape[1]
        else:
            N_vert = len(fin.variables[quantityName])
        x = ravel(array(fin.variables['x'], Float))
        y = ravel(array(fin.variables['y'], Float))
        if dynamic is True:
            q = array(fin.variables[quantityName][frameNumber], Float)
        else:
            q = ravel(array(fin.variables[quantityName], Float))

        q *= self.height_zScales[quantityName]
        q += self.height_offset[quantityName]

        for v in range(N_vert):
            points.InsertNextPoint(x[v], y[v], q[v])
        fin.close()
        return points

    def precache_height_quantities(self):
        """Precache any height-based quantities. Call before rendering
        beigns."""
        for q in self.height_quantities:
            if self.height_dynamic[q] is True:
                print 'Precaching %s' % q
                for i in range(self.maxFrameNumber + 1): # maxFrameNumber is zero-indexed
                    print ' - Frame %d of %d' % (i, self.maxFrameNumber)
                    self.vtk_heightQuantityCache[i][q]\
                        = self.read_height_quantity(q, True, i)

    def build_quantity_dict(self):
        quantities = {}
        fin = NetCDFFile(self.source, 'r')
        for q in filter(lambda n:n != 'x' and n != 'y' and n != 'z' and n != 'time' and n != 'volumes', fin.variables.keys()):
            if len(fin.variables[q].shape) == 1: # Not a time-varying quantity
                quantities[q] = ravel(array(fin.variables[q], Float))
            else: # Time-varying, get the current timestep data
                quantities[q] = array(fin.variables[q][self.frameNumber], Float)
        fin.close()
        return quantities

    def setup_gui(self):
        Visualiser.setup_gui(self)
        self.tk_quit.grid(row=0, column=0, columnspan=6, sticky=W+E)
        self.tk_restart = Button(self.tk_controlFrame, text="<<<", command=self.restart)
        self.tk_restart.grid(row=1, column=0, sticky=W+E)
        self.tk_back10 = Button(self.tk_controlFrame, text="<<", command=self.back10)
        self.tk_back10.grid(row=1, column=1, sticky=W+E)
        self.tk_back = Button(self.tk_controlFrame, text="<", command=self.back)
        self.tk_back.grid(row=1, column=2, sticky=W+E)
        self.tk_pauseResume = Button(self.tk_controlFrame, text="Pause", command=self.pauseResume, width=15)
        self.tk_pauseResume.grid(row=1, column=3, sticky=W+E)
        self.tk_forward = Button(self.tk_controlFrame, text=">", command=self.forward)
        self.tk_forward.grid(row=1, column=4, sticky=W+E)
        self.tk_forward10 = Button(self.tk_controlFrame, text=">>", command=self.forward10)
        self.tk_forward10.grid(row=1, column=5, sticky=W+E)

        # Make the buttons stretch to fill all available space
        for i in range(6):
            self.tk_controlFrame.grid_columnconfigure(i, weight=1)

    def run(self):
        self.tk_root.after(100, self.animateForward)
        Visualiser.run(self)

    def restart(self):
        self.frameNumber = 0
        self.redraw_quantities(True)
        self.pause()

    def back10(self):
        if self.frameNumber - 10 >= 0:
            self.frameNumber -= 10
        else:
            self.frameNumber = 0
        self.redraw_quantities(True)
        self.pause()

    def back(self):
        if self.frameNumber > 0:
            self.frameNumber -= 1
            self.redraw_quantities(True)
            self.pause()

    def pauseResume(self):
        if self.paused is True:
            self.resume()
        else:
            self.pause()

    def pause(self):
        self.paused = True
        self.tk_pauseResume.config(text="Resume")

    def resume(self):
        self.paused = False
        self.tk_pauseResume.config(text="Pause")
        self.tk_root.after(100, self.animateForward)

    def forward(self):
        self.forward_step()
        self.pause()

    def forward_step(self):
        if self.frameNumber < self.maxFrameNumber:
            self.frameNumber += 1
            self.redraw_quantities(True)
        else:
            self.pause()

    def forward10(self):
        if self.frameNumber + 10 <= self.maxFrameNumber:
            self.frameNumber += 10
        else:
            self.frameNumber = self.maxFrameNumber
        self.redraw_quantities(True)
        self.pause()

    def animateForward(self):
        if self.paused is not True:
            self.forward_step()
            self.tk_root.after(100, self.animateForward)
