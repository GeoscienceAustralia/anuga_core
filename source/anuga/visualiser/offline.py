from Numeric import array, Float, zeros
from Scientific.IO.NetCDF import NetCDFFile
from Tkinter import Button, E, W
from visualiser import Visualiser
from vtk import vtkCellArray, vtkPoints, vtkPolyData

class OfflineVisualiser(Visualiser):
    """A VTK-powered offline visualiser which runs in its own thread.
    """
    def __init__(self, source):
        """The source parameter is assumed to be a NetCDF sww file.
        """
        Visualiser.__init__(self, source)

        self.frame_number = 0
        fin = NetCDFFile(self.source, 'r')
        self.max_frame_number = fin.variables['time'].shape[0] - 1
        fin.close()

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
        fin = NetCDFFile(self.source, 'r')
        if(fin.variables.has_key(quantityName)):
            points = vtkPoints()
            if dynamic:
                N_vert = fin.variables[quantityName].shape[1]
            else:
                N_vert = len(fin.variables[quantityName])
            x = array(fin.variables['x'], Float)
            y = array(fin.variables['y'], Float)
            if dynamic is True:
                q = array(fin.variables[quantityName][self.frame_number], Float)
            else:
                q = array(fin.variables[quantityName], Float)

            q *= self.height_zScales[quantityName]
            q += self.height_offset[quantityName]

            for v in range(N_vert):
                points.InsertNextPoint(x[v], y[v], q[v])
            polydata = self.vtk_polyData[quantityName] = vtkPolyData()
            polydata.SetPoints(points)
            polydata.SetPolys(self.vtk_cells)
        else:
            self.height_quantities.remove(quantityName)
        fin.close()

    def setup_gui(self):
        Visualiser.setup_gui(self)
        self.tk_renderWidget.grid(row=0, column=0, columnspan=6)
        self.tk_quit.grid(row=2, column=0, columnspan=6, sticky=W+E)
        self.tk_restart = Button(self.tk_root, text="<<<", command=self.restart)
        self.tk_restart.grid(row=1, column=0, sticky=W+E)
        self.tk_back10 = Button(self.tk_root, text="<<", command=self.back10)
        self.tk_back10.grid(row=1, column=1, sticky=W+E)
        self.tk_back = Button(self.tk_root, text="<", command=self.back)
        self.tk_back.grid(row=1, column=2, sticky=W+E)
        self.tk_pauseResume = Button(self.tk_root, text="Pause", command=self.pauseResume)
        self.tk_pauseResume.grid(row=1, column=3, sticky=W+E)
        self.tk_forward = Button(self.tk_root, text=">", command=self.forward)
        self.tk_forward.grid(row=1, column=4, sticky=W+E)
        self.tk_forward10 = Button(self.tk_root, text=">>", command=self.forward10)
        self.tk_forward10.grid(row=1, column=5, sticky=W+E)

    def restart(self):
        self.frame_number = 0
        self.redraw_quantities(True)

    def back10(self):
        if self.frame_number - 10 >= 0:
            self.frame_number -= 10
        else:
            self.frame_number = 0
        self.redraw_quantities(True)

    def back(self):
        if self.frame_number > 0:
            self.frame_number -= 1
            self.redraw_quantities(True)

    def pauseResume(self):
        print "Pause/Resume"

    def forward(self):
        if self.frame_number < self.max_frame_number:
            self.frame_number += 1
            self.redraw_quantities(True)

    def forward10(self):
        if self.frame_number + 10 <= self.max_frame_number:
            self.frame_number += 10
        else:
            self.frame_number = self.max_frame_number
        self.redraw_quantities(True)
