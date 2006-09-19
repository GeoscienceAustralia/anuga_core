from Numeric import Float, zeros
from Tkinter import Button, E, W
from threading import Event
from visualiser import Visualiser
from vtk import vtkCellArray, vtkPoints, vtkPolyData

class RealtimeVisualiser(Visualiser):
    """A VTK-powered realtime visualiser which runs in its own thread.
    """
    def __init__(self, source):
        """The source parameter is assumed to be a Domain.
        """
        Visualiser.__init__(self, source)

        self.running = True

        # Synchronisation Constructs
        self.sync_idle = Event()
        self.sync_idle.clear()
        self.sync_unpaused = Event()
        self.sync_unpaused.set()
        self.sync_redrawReady = Event()
        self.sync_redrawReady.clear()

    def run(self):
        self.tk_root.after(100, self.sync_idle.set)
        Visualiser.run(self)

    def setup_grid(self):
        self.vtk_cells = vtkCellArray()
        # Also build vert_index - a list of the x & y values of each vertex
        N_vert = len(self.source.vertex_coordinates)
        self.vert_index = zeros((N_vert,2), Float)
        for n in range(N_vert):
            self.vtk_cells.InsertNextCell(3)
            for v in range(3):
                self.vert_index[self.source.triangles[n][v]] = self.source.vertex_coordinates[n][i*2:i*2+2]
                self.vtk_cells.InsertCellPoint(self.source.triangles[n][v])

    def update_height_quantity(self, quantityName, dynamic=True):
        N_vert = len(self.source.vertex_coordinates)
        qty_index = zeros(N_vert, Float)

        for n in range(len(self.source.triangles)):
            for v in range(3):
                qty_index[self.source.triangles[n][v]] = self.source.quantities[quantityName].vertex_values[n][v]

        points = vtkPoints()
        for v in range(N_vert):
            points.InsertNextPoint(self.vert_index[v][0],
                                   self.vert_index[v][1],
                                   qty_index[v] * self.height_zScales[quantityName]
                                   + self.height_offset[quantityName])

        polydata = self.vtk_polyData[quantityName] = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(self.vtk_cells)
        
    def build_quantity_dict(self):
        N_vert = len(self.source.vertex_coordinates)
        quantities = {}
        for q in self.source.quantities.keys():
            quantities[q] = zeros(N_vert, Float)
            for n in range(len(self.source.triangles)):
                for v in range(3):
                    quantities[q][self.source.triangles[n][v]] = self.source.quantities[q].vertex_values[n][v]
        return quantities

    def setup_gui(self):
        Visualiser.setup_gui(self)
        self.tk_pauseResume = Button(self.tk_controlFrame, text="Pause", command=self.pauseResume)
        self.tk_pauseResume.grid(row=1, column=0, sticky=E+W)

    def pauseResume(self):
        if self.sync_unpaused.isSet():
            self.sync_unpaused.clear()
            self.tk_pauseResume.config(text="Resume")
        else:
            self.sync_unpaused.set()
            self.tk_pauseResume.config(text="Pause")

    def shutdown(self):
        Visualiser.shutdown(self)
        self.running = False
        self.sync_idle.set()
        self.sync_unpaused.set()

