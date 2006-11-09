from Numeric import Float, zeros
from Tkinter import Button, E, Tk, W
from threading import Event
from visualiser import Visualiser
from vtk import vtkCellArray, vtkPoints, vtkPolyData

class RealtimeVisualiser(Visualiser):
    """A VTK-powered realtime visualiser which runs in its own thread.
    In addition to the functions provided by the standard visualiser,
    the following additional functions are provided:

    update() - Sync the visualiser to the current state of the model.
    Should be called inside the evolve loop.

    evolveFinished() - Clean up synchronisation constructs that tie the
    visualiser to the evolve loop. Call this after the evolve loop finishes
    to ensure a clean shutdown.
    """
    def __init__(self, source):
        """The source parameter is assumed to be a Domain.
        """
        Visualiser.__init__(self, source)

        self.running = True

        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None

        # Synchronisation Constructs
        self.sync_idle = Event()
        self.sync_idle.clear()
        self.sync_unpaused = Event()
        self.sync_unpaused.set()
        self.sync_redrawReady = Event()
        self.sync_redrawReady.clear()

    def run(self):
        self.alter_tkroot(Tk.after, (100, self.sync_idle.set))
        Visualiser.run(self)

    def setup_grid(self):
        self.vtk_cells = vtkCellArray()
        triangles = self.source.triangles
        N_tri = self.source.number_of_triangles
        verticies = self.source.get_vertex_coordinates()
        N_vert = len(verticies)
        # Also build vert_index - a list of the x & y values of each vertex
        self.vert_index = zeros((N_vert,2), Float)
        for n in range(N_tri):
            self.vtk_cells.InsertNextCell(3)
            for v in range(3):
                self.vert_index[triangles[n][v]] = verticies[n * 3 + v]
                self.vtk_cells.InsertCellPoint(triangles[n][v])

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
            if self.xmin == None or self.xmin > self.vert_index[v][0]:
                self.xmin = self.vert_index[v][0]
            if self.xmax == None or self.xmax < self.vert_index[v][0]:
                self.xmax = self.vert_index[v][0]
            if self.ymin == None or self.ymin > self.vert_index[v][1]:
                self.ymin = self.vert_index[v][1]
            if self.ymax == None or self.ymax < self.vert_index[v][1]:
                self.ymax = self.vert_index[v][1]
            if self.zmin == None or self.zmin > qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]:
                self.zmin = qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]
            if self.zmax == None or self.zmax < qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]:
                self.zmax = qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]

        polydata = self.vtk_polyData[quantityName] = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(self.vtk_cells)

    def get_3d_bounds(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax]
        
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

    def redraw(self):
        if self.running and self.sync_unpaused.isSet():
            self.sync_redrawReady.wait()
            self.sync_redrawReady.clear()
            self.redraw_quantities()
            self.sync_idle.set()
        Visualiser.redraw(self)

    def update(self):
        """Sync the visualiser to the domain. Call this in the evolve loop."""
        if self.running:
            self.sync_redrawReady.set()
            self.sync_idle.wait()
            self.sync_idle.clear()
            self.sync_unpaused.wait()

    def evolveFinished(self):
        """Stop the visualiser from waiting on signals from the evolve loop.
        Call this just after the evolve loop to ensure a clean shutdown."""
        self.running = False
        self.sync_redrawReady.set()
