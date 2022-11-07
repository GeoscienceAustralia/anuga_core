
import numpy as num
from tkinter import Button, E, Tk, W
from threading import Event
from .visualiser import Visualiser
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
        
    def setup_grid(self):
        self.vtk_cells = vtkCellArray()
        triangles = self.source.get_triangles()
        N_tri = len(self.source)
        verticies = self.source.get_vertex_coordinates()
        N_vert = len(verticies)

        # Also build vert_index - a list of the x & y values of each vertex
        self.vert_index = num.zeros((N_vert,2), float)
        for n in range(N_tri):
            self.vtk_cells.InsertNextCell(3)
            for v in range(3):
                self.vert_index[n * 3 + v] = verticies[n * 3 + v]
                self.vtk_cells.InsertCellPoint(n * 3 + v)

    def update_height_quantity(self, quantityName, dynamic=True):
        N_vert = len(self.source.get_vertex_coordinates())
        qty_index = num.zeros(N_vert, float)
        triangles = self.source.get_triangles()
        vertex_values, _ = self.source.get_quantity(quantityName).get_vertex_values(xy=False, smooth=False)

        
        for n in range(N_vert):
            qty_index[n] = vertex_values[n]

        points = vtkPoints()
        for v in range(N_vert):
            points.InsertNextPoint(self.vert_index[v][0],
                                   self.vert_index[v][1],
                                   qty_index[v] * self.height_zScales[quantityName]
                                   + self.height_offset[quantityName])
            if self.xmin is None or self.xmin > self.vert_index[v][0]:
                self.xmin = self.vert_index[v][0]
            if self.xmax is None or self.xmax < self.vert_index[v][0]:
                self.xmax = self.vert_index[v][0]
            if self.ymin is None or self.ymin > self.vert_index[v][1]:
                self.ymin = self.vert_index[v][1]
            if self.ymax is None or self.ymax < self.vert_index[v][1]:
                self.ymax = self.vert_index[v][1]
            if self.zmin is None or self.zmin > qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]:
                self.zmin = qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]
            if self.zmax is None or self.zmax < qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]:
                self.zmax = qty_index[v] * self.height_zScales[quantityName] + self.height_offset[quantityName]

        polydata = self.vtk_polyData[quantityName] = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(self.vtk_cells)

    def get_3d_bounds(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax]
        
    def build_quantity_dict(self):
        triangles = self.source.get_triangles()
        quantities = {}
        for q in self.source.get_quantity_names():
            #quantities[q], _ = self.source.get_quantity(q).get_vertex_values(xy=False)
            quantities[q]  = self.source.get_quantity(q).vertex_values.flatten()
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

    def update(self,pause=False):
        """Sync the visualiser to the domain. Call this in the evolve loop."""
            
        if self.running:
            self.sync_redrawReady.set()
            self.sync_idle.wait()
            self.sync_idle.clear()
            self.sync_unpaused.wait()

        if pause and self.running:
            if self.sync_unpaused.isSet():
                self.sync_unpaused.clear()
                self.tk_pauseResume.config(text="Resume")
                
                self.sync_redrawReady.set()
                self.sync_idle.wait()
                self.sync_idle.clear()
                self.sync_unpaused.wait()
            
        return self.running

    def evolveFinished(self):
        """Stop the visualiser from waiting on signals from the evolve loop.
        Call this just after the evolve loop to ensure a clean shutdown."""
        self.running = False
        self.sync_redrawReady.set()
