# VTK-Based visualiser
# intended to replace realtime_visualisation_new.py

import threading
import Tkinter
import vtk
from Numeric import *
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget
#from vtk_visualiser_ext import make_vtkpoints

class Visualiser(threading.Thread):

    """A VTK-powered visualiser, designed to replace the VPython one.
    Intended to be run in its own thread.

    Customisation options for the visualiser are as follows:

    setup: Dictionary mapping quantity name -> boolean.
    if setup[q] is true, draw the quantity when the visualiser starts.

    updating: Dictionary mapping quantity name -> boolean.
    if updating[q] is true, update the rendering of this quantity each
    timestep.

    qcolor: Dictionary mapping quantity name -> (float, float, float)
    if the name of a quantity is found in qcolor, the colour specified (in
    (r, g, b) from 0.0 to 1.0) is used for display of that quantity instead of
    (0.5, 0.5, 0.5) (the default)

    scale_z: Dictionary mapping quantity name -> float
    Override the z scaling of this quantity with the float specified.

    default_scale_z: float
    multiply all z coordinates by this amount unless an overriding value
    exists in scale_z.
    """

    def __init__(self, domain, default_scale_z=1.0, title='Test'):
        threading.Thread.__init__(self)
        # Initialise data structures. setup and updating are maps
        # quantity name -> boolean, setup means to render it,
        # update means to update with time.
        # qcolor maps quantity name -> (float, float, float): the colour (r, g, b)
        self.setup = {}
        self.updating = {}
        self.qcolor = {}
        self.scale_z = {}
        self.coloring = {}
        self.default_scale_z = default_scale_z
        self.domain = domain
        self.vertices = domain.vertex_coordinates
        
        self.idle = threading.Event()
        self.redraw_ready = threading.Event()
        self.unpaused = threading.Event()
        self.unpaused.set()

        # Internal use - storage of vtk objects
        self.grids = {}
        self.scalars = {}
        self.actors = {}
        self.polydata = {}
        self.mappers = {}

        self.running = True

        # Default options
        for x in self.domain.quantities:
            self.setup[x] = False
            self.updating[x] = False
            self.coloring[x] = False
        self.start()
            
    def run(self):
        self.initialise_gui()
        self.add_axes()
        self.setup_all()
        self.root.after(100, self.idle.set)
        self.root.mainloop()

    def initialise_gui(self):
        """Prepare the GUI for the Visualiser, and set up the
        renderer.

        """
        # Using the TK VTK widget allows extendability
        # should a gui need adding.
        self.root = Tkinter.Tk()

        # Message handling with after
        self.root.after(100, self.redraw)

        self.renderWidget = vtkTkRenderWidget(self.root, width=400, height=400)
        self.renderWidget.pack(expand='true', fill='both')
        self.renderWindow = self.renderWidget.GetRenderWindow()
        self.renderer = vtk.vtkRenderer()
        self.renderWindow.AddRenderer(self.renderer)

        self.quitButton = Tkinter.Button(self.root, text='Quit', command=self.shutdown)
        self.quitButton.pack(side=Tkinter.BOTTOM)

        self.pauseButton = Tkinter.Button(self.root, text='Pause', command=self.unpaused.clear)
        self.pauseButton.pack(side=Tkinter.LEFT)

        self.resumeButton = Tkinter.Button(self.root, text='Resume', command=self.unpaused.set)
        self.resumeButton.pack(side=Tkinter.LEFT)
        
        #self.w2if = vtk.vtkWindowToImageFilter()
        #self.w2if.SetInput(self.renderWindow)
        
        #self.pnmWr = vtk.vtkPNMWriter() 
        #self.pnmWr.SetInput(self.w2if.GetOutput())


    def add_axes(self):
        """Add axes to this Visualiser - TODO
        """
        pass

    def setup_all(self):
        """Draw in the data that is specified by setup or update
        """

        self.N_tri = len(self.domain.triangles)
        self.N_vert = len(self.vertices)
        self.cells = vtk.vtkCellArray()
        self.vertices = self.domain.get_vertex_coordinates()
        self.vert_index = zeros((self.N_vert,2), Float)
        for n in range(self.N_vert):
            for i in range(3):
                self.vert_index[self.domain.triangles[n][i]] = self.vertices[n][i*2:i*2+2]

        # Prepare the list of cells
        for t in self.domain.triangles:
            self.cells.InsertNextCell(3)
            for i in range(3):
                self.cells.InsertCellPoint(t[i])

        # Set up the rendering of each quantity
        for q in self.domain.quantities:
            if self.setup[q] | self.updating[q]:
                self.draw_quantity(q)

    def draw_quantity(self, q):
        if self.scale_z.has_key(q):
            scale = self.scale_z[q]
        else:
            scale = self.default_scale_z
            
        if q == 'elevation':
            shift = 0.001
        else:
            shift = 0.0
            
        #############################################################
        # Disabled the make_vtkpoints call because the extension is
        # difficult to get to link under windows
        #############################################################
        #self.grids[q] = make_vtkpoints(self.N_tri,
        #                           self.N_vert,
        #                           scale,
        #                           self.domain.quantities[q].vertex_values,
        #                           self.vert_index,
        #                           self.domain.triangles)
        #grid = self.grids[q]
        #############################################################

        qty_index = zeros(self.N_vert, Float)
        for n in range(self.N_tri):
            for v in range(3):
                qty_index[self.domain.triangles[n][v]] = self.domain.quantities[q].vertex_values[n][v]

        self.grids[q] = vtk.vtkPoints()
        self.scalars[q] = vtk.vtkFloatArray()
        grid = self.grids[q]
        scalars = self.scalars[q]

        for v in range(self.N_vert):
            grid.InsertNextPoint(self.vert_index[v][0],
                                 self.vert_index[v][1],
                                 shift + qty_index[v] * scale)
            scalars.InsertNextValue(qty_index[v]);

        # Can't recycle vtkPolyData objects: Apparently they behave
        # unusually if the points (i.e. vertex data) is set after
        # the polys (i.e. triangle data)
        self.polydata[q] = vtk.vtkPolyData()
        polydata = self.polydata[q]

        polydata.SetPoints(grid)
        if self.coloring[q]:
            polydata.GetPointData().SetScalars(scalars);
        polydata.SetPolys(self.cells)

        if self.mappers.has_key(q):
            mapper = self.mappers[q]
            mapper.SetInput(polydata)
            mapper.SetScalarRange(0.0,0.5)
            mapper.Update()
        else:
            self.mappers[q] = vtk.vtkPolyDataMapper()
            mapper = self.mappers[q]
            mapper.SetInput(polydata)

        if not self.actors.has_key(q):
            self.actors[q] = vtk.vtkActor()
            actor = self.actors[q]
            actor.SetMapper(mapper)

            if self.qcolor.has_key(q):
                actor.GetProperty().SetColor(self.qcolor[q])
            else:
                actor.GetProperty().SetColor(0.5, 0.5, 0.5)

            self.renderer.AddActor(actor)

    def redraw(self):
        if self.redraw_ready.isSet():
            self.redraw_ready.wait()
            self.redraw_ready.clear()
            for q in self.domain.quantities:
                if self.updating[q]:
                    self.draw_quantity(q)

            self.renderWindow.Render()
             
            self.root.update_idletasks()
            self.idle.set()
        self.root.after(100, self.redraw)     
    
    def shutdown(self):
        """Shutdown the visualiser
        """
        self.running = False
        self.idle.set()
        self.unpaused.set()
        self.root.withdraw()
        self.root.quit()

    def update(self):
        """Update the visualiser's display.
        Clients are expected to call this in their evolve() loop,
        to keep the visualiser in sync with the simulation.
        """
        if self.running:
            self.redraw_ready.set()
            self.idle.wait()
            self.idle.clear()
            self.unpaused.wait()
