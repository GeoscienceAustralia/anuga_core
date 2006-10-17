from threading import Thread
from Tkinter import Tk, Button, Frame, N, E, S, W
from types import FunctionType, TupleType
from vtk import vtkActor, vtkCubeAxesActor2D, vtkDelaunay2D, vtkFloatArray, vtkPoints, vtkPolyData, vtkPolyDataMapper, vtkRenderer
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget

class Visualiser(Thread):
    """Superclass of both the realtime and offline VTK visualisers
    """
    def __init__(self, source):
        Thread.__init__(self)

        self.source = source

        # Structures for Height Based quantities
        self.height_quantities = []
        self.height_zScales = {}
        self.height_dynamic = {}
        self.height_offset = {}

        # Structures for colouring quantities
        self.colours_height = {}

        # Structures used for VTK
        self.vtk_actors = {}
        self.vtk_axesSet = False
        self.vtk_drawAxes = False
        self.vtk_mappers = {}
        self.vtk_polyData = {}
        self.vtk_renderer = vtkRenderer()

        self.setup_gui()

    def run(self):
        self.setup_grid()
        # Draw Height Quantities
        for q in self.height_quantities:
            self.update_height_quantity(q, self.height_dynamic[q])
            self.draw_height_quantity(q)
        self.tk_root.mainloop()

    def redraw_quantities(self, dynamic_only=False):
        """Redraw all dynamic quantities, unless dynamic_only is True.
        """
        # Height quantities
        for q in self.height_quantities:
            if (dynamic_only is False) or (self.height_dynamic[q]):
                self.update_height_quantity(q, self.height_dynamic[q])
                self.draw_height_quantity(q)
        if self.vtk_drawAxes is True:
            self.vtk_axes.SetBounds(self.get_3d_bounds())
            if not self.vtk_axesSet:
                self.vtk_axesSet = True
                self.vtk_axes.SetCamera(self.vtk_renderer.GetActiveCamera())
                self.vtk_renderer.AddActor(self.vtk_axes)
        
    # --- Height Based Rendering --- #

    def render_axes(self):
        """Intstruct the visualiser to render cube axes around the render.
        """
        self.vtk_drawAxes = True
        self.vtk_axes = vtkCubeAxesActor2D()
        
    def get_axes(self):
        """Return the vtkCubeAxesActor2D object used to render the axes.
        This is to allow simple manipulation of the axes such as
        get_axes().SetNumberOfLabels(5) or similar.
        """
        return self.vtk_axes

    def setup_grid(self):
        """Create the vtkCellArray instance that represents the
        triangles. Subclasses are expected to override this function
        to read from their source as appropriate. The vtkCellArray should
        be stored to self.vtk_cells.
        """
        pass

    def render_quantity_height(self, quantityName, zScale=1.0, offset=0.0, dynamic=True):
        """Instruct the visualiser to render a quantity using the
        value at a point as its height.  The value at each point is
        multiplied by z_scale and is added to offset, and if
        dynamic=False, the quantity is not recalculated on each
        update.
        """
        self.height_quantities.append(quantityName)
        self.height_zScales[quantityName] = zScale
        self.height_offset[quantityName] = offset
        self.height_dynamic[quantityName] = dynamic

    def update_height_quantity(self, quantityName, dynamic=True):
        """Create a vtkPolyData object and store it in
        self.vtk_polyData[quantityName]. Subclasses are expected to override this
        function.
        """
        pass

    def get_3d_bounds(self):
        """Get the minimum and maximum bounds for the x, y and z directions.
        Return as a list of double in the order (xmin, xmax, ymin, ymax, zmin, zmax),
        suitable for passing to vtkCubeAxesActor2D::SetRanges(). Subclasses are expected
        to override this function.
        """
        pass

    def draw_height_quantity(self, quantityName):
        """Use the vtkPolyData and prepare/update the rest of the VTK
        rendering pipeline.
        """
        if self.vtk_mappers.has_key(quantityName):
            mapper = self.vtk_mappers[quantityName]
        else:
            mapper = self.vtk_mappers[quantityName] = vtkPolyDataMapper()
        mapper.SetInput(self.vtk_polyData[quantityName])
        mapper.Update()

        if not self.vtk_actors.has_key(quantityName):
            actor = self.vtk_actors[quantityName] = vtkActor()
            actor.SetMapper(mapper)
            self.vtk_renderer.AddActor(actor)

        if self.colours_height.has_key(quantityName):
            colour = self.colours_height[quantityName]
            if type(colour) == TupleType:
                if type(colour[0]) == FunctionType:
                    # It's a function, so take colour[1] as the
                    # lower bound on the scalar range and
                    # colour[2] as the upper bound on the scalar
                    # range.
                    scalars = vtkFloatArray()
                    map(scalars.InsertNextValue, colour[0](self.build_quantity_dict()))
                    self.vtk_polyData[quantityName].GetPointData().SetScalars(scalars)
                    mapper.SetScalarRange(colour[1:])
                    mapper.Update()
                else:
                    # It's a 3-tuple representing an RGB value.
                    actor.GetProperty().SetColor(colour)
            else:
                actor.GetProperty().SetColor(0.5, 0.5, 0.5)
        else:
            actor.GetProperty().SetColor(0.5, 0.5, 0.5)

    # --- Colour Coding --- #

    def build_quantity_dict(self):
        """Build and return a dictionary mapping quantity name->Numeric array of vertex
        values for that quantity. Subclasses are expected to override
        this function."""
        pass

    def colour_height_quantity(self, quantityName, colour=(0.5, 0.5, 0.5)):
        """Add colouring to a height based quantity.

        The colour parameter can be one of the following:
        - a 3-tuple of values in [0,1] to specify R, G, B values
        - a 3-tuple of values:
          - a function that takes a dictionary mapping quantity name->Numeric array of vertex values.
            This function returns a list of vertex values to be used in the colour coding.
          - a float for the lower bound on the colouring
          - a float for the upper bound on the colouring
        """
        self.colours_height[quantityName] = colour

    # --- Overlaid Polygons --- #

    def overlay_polygon(self, coords, height=0.0, colour=(1.0, 0.0, 0.0)):
        """Add a polygon to the output of the visualiser.

        coords is a list of 2-tuples representing x and y coordinates.
        These are triangulated by vtkDelaunay2D.

        height is the z-value given to all points.

        colour is the colour of the polygon, as a 3-tuple representing
        r, g, b values between 0 and 1."""
        points = vtkPoints()
        for coord in coords:
            points.InsertNextPoint(coord[0], coord[1], height)
        profile = vtkPolyData()
        profile.SetPoints(points)
        delny = vtkDelaunay2D()
        delny.SetInput(profile)
        mesh = vtkPolyDataMapper()
        mesh.SetInput(delny.GetOutput())
        actor = vtkActor()
        actor.SetMapper(mesh)
        actor.GetProperty().SetColor(colour)
        self.vtk_renderer.AddActor(actor)
        
    # --- Vector Fields --- #

    # --- GUI Setup --- #

    def setup_gui(self):
        self.tk_root = Tk()
        self.tk_root.title("Visualisation")
        self.tk_root.after(100, self.redraw)
        self.tk_root.bind("<Destroy>", self.destroyed)
        self.tk_root.grid_rowconfigure(0, weight=1)
        self.tk_root.grid_columnconfigure(0, weight=1)

        self.tk_renderWidget = vtkTkRenderWidget(self.tk_root, width=400, height=400)
        self.tk_renderWidget.grid(row=0, column=0, sticky=N+S+E+W)
        self.tk_controlFrame = Frame(self.tk_root)
        self.tk_controlFrame.grid(row=1, column=0, sticky=E+W)
        self.tk_controlFrame.grid_rowconfigure(0, weight=1)
        self.tk_controlFrame.grid_columnconfigure(0, weight=1)
        
        self.tk_quit = Button(self.tk_controlFrame, text="Quit", command=self.shutdown)
        self.tk_quit.grid(row=0, column=0, sticky=E+W)
        self.tk_renderWidget.GetRenderWindow().AddRenderer(self.vtk_renderer)

    # --- GUI Events --- #

    def destroyed(self, event):
        if event.widget == self.tk_root:
            self.shutdown()

    def redraw(self):
        self.tk_renderWidget.GetRenderWindow().Render()
        self.tk_root.update_idletasks()
        self.tk_root.after(100, self.redraw)

    def shutdown(self):
        self.tk_root.withdraw()
        self.tk_root.destroy()
