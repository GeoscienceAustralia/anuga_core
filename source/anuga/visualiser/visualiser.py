from threading import Event, Thread
from Tkinter import Tk, Button, N, E, S, W
from types import FunctionType, TupleType
from vtk import vtkActor, vtkFloatArray, vtkPolyDataMapper, vtkRenderer
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
        self.vtk_mappers = {}
        self.vtk_polyData = {}

    def run(self):
        self.setup_gui()
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

    # --- Height Based Rendering --- #

    def setup_grid(self):
        """Create the vtkCellArray instance that represents the
        triangles. Subclasses are expected to override this function
        to read from their source as appropriate.
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
        self.vtk_polyData[q]. Subclasses are expected to override this
        function.
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
        """Build a dictionary mapping quantity name->list of vertex
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
            
    # --- Vector Fields --- #

    # --- GUI Setup --- #

    def setup_gui(self):
        self.tk_root = Tk()
        self.tk_root.title("Visualisation")
        self.tk_root.after(100, self.redraw)
        self.tk_root.bind("<Destroy>", self.destroyed)

        self.tk_renderWidget = vtkTkRenderWidget(self.tk_root, width=400, height=400)
        self.tk_renderWidget.grid(row=0, column=0, sticky=N+S+E+W)
        self.tk_quit = Button(self.tk_root, text="Quit", command=self.shutdown)
        self.tk_quit.grid(row=2, column=0, sticky=E+W)
        self.vtk_renderer = vtkRenderer()
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
        self.tk_root.quit()
