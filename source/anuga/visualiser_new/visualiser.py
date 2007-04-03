from axes import Axes
from enthought.traits.api import false, HasTraits, Instance, List
from features import Feature
from threading import Thread
from Tkinter import Tk, Button, Frame, N, E, S, W
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget

class Visualiser(HasTraits, Thread):
    '''Generic visualiser class'''
    source = Any
    features = List(Feature)
    axes = Instance(Axes)
    recording = false

    def __init__(self, *args, **kwargs):
        HasTraits.__init__(self, *args, **kwargs)
        Thread.__init__(self)

    def run(self):
        self.setup_gui()
        self.tk_root.mainloop()

    def setup_gui(self):
        self.tk_root = Tk()
        self.tk_root.title("Visualisation")
        
        self.tk_renderWidget = vtkTkRenderWidget(self.tk_root, width=400, height=400)
        self.tk_renderWidget.grid(row=0, column=0, sticky=N+S+E+W)
        self.tk_root.grid_rowconfigure(0, weight=1)
        self.tk_root.grid_columnconfigure(0, weight=3)

        self.tk_featureFrame = Frame(self.tk_root)
        self.tk_featureFrame.grid(row=0, column=1)

        self.tk_controlFrame = Frame(self.tk_root)
        self.tk_controlFrame.grid(row=1, column=0, columnspan=2)
        self.tk_quit = Button(self.tk_controlFrame, text="Quit", command=self.shutdown)
        self.tk_quit.grid(row=0, column=0)
        self.tk_pause = Button(self.tk_controlFrame, text="Pause")
        self.tk_pause.grid(row=0, column=1)
        #TODO: Add record binding here
        self.tk_record = Button(self.tk_controlFrame, text="Record")
        self.tk_record.grid(row=0, column=2)
        if self.recording == True:
            self.tk_record.config(text="Stop Recording")

        self.tk_customControlFrame = Frame(self.tk_root)
        self.tk_customControlFrame.grid(row=2, column=0, columnspan=2)

        self.tk_root.after(100, self.redraw)
        self.tk_root.bind("<Destroy>", self.destroyed)

    def destroyed(self, event):
        if event.widget == self.tk_root:
            self.shutdown()

    def shutdown(self):
        self.tk_root.withdraw()
        self.tk_root.destroy()

    def redraw(self):
        self.tk_renderWidget.GetRenderWindow().Render()
        self.tk_root.update_idletasks()
        self.tk_root.after(100, self.redraw)
