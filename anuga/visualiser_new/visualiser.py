
from os.path import splitext
from tkinter import Tk, Button, Frame, Label, Scale
from tkinter import N, E, S, W, HORIZONTAL, VERTICAL
from vtk import vtkCellArray, vtkRenderer, vtkWindowToImageFilter
from vtk import vtkJPEGWriter, vtkPNGWriter, vtkPNMWriter
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget

class Visualiser(object):
    '''
    Generic Offline Visualiser. Subclasses need to be used that specify
    how to handle a data source.
    '''
    
    ### Public functions ###

    def __init__(self,
                 title="Visualisation",
                 width=400,
                 height=400,
                 recording=False,
                 recordPattern=None,
                 paused=False,
                 source=None):
        '''
        Constructor.

        Params:
        title: string - Title of the visualisation window
        width: int - Width of the visualisation window
        height: int - Height of the visualisation window
        recording: boolean - Start with recording enabled?
        recordPattern: string - Pattern for recorded images,
          e.g., cylinders%05g.png
        paused: boolean - Start with playback paused?
        source:- The data source to read.
          What is required here varies by visualiser.
        '''

        # Visualisation options
        self.vis_features = []
        self.vis_frame = 0
        self.vis_frameStep = 1
        self.vis_jumping = False
        self.vis_recording = recording
        self.vis_recordPattern = recordPattern
        self.vis_paused = paused
        self.vis_source = source

        # VTK structures
        self.vtk_cells = vtkCellArray()
        self.vtk_renderer = vtkRenderer()

        # Tk structures
        self.tk_root = Tk()
        self.tk_root.title(title)
        self.tk_root.grid_rowconfigure(0, weight=1)
        self.tk_root.grid_columnconfigure(0, weight=3)
        self.tk_root.bind('<Destroy>', self.destroyed)
        if not self.vis_paused: self.tk_root.after(100, self.animate)            

        self.tk_renderWidget = vtkTkRenderWidget(self.tk_root,
                                                 width=width,
                                                 height=height)
        self.tk_renderWidget.grid(row=0, column=0, sticky=N+S+E+W)
        self.tk_renderWidget.GetRenderWindow().AddRenderer(self.vtk_renderer)

        self.tk_featureFrame = Frame(self.tk_root)
        self.tk_featureFrame.grid(row=0, column=1, rowspan=2)
        Label(self.tk_featureFrame, text='Features:').grid(row=0, column=0)

        self.tk_controlFrame = Frame(self.tk_root)
        self.tk_controlFrame.grid(row=1, column=0)

        self.tk_quit = Button(self.tk_controlFrame, text="Quit",
                              command=self.shutdown)
        self.tk_quit.grid(row=0, column=0, columnspan=2)

        def pause():
            if self.vis_paused:
                self.tk_pause.config(text='Pause')
                self.tk_root.after(100, self.animate)
            else:
                self.tk_pause.config(text='Resume')
            self.vis_paused ^= True
        self.tk_pause = Button(self.tk_controlFrame, text="Pause", command=pause)
        self.tk_pause.grid(row=0, column=2, columnspan=2)

        if self.vis_recordPattern is not None:
            def record():
                if self.vis_recording:
                    self.tk_record.config(text='Start Recording')
                else:
                    self.tk_record.config(text='Stop Recording')
                self.vis_recording ^= True
            self.tk_record = Button(self.tk_controlFrame, text="Start Recording", command=record)
            self.tk_record.grid(row=0, column=4, columnspan=2)
            if self.vis_recording:
                self.tk_record.config(text="Stop Recording")

        def make_seek_button(label, column, frame):
            def jump():
                self.jumpTo(frame)
            b = Button(self.tk_controlFrame,
                       text=label,
                       command=jump)
            b.grid(row=1, column=column, sticky=W+E)
            return b
        self.tk_seek_start = make_seek_button("|<", 0, 0)
        self.tk_seek_back10 = make_seek_button("<<", 1,
                                               lambda: self.vis_frame - 10)
        self.tk_seek_back1 = make_seek_button("<", 2,
                                              lambda: self.vis_frame - 1)
        self.tk_seek_forward1 = make_seek_button(">", 3,
                                                 lambda: self.vis_frame + 1)
        self.tk_seek_forward10 = make_seek_button(">>", 4,
                                                  lambda: self.vis_frame + 10)
        self.tk_seek_end = make_seek_button(">|", 5,
                                            self.getMaxFrameNumber)

        Label(self.tk_controlFrame, text='Frame').grid(row=2, column=0,
                                                       sticky=W+E)
        def changeFrame(frame):
            if not self.vis_jumping:
                self.vis_jumping = True
                self.jumpTo(self.tk_frame.get())
                self.vis_jumping = False
        self.tk_frame = Scale(self.tk_controlFrame, command=changeFrame,
                              from_=0, to=0, orient=HORIZONTAL)
        self.tk_frame.grid(row=2, column=1, columnspan=2, sticky=W+E)

        Label(self.tk_controlFrame, text='Step').grid(row=2, column=3,
                                                     sticky=W+E)
        def changeFrameStep(step):
            self.vis_frameStep = int(step)
        self.tk_frameStep = Scale(self.tk_controlFrame, command=changeFrameStep,
                                  from_=1, to=1, orient=HORIZONTAL)
        self.tk_frameStep.grid(row=2, column=4, columnspan=2, sticky=W+E)

        self.setupGrid()

    def add_feature(self, feature):
        '''Add a feature to this visualiser'''
        self.vis_features.append(feature)
        feature.button(self.tk_featureFrame).grid(row=len(self.vis_features),
                                                  column=0, sticky=W+E)
        feature.visualiser = self

    def run(self):
        '''Start the visualiser'''
        self.redraw()
        self.tk_root.mainloop()

    ### Private funcitions ###

    def resetSliders(self):
        '''
        Recalculate the upper bound on the frame and frameStep
        sliders.
        '''
        maxFrame = self.getMaxFrameNumber()
        self.tk_frame.config(to=maxFrame)
        self.tk_frameStep.config(to=maxFrame)

    def jumpTo(self, frame):
        '''
        Jump to a given frame. If frame is a function, jump to the
        return value of frame().
        '''
        oldFrame = self.vis_frame
        if hasattr(frame, '__call__'):
            self.vis_frame = frame()
        else:
            self.vis_frame = frame

        maxFrame = self.getMaxFrameNumber()

        if self.vis_frame < 0:
            self.vis_frame = 0
        elif self.vis_frame > maxFrame:
            self.vis_frame = maxFrame
            self.vis_paused = True
            self.tk_pause.config(text='Resume')

        self.tk_frame.set(self.vis_frame)
        self.redraw(oldFrame != self.vis_frame)

    def redraw(self, update=False):
        self.resetSliders()
        for feature in [ f for f in self.vis_features if not f.dynamic ]:
            feature.draw(self.vtk_renderer)
        if update:
            for feature in [ f for f in self.vis_features if f.dynamic]:
                f.redraw(self.vtk_renderer)
        self.tk_renderWidget.GetRenderWindow().Render()
        self.tk_root.update_idletasks()

    ### Gui events ###

    def destroyed(self, event):
        if event.widget == self.tk_root:
            self.shutdown()

    def shutdown(self):
        self.tk_root.withdraw()
        self.tk_root.destroy()

    def animate(self):
        if not self.vis_paused:
            self.jumpTo(self.vis_frame + self.vis_frameStep)
            if self.vis_recording and self.vis_recordPattern is not None:
                self.save_image()
            self.tk_root.after(100, self.animate)

    def save_image(self):
        extmap = {'.jpg' : vtkJPEGWriter,
                  '.jpeg' : vtkJPEGWriter,
                  '.png' : vtkPNGWriter,
                  '.pnm' : vtkPNMWriter}
        _, ext = splitext(self.vis_recordPattern)
        try: writer = extmap[ext.lower()]()
        except KeyError:
            print('ERROR: Can\'t handle %s extension. Recording disabled.' % ext)
            self.vis_recordPattern = None
            return

        win = self.vtk_renderer.GetRenderWindow()
        w2i = vtkWindowToImageFilter()
        w2i.SetInput(win)
        w2i.Update()
        writer.SetInput(w2i.GetOutput())
        writer.SetFileName(self.vis_recordPattern % self.vis_frame)
        win.Render()
        writer.Write()

    ### Things subclasses need to override ###

    def setupGrid(self):
        '''
        Populate the vtkCellArray instance at
        self.vtk_cells. Subclasses are required to override this
        function to read from their source as appropriate.
        '''
        raise NotImplementedError('Subclass needs to override Visualiser::setupGrid!')

    def getMaxFrameNumber(self):
        '''
        Return the maximum frame number. This will need to be defined
        by a subclass.
        '''
        raise NotImplementedError('Subclass needs to override Visualiser::getMaxFrameNumber!')

    def getQuantityPoints(self, quantityName, dynamic=True, frameNumber=0):
        '''
        Return the points of a quantity at a given frame as a list
        [float]. Subclasses need to override this.
        '''
        raise NotImplementedError('Subclass needs to override Visualiser::getQuantityPoints!')

    def getQuantityDict(self):
        '''
        Return the values of all quantities at a given time as a
        dictionary. Sublclasses need to override this.
        '''
        raise NotImplementedError('Subclass needs to override Visualiser::getQuantityDict!')
