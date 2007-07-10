from enthought.traits.api import Int, String
from Scientific.IO.NetCDF import NetCDFFile
from Tkinter import Label, Scale, W, E, HORIZONTAL
from visualiser import Visualiser

class SWWVisualiser(Visualiser):
    '''A Visualiser that views SWW files. The source for this
    visualiser is a string containing the sww file name.
    '''
    source = String
    frameDelay = Int(100)
    frameStep = Int(1)

    def __init__(self, *args, **kwargs):
        Visualiser.__init__(self, *args, **kwargs)

        fin = NetCDFFile(self.source, 'r')
        self.frame = 0
        self.maxFrame = fin.variables['time'].shape[0] - 1
        self.paused = False
        self.heightFeatureCache = []
        for i in range(self.maxFrame + 1): # Zero indexed
            self.heightFeatureCache.append({})

    def setup_grid(self):
        fin = NetCDFFile(self.source, 'r')
        N_tri = fin.variables['volumes'].shape[0]
        for v in range(N_tri):
            self.vtk_cells.InsertNextCell(3)
            for i in range(3):
                self.vtk_cells.InsertCellPoint(fin.variables['volumes'][v][i])
        fin.close()
        
    def setup_gui(self):
        Visualiser.setup_gui(self)

        Label(self.tk_customControlFrame, text='Frame').grid(row=0, column=0, sticky=W+E)
        self.tk_frame = Scale(self.tk_customControlFrame, from_=0, to=self.maxFrame, orient=HORIZONTAL)
        self.tk_frame.grid(row=0, column=1, sticky=W+E)
        Label(self.tk_customControlFrame, text='Step').grid(row=0, column=2, sticky=W+E)
        self.tk_frameStep = Scale(self.tk_customControlFrame, from_=0, to=self.maxFrame, orient=HORIZONTAL)
        self.tk_frameStep.grid(row=0, column=3, sticky=W+E)
        
