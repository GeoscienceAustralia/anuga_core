
import numpy as num
from anuga.file.netcdf import NetCDFFile
from tkinter import Button, E, Tk, W, Label, StringVar, Scale, HORIZONTAL
from .visualiser import Visualiser
from vtk import vtkCellArray, vtkPoints, vtkPolyData

class OfflineVisualiser(Visualiser):
    """A VTK-powered offline visualiser which runs in its own thread.
    In addition to the functions provided by the standard visualiser,
    the following additional functions are provided:

    precache_height_quantities() - Precache all the vtkpoints
    structures for any dynamic height based quantities to render.
    """
    def __init__(self, source, frameDelay=100, frameStep=1):
        """The source parameter is assumed to be a NetCDF sww file.
        The frameDelay parameter is the number of milliseconds waited between frames.
        """
        Visualiser.__init__(self, source)

        self.frameNumber = 0
        fin = NetCDFFile(self.source, 'r')
        self.maxFrameNumber = fin.variables['time'].shape[0] - 1
        fin.close()
        
        #self.frameNumberTkVariable = StringVar()
        #self.frameNumberTkVariable.set('Frame - %05g'%self.framNumber)

        self.frameDelay = frameDelay

        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None

        self.frameStep= frameStep 

        self.vtk_heightQuantityCache = []
        for i in range(self.maxFrameNumber + 1): # maxFrameNumber is zero indexed.
            self.vtk_heightQuantityCache.append({})

        self.paused = False
        self.movie = False
        
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
        polydata = self.vtk_polyData[quantityName] = vtkPolyData()
        if dynamic is True:
            #print ' - Frame',self.frameNumber,'of',self.maxFrameNumber
            if quantityName not in self.vtk_heightQuantityCache[self.frameNumber]:
                self.vtk_heightQuantityCache[self.frameNumber][quantityName]\
                    = self.read_height_quantity(quantityName, True, self.frameNumber);
            polydata.SetPoints(self.vtk_heightQuantityCache[self.frameNumber][quantityName])
        else:
            polydata.SetPoints(self.read_height_quantity(quantityName, False))
        polydata.SetPolys(self.vtk_cells)

    def get_3d_bounds(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax]
            
    def read_height_quantity(self, quantityName, dynamic=True, frameNumber=0):
        """Read in a height based quantity from the NetCDF source file
        and return a vtkPoints object. frameNumber is ignored if
        dynamic is false."""
        fin = NetCDFFile(self.source, 'r')
        points = vtkPoints()
        if dynamic is True:
            N_vert = fin.variables[quantityName].shape[1]
        else:
            N_vert = len(fin.variables[quantityName])
        x = num.ravel(num.array(fin.variables['x'], float))
        y = num.ravel(num.array(fin.variables['y'], float))
        if dynamic is True:
            q = num.array(fin.variables[quantityName][frameNumber], float)
        else:
            q = num.ravel(num.array(fin.variables[quantityName], float))

        q *= self.height_zScales[quantityName]
        q += self.height_offset[quantityName]

        for v in range(N_vert):
            points.InsertNextPoint(x[v], y[v], q[v])
            if self.xmin is None or self.xmin > x[v]:
                self.xmin = x[v]
            if self.xmax is None or self.xmax < x[v]:
                self.xmax = x[v]
            if self.ymin is None or self.ymin > y[v]:
                self.ymin = y[v]
            if self.ymax is None or self.ymax < y[v]:
                self.ymax = y[v]
            if self.zmin is None or self.zmin > q[v]:
                self.zmin = q[v]
            if self.zmax is None or self.zmax < q[v]:
                self.zmax = q[v]
        fin.close()
        return points

    def precache_height_quantities(self):
        """Precache any height-based quantities. Call before rendering
        beigns."""
        for q in self.height_quantities:
            if self.height_dynamic[q] is True:
                print('Precaching %s' % q)
                for i in range(self.maxFrameNumber + 1): # maxFrameNumber is zero-indexed
                    print(' - Frame %d of %d' % (i, self.maxFrameNumber))
                    self.vtk_heightQuantityCache[i][q]\
                        = self.read_height_quantity(q, True, i)

    def build_quantity_dict(self):
        quantities = {}
        fin = NetCDFFile(self.source, 'r')
        for q in [n for n in list(fin.variables.keys()) if n != 'x' and n != 'y' and n != 'z' and n != 'time' and n != 'volumes']:
            if len(fin.variables[q].shape) == 1: # Not a time-varying quantity
                quantities[q] = num.ravel(num.array(fin.variables[q], float))
            else: # Time-varying, get the current timestep data
                quantities[q] = num.array(fin.variables[q][self.frameNumber], float)
        fin.close()
        return quantities

    def setup_gui(self):
        Visualiser.setup_gui(self)
        self.tk_quit.grid(row=0, column=0, sticky=W+E)
        self.tk_movie_toggle = Button(self.tk_controlFrame, text="Movie off", command=self.movie_toggle)
        self.tk_movie_toggle.grid(row=0, column=6,  sticky=W+E)
        
                
        self.tk_restart = Button(self.tk_controlFrame, text="<<<", command=self.restart, width=5)
        self.tk_restart.grid(row=1, column=0, sticky=W+E)
        self.tk_back10 = Button(self.tk_controlFrame, text="<<", command=self.back10, width=5)
        self.tk_back10.grid(row=1, column=1, sticky=W+E)
        self.tk_back = Button(self.tk_controlFrame, text="<", command=self.back, width=5)
        self.tk_back.grid(row=1, column=2, sticky=W+E)
        self.tk_pauseResume = Button(self.tk_controlFrame, text="Pause", command=self.pauseResume, width=15)
        self.tk_pauseResume.grid(row=1, column=3, sticky=W+E)
        self.tk_forward = Button(self.tk_controlFrame, text=">", command=self.forward, width=5)
        self.tk_forward.grid(row=1, column=4, sticky=W+E)
        self.tk_forward10 = Button(self.tk_controlFrame, text=">>", command=self.forward10, width=5)
        self.tk_forward10.grid(row=1, column=5, sticky=W+E)
        self.tk_forwardEnd = Button(self.tk_controlFrame, text=">>>", command=self.forwardEnd, width=5)
        self.tk_forwardEnd.grid(row=1, column=6, sticky=W+E)
        

        self.tk_frameNumber = Label(self.tk_controlFrame, text='Frame')
        self.tk_frameNumber.grid(row=2, column=0, sticky=W+E)
        self.tk_gotoFrame = Scale(self.tk_controlFrame, from_=0, to=self.maxFrameNumber, orient=HORIZONTAL)
        self.tk_gotoFrame.grid(row=2, column=1, columnspan=2, sticky=W+E)
        self.tk_stepLabel = Label(self.tk_controlFrame, text='Step')
        self.tk_stepLabel.grid(row=2, column=4, sticky=W+E)        
        self.tk_frameStep = Scale(self.tk_controlFrame, from_=0, to=self.maxFrameNumber, orient=HORIZONTAL)
        self.tk_frameStep.grid(row=2, column=5, columnspan=2, sticky=W+E)
        
        # Make the buttons stretch to fill all available space
        for i in range(7):
            self.tk_controlFrame.grid_columnconfigure(i, weight=1)

    def run(self):
        self.alter_tkroot(Tk.after, (self.frameDelay, self.animateForward))
        Visualiser.run(self)

    def restart(self):
        self.frameNumber = 0
        self.redraw_quantities()
        self.update_labels()
        self.pause()
        
        if self.movie:
            self.save_image()
 
    def forwardEnd(self):
        self.frameNumber = self.maxFrameNumber
        self.redraw_quantities()
        self.update_labels()
        self.pause()
                
    def movie_toggle(self):
        if self.movie == True:
            self.movie = False
            self.tk_movie_toggle.config(text='Movie off')
        else:
            self.movie = True
            self.tk_movie_toggle.config(text='Movie on ')
            
            
        
        
    def save_image(self):
        
        from vtk import vtkJPEGWriter, vtkJPEGWriter, vtkPNGWriter
        from vtk import vtkPNMWriter, vtkWindowToImageFilter
        from os import path
         
        sourcebase, _ = path.splitext(self.source)
        fname = sourcebase+'%05g.png' % self.frameNumber
        #print fname
        
        extmap = {'.jpg' : vtkJPEGWriter,
                  '.jpeg' : vtkJPEGWriter,
                  '.png' : vtkPNGWriter,
                  '.pnm' : vtkPNMWriter,
                  }
        basename, ext = path.splitext(fname)
        try: Writer = extmap[ext.lower()]
        except KeyError:
            error_msg("Don't know how to handle %s files" % ext, parent=self)
            return
    
        renWin = self.vtk_renderer.GetRenderWindow()
        w2i = vtkWindowToImageFilter()
        writer = Writer()
        w2i.SetInput(renWin)
        w2i.Update()
        writer.SetInput(w2i.GetOutput())
        writer.SetFileName(fname)
        renWin.Render()
        writer.Write()        
    
    def back10(self):
        if self.frameNumber - 10 >= 0:
            self.frameNumber -= 10
        else:
            self.frameNumber = 0
        self.redraw_quantities()
        self.update_labels()
        self.pause()

    def back(self):
        if self.frameNumber > 0:
            self.frameNumber -= 1
            self.redraw_quantities()
            self.update_labels()
            self.pause()

    def pauseResume(self):
        if self.paused is True:
            self.resume()
        else:
            self.pause()

    def pause(self):
        self.paused = True
        self.tk_pauseResume.config(text="Resume")

    def resume(self):
        self.paused = False
        self.tk_pauseResume.config(text="Pause")
        self.frameNumber = self.tk_gotoFrame.get()
        self.frameStep = self.tk_frameStep.get()
        self.tk_root.after(self.frameDelay, self.animateForward)

    def forward(self):
        if self.frameNumber < self.maxFrameNumber:
            self.frameNumber += 1
            self.redraw_quantities()
            self.update_labels()
            self.pause()
            
    def forward_step(self):
        if self.frameNumber + self.frameStep <= self.maxFrameNumber:
            self.frameNumber += self.frameStep
            self.redraw_quantities()
            self.update_labels()
        else:
            self.frameNumber = self.maxFrameNumber            
            self.redraw_quantities()
            self.update_labels()    
            self.pause()
         
        if self.movie:
             self.save_image()
                

    def forward10(self):
        if self.frameNumber + 10 <= self.maxFrameNumber:
            self.frameNumber += 10
        else:
            self.frameNumber = self.maxFrameNumber
        self.redraw_quantities()
        self.update_labels()
        self.pause()

    def animateForward(self):
        if self.paused is not True:
            self.forward_step()
            self.tk_root.after(self.frameDelay, self.animateForward)
            
    def update_labels(self): 
        #self.tk_frameNumber.config(text='%05g of %05g'%(self.frameNumber,self.maxFrameNumber))
        self.tk_gotoFrame.set(self.frameNumber)
        self.tk_frameStep.set(self.frameStep)
               
    def shutdown(self):
        #self.pause()
        self.tk_root.withdraw()
        self.tk_root.destroy()
        #Visualiser.shutdown(self)
