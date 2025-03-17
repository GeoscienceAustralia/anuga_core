#!/usr/bin/env python3

import  Pmw, math, time, string, marshal

try:
    from . import AppShell
    from .toolbarbutton import ToolBarButton
    from . import visualmesh
    from . import mesh
    from .mesh import SEG_COLOUR
except:
    from anuga.pmesh import AppShell
    from anuga.pmesh.toolbarbutton import ToolBarButton
    import anuga.pmesh.visualmesh as visualmesh
    import anuga.pmesh.mesh as mesh
    from anuga.pmesh.mesh import SEG_COLOUR


import tkinter.filedialog
from   tkinter.simpledialog import Dialog


from tkinter import  FALSE,TRUE, Frame,X, LEFT,YES,BOTH,ALL,Widget,CURRENT, \
     Label,W, Entry, E, StringVar, END, Checkbutton, Radiobutton, IntVar, \
     DISABLED, NORMAL
#from cursornames import TLC,TRC, BLC, BRC, TS, RS, LS, BS
from tkinter.messagebox import showerror, _show, QUESTION,YESNOCANCEL, askyesno
import types

import os, sys
import profile
import anuga.load_mesh.loadASCII
from anuga.alpha_shape.alpha_shape import AlphaError
import anuga.utilities.log as log

# CONSTANTS
VERT_SELECT_ADDING_SEG_COLOR = 'orange'
SELECT_COLOR = 'red'
TRIANGLE_COLOUR = 'green'
APPLICATION_NAME = 'Pmesh'

#for alpha shapes
NO_SELECTION = 0
AUTO = 1
SET_ALPHA = 2

# pmesh home directory
from anuga.utilities.system_tools import get_pathname_from_package
HOME_DIR = get_pathname_from_package('anuga.pmesh')
#HOME_DIR = os.path.dirname(os.path.abspath(__file__))


class Draw(AppShell.AppShell):
    usecommandarea = 1
    appname        = APPLICATION_NAME
    frameWidth     = 840
    frameHeight    = 600




    def createButtons(self):
        """
        Add buttons to the bottom of the GUI
        """
        # self.buttonAdd('Postscript',
        #       helpMessage='Save current drawing (as PostScript)',
        #       statusMessage='',
        #       command=self.ipostscript)
        self.buttonAdd('Clear', helpMessage='Delete the mesh',
              statusMessage='', command=self.clearMesh)
        self.buttonAdd('Close', helpMessage='Close Screen',
              statusMessage='', command=self.close)

    def createBase(self):
        """
        Create the GUI framework.  Set up the GUI
        """
        self.toolbar = self.createcomponent('toolbar', (), None,
                  Frame, (self.interior(),), background="gray95")
        self.toolbar.pack(fill=X)

        self.scrolledcanvas =  self.createcomponent('ScrolledCanvas', (), None,
                                       Pmw.ScrolledCanvas, (self.interior(),)
                                       ,borderframe = 1
                                       ,labelpos = 'n'
                                       )
        self.scrolledcanvas.configure(hscrollmode = 'dynamic')
        self.scrolledcanvas.configure(vscrollmode = 'dynamic')
        self.scrolledcanvas.pack(side=LEFT, expand=YES, fill=BOTH)
        self.canvas = self.scrolledcanvas.component('canvas')
        self.canvas.configure( background="white" )

        self.canvas.pack(side=LEFT, expand=YES, fill=BOTH)

        Widget.bind(self.canvas, "<Button-1>", self.mouseDown)
        Widget.bind(self.canvas, "<Button3-ButtonRelease>", self.rightMouseUp)
        Widget.bind(self.canvas, "<Button2-ButtonRelease>",self.DeleteSelectedMeshObject)
        Widget.bind(self.canvas, "<Motion>", self.displayCoords)
        # "<Delete>" didn't work..
        #Widget.bind(self.canvas, "<Delete>", self.DeleteSelectedMeshObject)

        #self.root.bind("<KeyPress>", self.setRegular)
        #self.root.bind("<KeyRelease>", self.setRegular)

        self.scrolledcanvas.resizescrollregion()

#     def setRegular(self, event):
#         if event.type == '2' and event.keysym == 'Shift_L':
#             self.regular = TRUE
#         else:
#             self.regular = FALSE

    def createMenus(self):
        """
        Add menus to the top of the GUI
        """
        self.menuBar.deletemenuitems('File',0)
        self.menuBar.addmenuitem('File', 'command', 'New mesh',
                                 label='New', command=self.clearMesh)
        self.menuBar.addmenuitem('File', 'command', 'Open mesh',
                                 label='Open...', command=self.importFile)
        self.menuBar.addmenuitem('File', 'command', 'Save mesh',
                                 label='Save', command=self.saveDrawing)
        self.menuBar.addmenuitem('File', 'command', 'Save mesh',
                                 label='SaveAs...', command=self.saveAsDrawing)

        self.menuBar.addmenuitem('File', 'separator')
        self.menuBar.addmenuitem('File', 'command',
                                 'Add ungenerated file from arcGIS',
                                 label='Add ungenerated file...',
                                 command=self.ImportUngenerate)

        #self.menuBar.addmenuitem('File', 'command',
        #                         'Export ASCII obj',
        #                         label='Export ASCII obj',
        #                         command=self.exportObj)

        self.menuBar.addmenuitem('File', 'command',
                                 'Export ASCII segment outline',
                                 label='Export ASCII segment outline...',
                                 command=self.exportASCIIsegmentoutlinefile)

        self.menuBar.addmenuitem('File', 'command',
                                 'Export ASCII csv file',
                                 label='Export ASCII csv file...',
                                 command=self.exportPointsFile)

        self.menuBar.addmenuitem('File', 'command',
                                 'Export Postcript file',
                                 label='Export Postscript file...',
                                 command=self.ipostscript)

        # self.menuBar.addmenuitem('File', 'command',
        #                          'add Segments to connect all vertices'  ,
        #                          label='join vertices',
        #                          command=self.joinVertices)
        # self.menuBar.addmenuitem('File', 'command',
        #                          'add Segments to form alpha shape'  ,
        #                          label='Auto segment',
        #                          command=self.auto_segment)
        # self.menuBar.addmenuitem('File', 'command',
        #                    'modify the alpha boundary by applying filters',
        #                          label='filter alpha boundary',
        #                          command=self.auto_segmentFilter)
        #self.menuBar.addmenuitem('File', 'command', 'Normalise mesh',
        #                         label='Normalise mesh', command=self.normaliseMesh)
        #self.menuBar.addmenuitem('File', 'command', 'Normalise mesh for glutobj',
         #                        label='Normalise mesh for glutobj', command=self.normalise4ObjMesh)
        self.menuBar.addmenuitem('File', 'separator')
        self.menuBar.addmenuitem('File', 'command', '',
                                 label='Print geo reference', command=self.printGeoReference)
        self.menuBar.addmenuitem('File', 'separator')
        self.menuBar.addmenuitem('File', 'command', 'Exit program',
                                 label='Exit', command=self.quit)

    def createTools(self):
        """
        Add buttons to the top of the GUI
        """
        self.mouseDownFunc = {}
        self.modeClass = {}
        ToolBarButton(self, self.toolbar, 'sep', 'Separator.gif',
                      width=10, state='disabled',home_dir=HOME_DIR)
        for key, balloon, mouseDownFunc, Mode in [
            #('Pointer','Edit drawing eventually.  Right now this does nothing', self.drag, None)
            ('Add-Vertex',    'Add vertex', self.drawVertex, mesh.Vertex)
            ,('Segment', 'Join vertices to form a segment',self.selectSegmentPoint, mesh.Segment)
            ,('Add-Hole', 'Add hole',self.drawHole, mesh.Hole)
            ,('Add-Region', 'Add region',self.drawRegion, mesh.Region)
            ]:
            t = ToolBarButton(self, self.toolbar, key, '%s.gif' % key,
                          command=self.selectFunc, balloonhelp=balloon,
                               statushelp='', home_dir=HOME_DIR)
            t.cycle("DrawMode")
            if key == 'Add-Vertex': #FIXME- this is specified in line 1062 as well
                                 # self.selectFunc('pointer')
                self.curFunc  = self.drawVertex
                t.setInitialSunkenButton("DrawMode")
            self.modeClass[key] = Mode
            # for actions that occur when the mouse goes down
            self.mouseDownFunc[key] = mouseDownFunc

    def createZooms(self):
        """
        Add zoom buttons to the top of the GUI
        """
        ToolBarButton(self, self.toolbar, 'sep', 'Separator.gif', width=10,
                      state='disabled', home_dir=HOME_DIR)
        zoom = '0.5'
        ToolBarButton(self, self.toolbar, zoom, 'Zoom-Out.gif',
                      command=self.selectZoom,
                      balloonhelp='Zoom out',
                      statushelp='', home_dir=HOME_DIR)

        ToolBarButton(self, self.toolbar,'1.0', 'Zoom-Extents.gif',
                      command=self.ResizeToFitWrapper,
                      balloonhelp='Zooms to mesh size',
                      statushelp='', home_dir=HOME_DIR)
        zoom = '2'
        ToolBarButton(self, self.toolbar, zoom, 'Zoom-In.gif',
                      command=self.selectZoom,
                      balloonhelp='Zoom in',
                      statushelp='', home_dir=HOME_DIR)

    def createEdits(self):
        """
        Add Edit buttons to the top of the GUI
        """
        ToolBarButton(self, self.toolbar, 'sep', 'Separator.gif', width=10,
                      state='disabled', home_dir=HOME_DIR)
        for key, func, balloon in [
                ('Add-Vertex-Dialog', self.windowAddVertex, 'Add vertex dialog'),
                ('Delete', self.windowDelete, 'Delete selected object'),
                ('Edit', self.windowEdit, 'Edit selected object'),
                ('Tag-Segment', self.windowDefault, 'Set default tag value for new segments'),
                ('Join-Vertices', self.joinVerticesButton, 'Add Segments to connect all vertices'),
             #   ('autoSeg', self.auto_segmentButton, 'add Segments to form alpha shape'),
                ('Alpha', self.auto_segmentGiveAlphaButton, 'Add Segments to form alpha shape, specify alpha'),
                ('Mesh', self.windowMeshGen, 'Generate Mesh')]:
            ToolBarButton(self, self.toolbar, key, '%s.gif' % key,
                          command=func, balloonhelp=balloon,
                               statushelp='', home_dir=HOME_DIR)


    def createVisualiseIcons(self):
        """
        Add Edit buttons to the top of the GUI
        """
        ToolBarButton(self, self.toolbar, 'sep', 'Separator.gif', width=10,
                      state='disabled', home_dir=HOME_DIR)
        for key, func, balloon in [
                ('Mesh-Visibility-Toggle', self.visualise, 'Toggle visibility of mesh')]:
            ToolBarButton(self, self.toolbar, key, '%s.gif' %key,
                          command=func, balloonhelp=balloon,
                               statushelp='', home_dir=HOME_DIR)


    def clearSelection(self,parent):
    #FIXME looks like self.clearSelections - change name (Peter)
        """
        """
        self.canvas.delete(ALL)
        self.selSet = self.mesh.clearSelection()
        self.visualiseMesh(self.mesh)

    def visualise(self,parent):
        self.canvas.delete(ALL)
        if self.Visualise:
            self.Visualise = False
            self.visualiseMesh(self.mesh)
        else:
            self.Visualise = True
            self.visualiseMesh(self.mesh)

    def createMesh(self):
        """
        Build the data structures for storing the mesh objects
        """
        self.mesh = mesh.Mesh()

        self.Vertices = visualmesh.vPoints(mesh.Vertex)
        self.Segments = visualmesh.vSegments(mesh.Segment)
        self.Holes = visualmesh.vPoints(mesh.Hole)
        self.Regions = visualmesh.vRegions(mesh.Region)
        self.UserMesh = visualmesh.vMesh([self.Vertices,self.Segments,self.Holes,self.Regions])


    def deleteMesh(self):
        """
        Delete the data structures for storing the mesh objects
        """
        self.mesh = None
        self.Vertices = None
        self.Segments = None

    def addCylinders(self):
        """
        Automatically add some verts and segs to the mesh.Used in Debugging

        center and radius

        """
        from anuga.coordinate_transforms.geo_reference import Geo_reference, \
             DEFAULT_ZONE
        offset_x = 30
        offset_y = 40

        x_origin = 10-offset_x
        y_origin = 20-offset_y
        r = 10
        pi = math.pi
        num_of_cuts = 100
        cuts = []
        factor = 2* math.pi/num_of_cuts
        for cut in range(num_of_cuts):
             cuts.append(cut*factor)

        for radius in cuts:
            x = x_origin + r * math.cos(radius)
            y = y_origin + r * math.sin(radius)
            v = self.drawVertex(x,y,None)
            if not radius == 0.0:   # FIXME
                self.drawSegment(v,v_old)
            else:
                v_first = v
            v_old = v
        self.drawSegment(v,v_first)
        region = self.drawRegion(x_origin, y_origin, 0)
        region.setTag("setheight5")


        x_origin = 30-offset_x
        y_origin = 30-offset_y
        r = 5
        pi = math.pi
        num_of_cuts = 100
        cuts = []
        factor = 2* math.pi/num_of_cuts
        for cut in range(num_of_cuts):
             cuts.append(cut*factor)

        for radius in cuts:
            x = x_origin + r * math.cos(radius)
            y = y_origin + r * math.sin(radius)
            v = self.drawVertex(x,y,None)
            if not radius == 0.0:   # FIXME
                self.drawSegment(v,v_old)
            else:
                v_first = v
            v_old = v
        self.drawSegment(v,v_first)
        region = self.drawRegion(x_origin, y_origin, 0)
        region.setTag("setheight10")
        self.mesh.geo_reference = Geo_reference(zone=DEFAULT_ZONE,
                                                xllcorner=offset_x,
                                                yllcorner=offset_y)

        #Since the new vertex may be off screen
        self.scrolledcanvas.resizescrollregion()

    def selectFunc(self, tag):
        """
        Change the current mode class
        When changing from one mode to another
        """
        self.mouseDownCurFunc = self.mouseDownFunc[tag]
        self.curModeClass = self.modeClass[tag]
        self.clearSelections()
        self.canvas.config(cursor='arrow')
#        I can make it arrow, but it will change back to pointer, after
#        adding an object
#        if self.curFunc == self.func['pointer']:
#            self.canvas.config(cursor='arrow')
#         else:
#             self.canvas.config(cursor='crosshair')

    def clearSelections(self):
        """
        deselect objects that have been selected
        """
        if self.selMeshObject:
            self.deselectMeshObject(self.selMeshObject, self.selMeshTag)
        if self.selVertex:
            self.deselectVertex(self.selVertex, self.selVertexTag)


    def selectZoom(self, tag):
        """
        Zoom in or out of the current mesh view
        """

        from locale import atof
        if type(tag) is str:
            fraction = atof(tag)
        else:
            fraction = tag
        self.SCALE *= fraction
        self.scrolledcanvas.scale(ALL, 0, 0, fraction, fraction)

        # Redraw all of the vertices, holes and regions,
        #so the squares representing vertices
        # don't get bigger
        vertices = self.mesh.getUserVertices()
        holes = self.mesh.getHoles()
        regions = self.mesh.getRegions()
        MeshObjects  = vertices + holes + regions

        # make a list of tags to delete
        guiIDs = [getattr(MeshObjects[i],'guiID') for i in range(len(MeshObjects))]
        self.canvas.delete(*guiIDs)
        for obj in MeshObjects:
            if self.selVertex == obj:
                obj.draw(self.canvas,obj.guiID,  scale =self.SCALE ,colour= VERT_SELECT_ADDING_SEG_COLOR)
            elif self.selMeshObject == obj:
                obj.draw(self.canvas,obj.guiID,  scale =self.SCALE ,colour= SELECT_COLOR)
            else:
                obj.draw(self.canvas,obj.guiID,  scale =self.SCALE )
        top, bottom = self.scrolledcanvas.xview()
        xcenter  = (top + bottom)/2
        xdiff =  xcenter - top
        xcnew = xcenter - xdiff/fraction

        top, bottom = self.scrolledcanvas.yview()
        ycenter = (top + bottom)/2
        ydiff = ycenter - top
        ycnew = ycenter - ydiff/fraction

        self.scrolledcanvas.resizescrollregion()
        # update so the moveto calls will work...
        self.scrolledcanvas.update()
        # but calling update now does make things jerky
        self.canvas.xview_moveto(xcnew)
        self.canvas.yview_moveto(ycnew)


    def windowAddVertex (self, parent):
        """
        add a vertex using a window and entering x y values.

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """

        dialog = AddVertexDialog(self.canvas)
        if dialog.xyValuesOk:
            log.critical(str(dialog.x))
            log.critical(str(dialog.y))
            #self.drawVertex(dialog.x*self.SCALE,dialog.y*self.SCALE,None)
            self.drawVertex(dialog.x, dialog.y, None)
            #Since the new vertex may be off screen
            #self.ResizeToFit()
        else:
            log.critical("bad values")

    def windowDelete(self, parent):
        self.DeleteSelectedMeshObject(None)

    def windowDefault (self, parent):
        """

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        # self.UserMesh is a vMesh instance
        self.UserMesh.defaultWindow(self.canvas, self.curModeClass)

    def windowEdit (self, parent):
        """

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        if self.selMeshObject:
            self.UserMeshChanged = self.UserMesh.editWindow(self.canvas,
                                     self.selMeshObject,
                                     self.UserMeshChanged)

    def auto_segmentButton (self, parent):
        self.auto_segment()


    def auto_segmentGiveAlphaButton (self, parent):
        dialog = auto_segmentDialog(self.canvas, self.meshLastAlpha)
        if dialog.use_optimum.get() == SET_ALPHA:
            if dialog.alphaValueOk:
                self.auto_segment(alpha = dialog.alpha,
                                 raw_boundary=dialog.raw_boundary.get(),
                                 remove_holes=dialog.remove_holes.get(),
                                 smooth_indents=dialog.smooth_indents.get(),
                                 expand_pinch=dialog.expand_pinch.get())
            else:
                 showerror('pMesh',
                      'Bad alpha value.')
        else:
            self.auto_segment(raw_boundary=dialog.raw_boundary.get(),
                             remove_holes=dialog.remove_holes.get(),
                             smooth_indents=dialog.smooth_indents.get(),
                             expand_pinch=dialog.expand_pinch.get())


    def auto_segment (self, alpha = None,
                                 raw_boundary=True,
                                 remove_holes=False,
                                 smooth_indents=False,
                                 expand_pinch=False ):
        """
        add Segments to bound all vertices

        """
        if len(self.mesh.getUserVertices()) >= 3:
            try:
                newsegs, ObjectsToVisuallyDelete, self.meshLastAlpha = \
                     self.mesh.auto_segment(alpha=alpha,
                                           remove_holes=remove_holes,
                                           smooth_indents=smooth_indents,
                                           expand_pinch=expand_pinch)
            except AlphaError:
                showerror('pMesh',
                          'Unable to auto_segment.')
            else:

                for drawOb in ObjectsToVisuallyDelete:
                    self.UserMesh.unvisualise(drawOb, self.canvas)

                for segment in newsegs:
                    self.serial +=1
                    self.uniqueID = 'M*%d' % self.serial
                    self.Segments.visualise(segment,
                                            self.uniqueID,
                                            self.canvas,
                                            self.SCALE)

        else:
            showerror('pMesh',
                      'Three or more vetices are needed to be able to auto_segment.')


    def auto_segmentFilter (self):
        dialog = auto_segmentFilterDialog(self.canvas)
        newsegs, ObjectsToVisuallyDelete, self.meshLastAlpha = \
                 self.mesh.auto_segmentFilter(raw_boundary=dialog.raw_boundary.get(),
                             remove_holes=dialog.remove_holes.get(),
                             smooth_indents=dialog.smooth_indents.get(),
                             expand_pinch=dialog.expand_pinch.get())

        for drawOb in ObjectsToVisuallyDelete:
            self.UserMesh.unvisualise(drawOb, self.canvas)

        for segment in newsegs:
            self.serial +=1
            self.uniqueID = 'M*%d' % self.serial
            self.Segments.visualise(segment,
                                    self.uniqueID,
                                    self.canvas,
                                    self.SCALE)

    def joinVerticesButton (self, parent):
        ans = askyesno("", "This cannot be undone. Are you sure?")
        if ans:
            self.joinVertices()

    def joinVertices (self):
        """
        add Segments to connect all vertices

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        if len(self.mesh.getUserVertices()) >= 3:
            newsegs = self.mesh.joinVertices()
            for segment in newsegs:
                self.serial +=1
                self.uniqueID = 'M*%d' % self.serial
                self.Segments.visualise(segment,
                                        self.uniqueID,
                                        self.canvas,
                                        self.SCALE)
        else:
            showerror('pMesh',
                      'Three or more vetices are needed to be able to join vertices.')

    def windowMeshGen (self, parent):
        """
        The parent attribute isn't used by this function.
        need to understand toolbarbutton.py to know how to
        get rid of it.
        """
        # Put exceptions round things.
        #catch failure in self.mesh.generateMesh
        dialog = MeshGenDialog(self.canvas,
                               self.MeshMinAngle,
                               self.MeshMaxArea,
                               self.MeshnumTriangles,
                               self.MeshMaxAreaLast)

        if dialog.ValuesOk:
            log.critical(str(dialog.minAngle))
            log.critical(str(dialog.maxArea))

            self.clearSelections()
            self.canvas.delete(ALL)
            if dialog.goodMaxArea == True:
                self.MeshMinAngle = dialog.minAngle
                self.MeshMaxArea = dialog.maxArea
                self.MeshMaxAreaLast = True

                self.mesh = self.MeshGenAreaAngle (dialog.minAngle,
                                          dialog.maxArea,
                                          self.mesh)
            elif dialog.goodNumTriangles == True:
                self.MeshMinAngle = dialog.minAngle
                self.MeshnumTriangles  = dialog.numTriangles
                self.MeshMaxAreaLast = False

                self.mesh = self.MeshGenAreaNumTriangles (dialog.minAngle,
                                          dialog.numTriangles,
                                          self.mesh)
            else:
                pass
            log.critical("userMeshChanged = False")
            self.UserMeshChanged = False
            self.visualiseMesh(self.mesh)
            log.critical("Mesh Generation finished")

    def MeshGenAreaAngle (self, minAngle, maxArea, mesh):
        """
        Generate a mesh, given a minAngle and max area
        """
        tempMesh = mesh
        try:
            tempMesh.generateMesh(mode = "pzq"+str(minAngle)
                                  +"a"+str(maxArea)
                                  +"a") #So areas for regions will be used
        except AttributeError : # can't catch PyEval_RestoreThread
            # This doesn't catch tempMesh.generateMesh failing
            tempMesh = mesh
        return tempMesh


    def MeshGenAreaNumTriangles (self, minAngle, numTriangles, mesh):
        """
        Generate a mesh, given a minAngle and rough # of triangles
        """
        #get big triangles
        #calc area
        #calc max triangle area
        #
        tempMesh = mesh
        try:
            tempMesh.generateMesh("pzq1")
        except AttributeError : # can't catch PyEval_RestoreThread
            # This doesn't catch tempMesh.generateMesh failing
            pass
        meshArea = 0
        meshArea = tempMesh.tri_mesh.calc_mesh_area()
        maxArea = meshArea/numTriangles


        return self.MeshGenAreaAngle (minAngle,
                                      maxArea,
                                      self.mesh)

    def mouseDown(self, event):
        """
        On a mouse down event, depending on the current state,
        either add a vertex or a seg etc
        """
        self.curObject = None
        self.lastx = self.startx = self.canvas.canvasx(event.x)
        #The screen canvas has y 'flipped'.  -1* unflips it
        self.lasty = self.starty = -1*self.canvas.canvasy(event.y)
        log.critical("----------------------")
        log.critical(f"mouseDown x: {self.lastx}, y: {self.lasty}")
        self.mouseDownCurFunc( self.lastx,
                               self.lasty,event) #!!! remove the event?
                                                 # do last

    def displayCoords(self, event):
        messageBar = self.messageBar()
        disp_x = str(int(self.canvas.canvasx(event.x)))
        disp_y = str(int(-1*self.canvas.canvasy(event.y)))
        messageBar.helpmessage(disp_x + "," + disp_y)


    def rightMouseUp(self, event):
        """
        On a right mouse button up event select the nearest object.
        """
        found=False
        if event.widget.find_withtag(CURRENT): # if there's a widget with a tag
            tag = self.canvas.gettags(CURRENT)[0] # get a list of them
            log.critical("tag %s" % str(tag))  #tags ('M*1008', 'current')
            if tag[:2] == 'M*':   #!!! this can be removed when there are
                #    only mesh objects on screen
                #key, value = string.split(tag, '*')
                objectID = tag
                log.critical("Found!! objectID: %s" % str(objectID))

                meshObjects = self.getAllUserMeshObjects()
                # It may be a triangle, which is ignored
                if meshObjects.hasKey(objectID):
                    selMeshObject = meshObjects.getMeshObject(objectID)
                    found = True
                    log.critical("Found! selMeshObject: %s"
                                 % str(selMeshObject))
                    #Only select one object at a time
                    if self.selMeshObject:
                        self.deselectMeshObject(self.selMeshObject, self.selMeshTag)
                    self.selectMeshObject(selMeshObject,objectID)

    def getAllUserMeshObjects(self):
        return self.UserMesh

    def DeleteSelectedMeshObject(self, event):
        """
        if an object is selected, delete it.
        """
        if self.selMeshObject:
            #an object is selected
            #first deselect the vertex, for selecting a segment
            if self.selVertex:
                self.deselectVertex(self.selVertex, self.selVertexTag)
            ObjectsToVisuallyDelete = self.mesh.deleteMeshObject (self.selMeshObject)
            for drawOb in ObjectsToVisuallyDelete:
                self.UserMesh.unvisualise(drawOb, self.canvas)

            self.selMeshObject = None
            self.selMeshTag = None

    def selectMeshObject(self, meshObject, objectID):
        """
        selected a mesh object.
        """
        self.canvas.delete(objectID)
        self.selMeshObject = meshObject
        self.selMeshTag = objectID
        meshObject.draw(self.canvas,objectID, scale =self.SCALE ,colour = SELECT_COLOR)

    def deselectMeshObject(self, meshObject, objectID):
        """
        deselected a mesh object.
        """
        self.canvas.delete(objectID)
        self.selMeshObject = None
        self.selMeshTag = None
        if isinstance(meshObject, mesh.Segment):
            meshObject.draw(self.canvas,objectID,
                        scale =self.SCALE ,colour = SEG_COLOUR)
        else:
            meshObject.draw(self.canvas,objectID,
                        scale =self.SCALE )

    def drag(self,x,y,event):
        """
        Hack function.  called when in select and left mouse goes down
        """
        pass


    def drawEastingNorthingVertex(self,x,y,event):
        """
        draw a vertex object, plus add it to the mesh data structure

        event isn't used
        """
        self.serial +=1
        self.uniqueID = 'M*%d' % self.serial
        #x_scaled =  self.SCALE*x
        #y_scaled = -1*self.SCALE*y
        vert = self.Vertices.draw(x-self.mesh.geo_reference.get_xllcorner,
                                  y-self.mesh.geo_reference.get_yllcorner,
                                  self.mesh,
                                  self.uniqueID,
                                  self.SCALE,
                                  self.canvas,
                                  event) #FIXME why is event passed on.
        self.UserMeshChanged = True
        return vert

    def drawVertex(self,x,y,event):
        """
        draw a vertex object, plus add it to the mesh data structure

        event isn't used
        """
        self.serial +=1
        self.uniqueID = 'M*%d' % self.serial
        
        #x_scaled =  self.SCALE*x
        #y_scaled = -1*self.SCALE*y
        vert = self.Vertices.draw(x,y,self.mesh,self.uniqueID,self.SCALE,self.canvas,event)
        self.UserMeshChanged = True
        return vert

    def drawHole(self,x,y,event):
        """
        draw a hole object, plus add it to the mesh data structure

        event isn't used
        """
        self.serial +=1
        self.uniqueID = 'M*%d' % self.serial
        self.userMeshChanged = True
        hole = self.Holes.draw(x,y,self.mesh,self.uniqueID,self.SCALE,self.canvas,event)
        return hole

    def drawRegion(self,x,y,event):
        """
        draw a region object, plus add it to the mesh data structure

        event isn't used
        """
        self.serial +=1
        self.uniqueID = 'M*%d' % self.serial
        region = self.Regions.draw(x,y,self.mesh,self.uniqueID,self.SCALE,self.canvas,event)
        return region

    def selectSegmentPoint(self,x,y, event):
        """
        logic when selecting a vertex object to add a segment
        """
        found=False
        if event.widget.find_withtag(CURRENT): # if there's a widget with a tag
            [tag,string] = self.canvas.gettags(CURRENT) # get a list of them
            log.critical("tag %s" % str(tag))  #tags ('M*1008', 'current')
            objectID = tag
            if self.Vertices.hasKey(objectID): #isinstance(self.meshObjects[objectID],mesh.Vertex):
                vertex = self.Vertices.getMeshObject(objectID)
                found = True
                log.critical("Found! vertex: %s" % str(vertex))

        if found and self.selVertex == vertex:
            log.critical("The selected vertex has already been selected")
            #The selected vertex has already been selected
            # therefore deselect it
            self.deselectVertex(self.selVertex, self.selVertexTag)
            found = False

        if found:
            #A vertex has been selected!
            if self.selVertex:
                if self.mesh.isUserSegmentNew(self.selVertex,vertex):
                    #vertex is the 2nd vertex
                    self.drawSegment(vertex,self.selVertex)
                    self.deselectVertex(self.selVertex, self.selVertexTag)
                    self.selectVertex(vertex,objectID)
            else:
                log.critical("vertex is the 1st vertex")
                #vertex is the 1st vertex
                self.selectVertex(vertex,objectID)
        else:
            log.critical(" There are no widgets.  This happen's too much")


    def selectVertex(self, vertex,objectID):
        """
        select a vertex object when adding a segment
        """
        self.canvas.delete(objectID)
        self.selVertex = vertex
        self.selVertexTag = objectID
        vertex.draw(self.canvas,objectID, scale =self.SCALE ,colour = VERT_SELECT_ADDING_SEG_COLOR)

    def deselectVertex(self, vertex,objectID):
        """
        deselect a vertex object when adding a segment
        """
        self.canvas.delete(objectID)
        self.selVertex = None
        self.selVertexTag = None
        vertex.draw(self.canvas,objectID,  scale =self.SCALE )

    def drawSegment(self,v1,v2):
        """
        Create a seg object, draw it and add it to the mesh data structure
        """
        self.serial +=1
        self.uniqueID = 'M*%d' % self.serial
        self.userMeshChanged = True
        seg = self.Segments.draw(v1,v2,self.mesh,self.uniqueID,self.SCALE,self.canvas,None)
        return seg
    def printGeoReference(self):
        try:
            log.critical("geo reference %s" % str(self.mesh.geo_reference))
        except:
            log.critical("no geo reference")

    def visualiseMesh(self,mesh):
        """
        visualise vertices, segments, triangulation, holes
        """
        if self.Visualise:
            try:
                mesh.tri_mesh.draw_triangulation(self.canvas,
                                             scale = self.SCALE)
            except:
                pass

        for segment in mesh.getUserSegments():
            self.serial +=1
            self.uniqueID = 'M*%d' % self.serial
            self.Segments.visualise(segment,
                                    self.uniqueID,
                                    self.canvas,
                                    self.SCALE)
        for vertex in mesh.getUserVertices():
            self.serial +=1
            self.uniqueID = 'M*%d' % self.serial
            self.Vertices.visualise(vertex,
                                    self.uniqueID,
                                    self.canvas,
                                    self.SCALE)

        for hole in mesh.getHoles():
            self.serial +=1
            self.uniqueID = 'M*%d' % self.serial
            self.Holes.visualise(hole,
                                    self.uniqueID,
                                    self.canvas,
                                    self.SCALE)
        for region in mesh.getRegions():
            self.serial +=1
            self.uniqueID = 'M*%d' % self.serial
            self.Regions.visualise(region,
                                    self.uniqueID,
                                    self.canvas,
                                    self.SCALE)

    def obsolete_normalise4ObjMesh(self):
        if self.mesh:
            self.clearSelections()
            self.canvas.delete(ALL)
            self.mesh.normaliseMesh(400,-200,20)
            self.visualiseMesh(self.mesh)
            self.ResizeToFit()
            self.ResizeToFit()

    def obsolete_nnormaliseMesh(self):
        if self.mesh:
            self.clearSelections()
            self.canvas.delete(ALL)
            self.mesh.normaliseMesh(1,0,1)
            self.visualiseMesh(self.mesh)
            self.ResizeToFit()
            self.ResizeToFit()


    def clearMesh(self):
        """Clear the current mesh object, and the canvas """

        self.clearSelections()
        self.canvas.delete(ALL)
        self.deleteMesh()
        self.initData()
        self.createMesh()

    def exportObj(self):
        fileType = "obj"
        fileTypeDesc = "obj mesh"

        ofile = tkinter.filedialog.asksaveasfilename(filetypes=[(fileTypeDesc,
                                                           fileType),
                                                          ("All Files", "*")])
        if ofile:
            addOn = "." + fileType
            jumpback = - len(addOn)
            if ofile[jumpback:] != addOn:
                ofile = ofile + addOn
            try:
                self.mesh.exportASCIIobj(ofile)
            except IOError:
                showerror('Export ASCII file',
                                   'Can not write to file.')
            except RuntimeError:
                showerror('Export ASCII file',
                                   'No triangulation to export.')



    def ImportUngenerate(self):
        ofile = tkinter.filedialog.askopenfilename(initialdir=self.currentPath,
                   filetypes=[ ("ungenerated polygon information", "txt"),
                                           ("All Files", "*")])
        if ofile == "":
            # The user cancelled the loading action
            return

        try:
            self.clearSelections()
            self.canvas.delete(ALL)
            dict = mesh.importUngenerateFile(ofile)
            self.mesh.addVertsSegs(dict)

        except SyntaxError:
            # This is assuming that the SyntaxError is thrown in
            # importUngenerateFile
            showerror('File error',
                      ofile + ' is not in the correct format.')
        except IOError:
            #!!! this error type can not be thrown?
            showerror('File error',
                      'file ' + ofile + ' could not be found.')
        except RuntimeError:
            showerror('File error',
                  'file ' + ofile + ' has an unknown file type.')

        self.visualiseMesh(self.mesh)
        self.ResizeToFit()

    def exportASCIIsegmentoutlinefile(self):

        ofile = tkinter.filedialog.asksaveasfilename(initialdir=self.currentPath,
                                               filetypes=[("mesh", "*.tsh *.msh"),
                                             ("All Files", "*")])

        if ofile:
            # .tsh is the default file format
            if (ofile[-4:] == ".tsh" or ofile[-4:] == ".msh"):
                self.currentFilePathName = ofile
            else:
                self.currentFilePathName = ofile + ".tsh"

            try:
                self.mesh.exportASCIIsegmentoutlinefile(ofile)
            except IOError: #FIXME should this function be throwing any errors?
                showerror('Export ASCII file',
                                   'Can not write to file.')

    def exportPointsFile(self):
        ofile = tkinter.filedialog.asksaveasfilename(initialdir=self.currentPath,
                                         filetypes=[("point files",
                                                     "*.csv *.txt *.pts"),
                                                ("All Files", "*")])
        if ofile:
            # .csv is the default file format
            if (ofile[-4:] == ".csv" or ofile[-4:] == ".pts"):
                self.currentFilePathName = ofile
            else:
                self.currentFilePathName = ofile + ".csv"

            try:
                self.mesh.exportPointsFile(ofile)
            except IOError:
                showerror('Export ASCII file',
                                   'Can not write to file.')

    def importFile(self):
        """
        import mesh data from a variety of formats (currently 2!)
        """
        log.critical("self.currentPath %s" % str(self.currentPath))
        ofile = tkinter.filedialog.askopenfilename(initialdir=self.currentPath,
                                             filetypes=[ ("Mesh",
                                                          "*.tsh *.msh"),
                                                         ("Points",
                                                          "*.csv *.txt *.pts"),
                                                         ("All Files", "*")])
        if ofile == "":
            # The user cancelled the loading action
            return

        try:
            newmesh = mesh.importMeshFromFile(ofile)
            self.currentPath, dummy = os.path.split(ofile)
            self.currentFilePathName = ofile
            self.clearMesh()
            self.mesh = newmesh

            #FIXME - to speed things up, don't visualise the mesh
            # use ResizeToFitWrapper
            self.visualiseMesh(self.mesh)
            self.ResizeToFit()

        except IOError:
            #!!! this error type can not be thrown?
            showerror('File error',
                      'file ' + ofile + ' could not be loaded.')

        except RuntimeError:
            showerror('File error',
                  'file ' + ofile + ' has an unknown file type.')
        # Could not get the file name to showup in the title
        #appname =  ofile + " - " + APPLICATION_NAME

        except anuga.load_mesh.loadASCII.TitleAmountError:
            showerror('File error',
                  'file ' + ofile + ' has a bad title line (first line).')


    def ResizeToFitWrapper(self, Parent):
        """
        The parent attribute isn't used by this function.
        need to understand toolbarbutton.py to know how to
        get rid of it.
        """
        self.ResizeToFit()

    def ResizeToFit(self):
        """Visualise the mesh so it fits in the window"""
        if self.mesh.getUserVertices() == []:
            return #There are no vertices!
        # Resize the window
        self.scrolledcanvas.resizescrollregion()
        # I need this so the xview values are correct
        self.scrolledcanvas.update()

        xtop, xbottom = self.scrolledcanvas.xview()
        ytop, ybottom = self.scrolledcanvas.yview()
        xdiff = xbottom-xtop
        ydiff = ybottom-ytop
        if xdiff == 1 and xdiff == 1:
            #The mesh might be too small.
            #Make it too large, then resize
            #!!! Recursive hack.  Should be a better way
            fraction = 50
            self.SCALE *= fraction
            self.scrolledcanvas.scale(ALL, 0, 0, fraction, fraction)
            self.ResizeToFit()
        else:
            # without 0.99 some of the mesh may be off screen
            fraction = 0.99*min(xdiff,ydiff)
            self.selectZoom(fraction)

    def saveDrawing(self):
        """
        Save the current drawing
        """
        if (self.currentFilePathName[-4:] != ".tsh" or
            self.currentFilePathName[-4:] != ".msh"):
            # force user to choose a name
            self.saveAsDrawing()
        else:
            self.exportASCIItriangulationfile(self.currentFilePathName)

    def saveAsDrawing(self):
        """
        Save the current drawing, prompting for a file name
        """
        ofile = tkinter.filedialog.asksaveasfilename(initialdir=self.currentPath,
                                               filetypes=[("mesh", "*.tsh *.msh"),
                                             ("All Files", "*")])

        if ofile:
            # .tsh is the default file format
            if (ofile[-4:] == ".tsh" or ofile[-4:] == ".msh"):
                self.currentFilePathName = ofile
            else:
                self.currentFilePathName = ofile + ".tsh"
            self.exportASCIItriangulationfile(self.currentFilePathName)

    def exportASCIItriangulationfile(self,currentFilePathName):
        """
        Have a warning prompt when saving a mesh where the generated mesh is
        different from the user mesh - eg boundary tags that aren't carried
        thru. Warning ~"Mesh not generated after changes.  Generate mesh?  "
        - cancel, don't gen, don't save.  Yes - generate mesh, go to save
        screen.  No - goto save screen.  To implement this need to know when
        the user has done a change, and the mesh hasn't been generated.  If
        there is no generated mesh do not prompt.
        """
        if (self.UserMeshChanged) and self.mesh.isTriangulation():

            m = _show("Warning",
                                   "A triangulation has not been generated, after mesh changes.  Generate triangulation before saving?",
                                   icon=QUESTION,
                                   type=YESNOCANCEL)
            if m == "no":
                self.mesh.export_mesh_file(currentFilePathName)
                self.UserMeshChanged = False
            elif m == "cancel":
                pass
            elif m == "yes":
                self.windowMeshGen(None)
                self.mesh.export_mesh_file(currentFilePathName)
        else:
            self.mesh.export_mesh_file(currentFilePathName)
            self.UserMeshChanged = False

    def initData(self):
        """
        Initialise various lists and flags
        """
        self.serial      = 1000
        self.currentFilePathName = 'untitled'

        # these are attributes I've added
        self.SCALE       = 1
        self.selVertex   = None     #The last vertex selected, in draw seg mode
        self.selVertexTag   = None     # The selected vertex drawn object tag
        self.mesh        = None
        self.MeshMinAngle = 20
        self.MeshMaxArea = 200
        self.MeshnumTriangles = 20
        self.MeshMaxAreaLast = False
        self.selMeshObject = None  # The mesh object selected in the current mode
        self.selMeshTag = None
        mesh.Segment.set_default_tag("")
        self.UserMeshChanged = False
        self.meshLastAlpha = None

        self.Visualise = True #Is the mesh shown or not?

    def ipostscript(self):
        """
        Print the canvas as a postscript file
        """
        ofile = tkinter.filedialog.asksaveasfilename(initialdir=self.currentPath,
                                               filetypes=[("postscript", "ps"),
                                                          ("All Files", "*")])
        if ofile:
            if ofile[-3:] != ".ps":
                ofile = ofile + ".ps"
            postscript = self.canvas.postscript()
            fd = open(ofile, 'w')
            fd.write(postscript)
            fd.close()

    def close(self):
        self.quit()

    def createInterface(self):
        """
        Call all functions that create the GUI interface
        """
        self.initData()
        self.createMesh()
        AppShell.AppShell.createInterface(self)
        self.createButtons()
        self.createMenus()
        self.createBase()
        self.createTools()
        self.createZooms()
        self.createEdits()
        self.createVisualiseIcons()
        #print "FIX THIS BEFORE "
        #self.addCylinders() # !!!DSG start pmesh with a triangle
        #self.selectFunc('Pointer')
        self.selectFunc('Add-Vertex')
        self.currentPath = os.getcwd()

    def loadtestmesh(self,ofile):
        """
        debugging script to test loading a file
        """
        fd = open(ofile)
        a = mesh.Vertex (-10.0, 0.0)
        d = mesh.Vertex (0.0, 4.0)
        f = mesh.Vertex (4.0,0.0)
        g = mesh.Vertex (-5.0,5.0)

        s1 = mesh.Segment(a,d)
        s2 = mesh.Segment(d,f)
        s3 = mesh.Segment(a,f)

        r1 = mesh.Region(0.3, 0.3)

        m = mesh.Mesh(userVertices=[a,d,f,g], userSegments=[s1,s2,s3], regions=[r1] )

        fd.close()
        log.critical('returning m')
        return oadtestmesh(ofile)

class  AddVertexDialog(Dialog):
    """
    Dialog box for adding a vertex by entering co-ordindates
    """
    def body(self, master):
        """
        GUI description
        """
        self.title("Add New Vertex")

        Label(master, text='X position:').grid(row=0, sticky=W)
        Label(master, text='Y position:').grid(row=1, sticky=W)

        self.xstr   = Entry(master, width = 16, name ="entry")
        self.ystr  = Entry(master, width = 16)

        self.xstr.grid(row=0, column=1, sticky=W)
        self.ystr.grid(row=1, column=1, sticky=W)
        self.xstr.focus_force()
        self.x  = 0
        self.y  = 0
        self.xyValuesOk = False


    def apply(self):
        """
        check entered values
        """
        try:
            self.x = float(self.xstr.get())
            self.y = float(self.ystr.get())
            self.xyValuesOk = True

        except ValueError:
            showerror('Bad Vertex values',
                                   'X Y values are not numbers.')


class  auto_segmentDialog(Dialog):
    """
    Dialog box for adding segments
    """
    def __init__(self, parent, alpha):
        self.alpha = alpha
        Dialog.__init__(self, parent)

    def body(self, master):
        """
        GUI description
        """
        self.title("Automatically Add Segments")

        self.use_optimum = IntVar()
        self.use_optimum.set(AUTO) # should initialise the radio buttons.
                                   #  It doesn't

        #self.use_optimum.set(NO_SELECTION)
        self.ck = Radiobutton(master, value = AUTO, variable=self.use_optimum)
        self.ck.grid(row=1, column=0)
        Label(master, text='Use optimum alpha').grid(row=1, column=1, sticky=W)

        self.ck2 = Radiobutton(master, value = SET_ALPHA,
                               variable=self.use_optimum)
        self.ck2.grid(row=2, column=0)

        Label(master, text='alpha:').grid(row=2, column=1, sticky=W)
        if (self.alpha):
            alphaVar = StringVar()
            alphaVar.set(self.alpha)
            self.alpha_str  = Entry(master,
                                     textvariable = alphaVar,
                                     width = 16, name ="entry")
        else:
            self.alpha_str = Entry(master, width = 16, name ="entry")

        self.alpha_str.grid(row=2, column=3, sticky=W)

        #boundary type buttons
        self.raw_boundary = IntVar()
        self.remove_holes = IntVar()
        self.smooth_indents = IntVar()
        self.expand_pinch = IntVar()
        self.ck3 = Checkbutton(master, state=NORMAL,
                               variable=self.raw_boundary)
        self.ck3.grid(row=3, column=0)
        Label(master, text='Raw boundary').grid(row=3, column=1, sticky=W)
        #
        self.ck4 = Checkbutton(master, state=NORMAL,
                               variable=self.remove_holes)
        self.ck4.grid(row=4, column=0)
        Label(master, text='Remove small holes').grid(row=4,column=1, sticky=W)
        #
        self.ck5 = Checkbutton(master,state=NORMAL,
                               variable=self.smooth_indents)
        self.ck5.grid(row=5, column=0)
        Label(master,
              text='Remove sharp indents').grid(row=5, column=1, sticky=W)
        #
        self.ck6 = Checkbutton(master,state=NORMAL,
                               variable=self.expand_pinch)
        self.ck6.grid(row=6, column=0)
        Label(master,
              text='Remove pinch off').grid(row=6, column=1,  sticky=W)


        self.alpha  = 0
        self.alphaValueOk = False


    def apply(self):
        """
        check entered values
        """
        try:
            self.alpha = float(self.alpha_str.get())
            self.alphaValueOk = True

        except ValueError:
            pass
            #showerror('Bad Alpha value',
            #                       'Alpha is negative.')


class  auto_segmentFilterDialog(Dialog):
    """
    Dialog box for adding segments
    """
    def __init__(self, parent):
        Dialog.__init__(self, parent)

    def body(self, master):
        """
        GUI description
        """
        self.title("Automatically Add Segments")

        self.use_optimum = IntVar()
        self.use_optimum.set(AUTO) # should initialise the radio buttons.
                                   #  It doesn't
        self.boundary_type = IntVar()

        #boundary type buttons
        self.raw_boundary = IntVar()
        self.remove_holes = IntVar()
        self.smooth_indents = IntVar()
        self.expand_pinch = IntVar()

        self.ck3 = Checkbutton(master, state=NORMAL,
                               variable=self.raw_boundary)
        self.ck3.grid(row=3, column=0)
        Label(master, text='Raw boundary').grid(row=3, column=1, sticky=W)
        #
        self.ck4 = Checkbutton(master, state=NORMAL,
                               variable=self.remove_holes)
        self.ck4.grid(row=4, column=0)
        Label(master, text='Remove small holes').grid(row=4,column=1, sticky=W)
        #
        self.ck5 = Checkbutton(master,state=NORMAL,
                               variable=self.smooth_indents)
        self.ck5.grid(row=5, column=0)
        Label(master,
              text='Remove sharp indents').grid(row=5, column=1, sticky=W)
        #
        self.ck6 = Checkbutton(master,state=NORMAL,
                               variable=self.expand_pinch)
        self.ck6.grid(row=6, column=0)
        Label(master,
              text='Remove pinch off').grid(row=6, column=1,  sticky=W)




class  MeshGenDialog(Dialog):
    """
    Dialog box for generating a mesh
    """
    # initial values, hard coded.
    # should be values associated with the current mesh
    lastMinAngle = 20
    lastMaxArea  = 100


    def __init__(self,
                 parent,
                 minAngle,
                 maxArea,
                 numTriangles,
                 MeshMaxAreaLast):
        self.minAngle = minAngle
        self.maxArea = maxArea
        self.numTriangles = numTriangles
        self.MeshMaxAreaLast = MeshMaxAreaLast

        Dialog.__init__(self, parent)


    def body(self, master):
        """
        GUI description
        """
        self.title("Generate Mesh")

        Label(master,
              text='Minimum Angle(0 - 40):').grid(row=0, sticky=W)
        Label(master,
              text='Angles>33 may not converge').grid(row=1, sticky=W)
        Label(master, text='Maximum Area:').grid(row=2, sticky=W)
        Label(master, text='OR # of triangles:').grid(row=3, sticky=W)


        minAngleVar = StringVar()
        minAngleVar.set(self.minAngle)
        self.minAnglestr   = Entry(master,
                                   width = 16,
                                   textvariable = minAngleVar,
                                   takefocus=1)
        if (self.MeshMaxAreaLast):
            maxAreaVar = StringVar()
            maxAreaVar.set(self.maxArea)
            self.maxAreastr  = Entry(master,
                                     textvariable = maxAreaVar,
                                     width = 16)
            self.numTrianglesstr  = Entry(master,
                                          width = 16)
        else:
            self.maxAreastr  = Entry(master,
                                     width = 16)
            self.maxAreastr.focus_force()
            numTrianglesVar = StringVar()
            numTrianglesVar.set(self.numTriangles)
            self.numTrianglesstr  = Entry(master,
                                          textvariable = numTrianglesVar,
                                          width = 16)


        self.minAnglestr.grid(row=0, column=1, sticky=W)
        self.maxAreastr.grid(row=2, column=1, sticky=W)
        self.numTrianglesstr.grid(row=3, column=1, sticky=W)

        self.numTriangles = 0
        self.minAngle  = 0
        self.maxArea  = 0
        self.ValuesOk = False


    def apply(self):
        """
        check entered values
        """
        self.goodMaxArea = self.goodNumTriangles = True
        self.ValuesOk = True
        try:
            self.minAngle = float(self.minAnglestr.get())
            MeshGenDialog.lastMinAngle =self.minAngle
        except ValueError:
            self.ValuesOk = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')

        try:
            self.maxArea = float(self.maxAreastr.get())
            MeshGenDialog.lastMaxArea =self.maxArea
        except ValueError:
            self.goodMaxArea = False

        try:
            self.numTriangles = int(self.numTrianglesstr.get())
            MeshGenDialog.lastNumTriangles =self.numTriangles
        except ValueError:
            self.goodNumTriangles= False

        if self.goodMaxArea == False and self.goodNumTriangles == False:
            self.ValuesOk = False
            showerror('Bad mesh generation values',
                      'Values are not numbers.')

        if self.goodMaxArea == True and self.goodNumTriangles == True:
            self.ValuesOk = False
            showerror('Bad mesh generation values',
                      'Give a maximum area OR number of triangles, not both.')

        try:
            # value checking
            if self.minAngle <0.0 or self.minAngle >40.0:
                raise IOError
            if self.goodMaxArea == True and self.maxArea <0.0:
                raise IOError
            if self.goodNumTriangles == True and self.numTriangles <=0:
                raise IOError

        except IOError:
            self.ValuesOk = False
            showerror('Bad mesh generation values',
                                   'Values are out of range.')

if __name__ == '__main__':
    draw = Draw()
    draw.run()
    #profile.run('draw.run()', 'pmeshprof')
    
