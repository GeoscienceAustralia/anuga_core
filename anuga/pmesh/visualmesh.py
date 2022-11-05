

try:
    from . import mesh
except:
    import mesh

from tkinter.simpledialog import Dialog,askfloat, askinteger, askstring
from tkinter import  FALSE,TRUE, Frame,X, LEFT,YES,BOTH,ALL,Widget,CURRENT, Label,W, Entry, E, ACTIVE, NORMAL, StringVar
from tkinter.messagebox import showerror

class vAbstract(object):
    """
    Define the visualisation of a mesh object
    """
    def __init__(self, association):
        # The association is the class oncanvas will hold ie mesh.Vertex
        self.association = association 
        self.oncanvas = {} 

    def visualise(self, meshobject, uniqueID, canvas, scale):
        """
        draws a mesh object and add it the the canvas data structure.
        Note: a guiID attribute is added to each point object thats visualised.
        This is used in zooming.

        meshobject is from mesh.py, eg vertex, hole
        """
        meshobject.guiID = uniqueID 
        meshobject.draw(canvas, uniqueID, scale =scale)
        self.oncanvas[uniqueID] = meshobject
        return meshobject

    def unvisualise(self, meshobject, canvas):
        """
        delete a visual mesh object and
        delete it from the canvas mesh data structure
        """
        canvas.delete(meshobject.guiID)
        del self.oncanvas[meshobject.guiID]

    def hasKey(self, key):
        """
        returns true if the key is in the canvas list
        """
        if key in self.oncanvas:
            return True
        else:
            return False
 
    def getMeshObject(self, key):
        """
        returns the mesh object,
        given it's key.

        precon:  The key is valid
        """
        return self.oncanvas[key]

    def editWindow(self, canvas, selMeshObject,userMeshChanged):
        print("Not implimented")
        return userMeshChanged

    def defaultWindow(self, canvas):
        print("Not implimented")

    def draw (self):
        """
        This method must be overriden
        """
        pass

    def __add__(self, other):
        return VMesh([self,other])

    def __repr__(self):
        return str(self.__class__)
    
class vPoints(vAbstract):
    """
    Define the visualisation of a point
    """

    def draw(self,x,y,mesh,uniqueID,scale,canvas,event):
        """
        draw a point object, plus add it to the mesh data structure

        event isn't used
        """
        point = mesh.addUserPoint(self.association,x/scale,y/scale)
        self.visualise(point, uniqueID, canvas, scale)
        return point
 
class vRegions(vPoints):
    """
    Define the visualisation of a region
    """ 
    def editWindow(self, canvas, selMeshObject,userMeshChanged):
        dialog = EditRegionDialog (canvas,
                                   selMeshObject.getTag(),
                                   selMeshObject.getMaxArea() )
        if dialog.ValueOk:
            selMeshObject.setTag(dialog.tag)
            if dialog.maxArea != 0.0: #using 0.0 to mean no max area (hacky)
                selMeshObject.setMaxArea(dialog.maxArea)
            else:
                selMeshObject.deleteMaxArea()
            userMeshChanged = True
        return userMeshChanged
            
class EditRegionDialog(Dialog):
    """
    Dialog box for editing region info
    """
    # initial values, hard coded.
    # should be values associated with the current mesh


    def __init__(self,
                 parent,
                 tag,
                 maxArea):


        self.tag = tag
        self.maxArea = maxArea

        Dialog.__init__(self, parent)

        
    def body(self, master):
        """
        GUI description
        """
        self.title("Edit Region")
        
        Label(master,
              text='tag:').grid(row=0, sticky=W)
        Label(master,
              text='Maximum Triangle Area').grid(row=1, sticky=W)


        tagVar = StringVar()
        tagVar.set(self.tag)
        self.tagstr   = Entry(master,
                                   width = 16,
                                   textvariable = tagVar,
                                   takefocus=1)

    
        maxAreaVar = StringVar()
        if self.maxArea is None:
          self.maxArea = ""  
        maxAreaVar.set(self.maxArea)
        self.maxAreastr   = Entry(master,
                                   width = 16,
                                   textvariable = maxAreaVar,
                                   takefocus=1)

        self.tagstr.grid(row=0, column=1, sticky=W)
        self.maxAreastr.grid(row=1, column=1, sticky=W)

        self.tag  = -1
        self.maxArea  = -1
        self.ValueOk = False


    def apply(self):
        """
        check entered values
        """
        self.goodMaxArea = self.goodtag = True
        self.ValueOk = True
        try:
            self.tag = self.tagstr.get()
        except ValueError:
            self.ValueOk = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')

        if not self.maxAreastr.get() == "":
            try:    
                self.maxArea = float(self.maxAreastr.get())
                #MeshGenDialog.lastMaxArea =self.maxArea
            except ValueError:
                self.ValueOk = False
                showerror('Bad mesh generation values',
                                   ' Values are not numbers.')
        else:
            self.maxArea = 0.0

        try: 
            # value checking
            if self.goodMaxArea == True and self.maxArea < 0.0:
                raise IOError
            
        except IOError:
            self.ValueOk = False
            showerror('Bad mesh generation values',
                                   'Values are out of range.')

       
class vSegments(vAbstract):
    """
    Define the visualisation of a segment
    """

    def draw(self,v1,v2,mesh,uniqueID,scale,canvas,event):
        """
        draw a segment object, plus add it to the mesh data structure

        event isn't used
        """
        segment = mesh.addUserSegment(v1,v2)
        self.visualise(segment, uniqueID, canvas, scale)
        return segment    
  
    def editWindow_int(self, canvas, selMeshObject):
        ask_int = askinteger("Edit Segment", "Tag",
                             minvalue=0, initialvalue= selMeshObject.tag)
        if ask_int >= 0:
            selMeshObject.tag = ask_int
            
    def editWindow(self, canvas, selMeshObject, userMeshChanged):
        ask_string = askstring("Edit Segment Tag", "Tag",
                              initialvalue= str(selMeshObject.tag))
        if ask_string != None:
            selMeshObject.tag = ask_string
            userMeshChanged = True
        return userMeshChanged
    
    def defaultWindow(self, canvas):
        ask_string = askstring("Edit Default Segment Tag", "Tag",
                     initialvalue= str(self.association.get_default_tag()))
        if ask_string != None:
            self.association.set_default_tag(ask_string)
        
class vTriangles(vAbstract):
    """
    Define the visualisation of a triangle
    """

    def draw(self,mesh,uniqueID,scale,canvas,event):
        """
        currently triangles have no draw function
        """
        raise TypeError


class vMesh(vAbstract):
    """
    Define the visualisation of a collection of mesh object dictionaries (eg vSegments)
    """
    def __init__(self,MeshObjects):
        self.association = mesh.MeshObject
        self.MeshObjectsList = MeshObjects

    def visualise(self, meshobject, uniqueID, canvas, scale):
        """
        This should find the right vAbstract instance and call visualise
        """
        raise Exception

    def unvisualise(self, MeshObject, canvas):
        """
        delete a visual object and
        delete it from the canvas mesh data structure
        """
        for MeshObjects in self.MeshObjectsList:
            if isinstance(MeshObject,MeshObjects.association):
                MeshObjects.unvisualise(MeshObject, canvas)
                break
            
    def hasKey(self, key):
        """
        returns true if the key is in the canvas list
        """
        for MeshObjects in self.MeshObjectsList:
            if key in MeshObjects.oncanvas:
                return True
        return False

    def getMeshObject(self, key):
        """
        returns the mesh object,
        given it's key.

        precon:  The key is valid
        """
        for MeshObjects in self.MeshObjectsList:
            try:
                return MeshObjects.oncanvas[key]
            except KeyError:
                pass
        raise KeyError

    def editWindow(self, canvas, selMeshObject, userMeshChanged):
        for MeshObjects in self.MeshObjectsList:
            if isinstance(selMeshObject,MeshObjects.association):
                userMeshChanged =MeshObjects.editWindow(canvas, selMeshObject,userMeshChanged)
                break
        return userMeshChanged
    
    def defaultWindow(self, canvas, MeshClass):
        for MeshObjects in self.MeshObjectsList:
            if MeshClass == MeshObjects.association:
                MeshObjects.defaultWindow(canvas)
                break 
