
SET_COLOUR = 'red'
DEFAULT_ATTRIBUTE = 'elevation'



    def Out_createSetTools(self):
        """
        Add set tool buttons to the top of the GUI
        """
        ToolBarButton(self, self.toolbar, 'sep', 'sep.gif', width=10,
                      state='disabled')
        for key, func, balloon in [
                ('threshold', self.threshold, 'threshold the set'),
                ('Courant_threshold', self.Courant_threshold, 'Courant_threshold the set'),
                ('gradient_threshold', self.gradient_threshold, 'gradient_threshold the set'),
                ('smooth', self.smooth_polySet, 'smooth the polygons'),
                ('polyset', self.triangles_to_polySet, 'make a poly set out of selected triangles')]:      
                #('refineSet', self.refineSet, 'Refine the set')]:
            ToolBarButton(self, self.toolbar, key, '%s.gif' %key,
                          command=func, balloonhelp=balloon,
                               statushelp='' )

            

    def createSetIcons(self):
        """
        Add Edit buttons to the top of the GUI
        """
        ToolBarButton(self, self.toolbar, 'sep', 'sep.gif', width=10,
                      state='disabled')
        for key, func, balloon in [
                ('selectAllTriangles', self.selectAllTriangles, 'select all'),
                ('none', self.clearSelection, 'clear selection')]:
            ToolBarButton(self, self.toolbar, key, '%s.gif' %key,
                          command=func, balloonhelp=balloon,
                               statushelp='' )

    def refineSet(self,parent):
        self.mesh.refineSet(self.selSet)
        self.visualiseMesh(self.mesh)


    def setStructureNumber(self,parent):
        dialog =  setStructureNumberDialog(self.canvas)
        if dialog.numberOK:
            self.structureSize = dialog.number

    def erode(self, parent):
#Not implimented
        self.canvas.delete(ALL)
        self.mesh.erode(self.selSet,structureSize=self.structureSize)
        self.visualiseMesh(self.mesh)


    def dilate(self, parent):
#Not implimented
        self.canvas.delete(ALL)
        self.mesh.dilate(self.selSet,structureSize=self.structureSize)
        self.visualiseMesh(self.mesh)

    def general_threshold(self,parent,function,function_description):
        """
        add a vertex using a window and entering x y values.

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        if self.selSet == 'None':
            self.selectAllTriangles(parent)

        dialog = GeneralThresholdDialog(self.canvas,self.mesh.attributeTitles,function_description)
        if dialog.minmaxValuesOk:
            self.canvas.delete(ALL)
            min = dialog.min
            max = dialog.max
            attribute_name = dialog.attribute_name
            self.mesh.general_threshold(self.selSet,min=min,max=max,attribute_name = attribute_name,function=function)
            self.visualiseMesh(self.mesh)

    def threshold(self, parent):
        """
        add a vertex using a window and entering x y values.

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        function = self.mesh.av_att
        function_description = 'average attribute of triangle'
        self.general_threshold(parent,function,function_description)


    def Courant_threshold(self, parent):
        """
        add a vertex using a window and entering x y values.

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        function = self.mesh.Courant_ratio
        function_description = 'average attribute/area of triangle'
        self.general_threshold(parent,function,function_description)

    def gradient_threshold(self, parent):
        """
        add a vertex using a window and entering x y values.

        the parent attribute isn't used by this function.
        need to userstand toolbarbutton.py to know how to
        get rid of it.
        """
        function = self.mesh.Gradient
        function_description = 'average gradient of triangle'
        self.general_threshold(parent,function,function_description)

    def smooth_polySet(self,parent):
        dialog = SmoothDialog(self.canvas)
        if dialog.valueOK:
            min_radius = dialog.min_radius
            self._smooth_polySet(min_radius)

    def _smooth_polySet(self,min_radius):
        userVertices,userSegments,alphaSegments = \
            self.mesh.smooth_polySet(min_radius=min_radius)

        self.mesh.userVertices=[]
        self.mesh.userSegments=[]
        self.mesh.alphaSegments=[]
        self.canvas.delete(ALL)
        event = None
        print 'len(userVertices.keys())'
        print len(userVertices.keys())
        print 'len(userSegments.keys())'
        print len(userSegments.keys())
        print 'len(alphaSegments.keys())'
        print len(alphaSegments.keys())

        #######
        point_keys = {}
        for vert in userVertices.keys():
            assert not point_keys.has_key((vert.x,vert.y))
            point_keys[(vert.x,vert.y)]=vert
        assert len(point_keys.keys())==len(userVertices.keys())
        #######

        for v in userVertices.keys():
            x = v.x*self.SCALE
            y = v.y*self.SCALE
            userVertices[(v.x,v.y)]=self.drawVertex(x,y,event)

        for line in userSegments.keys():
            v0 = userVertices[line[0]]
            v1 = userVertices[line[1]]
            segment = self.drawSegment(v0,v1)
            segment.set_tag(userSegments[line].tag)

        for line in alphaSegments.keys():
            v0 = userVertices[line[0]]
            v1 = userVertices[line[1]]
            segment = self.drawSegment(v0,v1)
            segment.set_tag(alphaSegments[line].tag)
        self.visualiseMesh(self.mesh)


    def triangles_to_polySet(self,parent):
        userVertices,userSegments,alphaSegments = \
            self.mesh.triangles_to_polySet(self.selSet)

        self.mesh.userVertices=[]
        self.canvas.delete(ALL)

        event = None
        print 'len(userVertices.keys())'
        print len(userVertices.keys())
        print 'len(userSegments.keys())'
        print len(userSegments.keys())
        print 'len(alphaSegments.keys())'
        print len(alphaSegments.keys())


        #######
        point_keys = {}
        for vert in userVertices.keys():
            assert not point_keys.has_key((vert.x,vert.y))
            point_keys[(vert.x,vert.y)]=vert
        assert len(point_keys.keys())==len(userVertices.keys())
        #######

        for v in userVertices.keys():
            if userVertices[v] is True:
                x = v.x*self.SCALE
                y = v.y*self.SCALE
                userVertices[(v.x,v.y)]=self.drawVertex(x,y,event)

        for line in userSegments.keys():
            v0 = userVertices[line[0]]
            v1 = userVertices[line[1]]
            segment = self.drawSegment(v0,v1)
            segment.set_tag(userSegments[line].tag)

        for line in alphaSegments.keys():
            v0 = userVertices[line[0]]
            v1 = userVertices[line[1]]
            segment = self.drawSegment(v0,v1)
            segment.set_tag(alphaSegments[line].tag)
        self.visualiseMesh(self.mesh)
        #self.smooth_polySet(parent)

    def selectTriangles(self,setName):
        """
        """
        self.canvas.delete(ALL)
        self.selSet = setName
        self.visualiseMesh(self.mesh)

    def selectAllTriangles(self,parent):
        """
        selected all triangles in the mesh
        """
        self.canvas.delete(ALL)
        self.selSet = self.mesh.selectAllTriangles()
        self.visualiseMesh(self.mesh)
     
class  GeneralThresholdDialog(Dialog):
    """
    Dialog box for thresholding a set by entering minimum 
    and maximum values
    """
    def __init__(self,
                 parent,
                 attribute_titles,
                 function_description):
        self.attribute_titles=attribute_titles
        self.function_description=function_description

        Dialog.__init__(self, parent)


    def body(self, master):
        """
        GUI description
        """
        self.title("Threshold selected set")
        blurb1 = 'Threshold selected set between minimum'
        blurb2 = 'and maximum ' + self.function_description

        Label(master,text=blurb1).grid(row=0, sticky=W)
        Label(master,text=blurb2).grid(row=1, sticky=W)

        Label(master, text='minimum attribute:').grid(row=2, sticky=W)
        Label(master, text='maximum attribute:').grid(row=3, sticky=W)
        Label(master, text='attribute name').grid(row=4, sticky=W)


        nameVar = StringVar()
        nameVar.set('elevation')

        self.minstr = Entry(master, width = 16, name ="entry")
        self.maxstr = Entry(master, width = 16)
        self.attstr = Entry(master, width = 16,textvariable = nameVar)
       
        self.minstr.grid(row=2, column=1, sticky=W)
        self.maxstr.grid(row=3, column=1, sticky=W)
        self.attstr.grid(row=4, column=1, sticky=W)
        self.minstr.focus_force()
        self.min  = 0
        self.max  = 0
        self.attribute_name = 'elevation'
        self.minmaxValuesOk = False
        
    def apply(self):
        self.minmaxValuesOk = True
        try:
            self.min = float(self.minstr.get())
            self.max = float(self.maxstr.get())
        except ValueError:
            self.minmaxValuesOk = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')
        try:
            self.attribute_titles.index(self.attstr.get())#dodgey.
            self.attribute_name = self.attstr.get()
        except ValueError:
            self.attribute_name = None
            showerror('Bad attribute name',
                                   'Using h = 1')

class  ThresholdDialog(Dialog):
    """
    Dialog box for thresholding a set by entering minimum 
    and maximum values
    """
    def __init__(self,
                 parent,
                 attribute_titles):
        self.attribute_titles=attribute_titles
        Dialog.__init__(self, parent)


    def body(self, master):
        """
        GUI description
        """
        self.title("Threshold selected set")
        
        Label(master, text='minimum attribute:').grid(row=0, sticky=W)
        Label(master, text='maximum attribute:').grid(row=1, sticky=W)
        Label(master, text='attribute name').grid(row=2, sticky=W)


        nameVar = StringVar()
        nameVar.set('elevation')

        self.minstr   = Entry(master, width = 16, name ="entry")
        self.maxstr   = Entry(master, width = 16)
        self.attstr   = Entry(master, width = 16,textvariable = nameVar)
       
        self.minstr.grid(row=0, column=1, sticky=W)
        self.maxstr.grid(row=1, column=1, sticky=W)
        self.attstr.grid(row=2, column=1, sticky=W)
        self.minstr.focus_force()
        self.min  = 0
        self.max  = 0
        self.attribute_name = 'elevation'
        self.minmaxValuesOk = False
        
    def apply(self):
        self.minmaxValuesOk = True
        try:
            self.min = float(self.minstr.get())
            self.max = float(self.maxstr.get())
        except ValueError:
            self.minmaxValuesOk = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')
        try:
            self.attribute_titles.index(self.attstr.get())#dodgey.
            self.attribute_name = self.attstr.get()
        except ValueError:
            self.minmaxValuesOk = False
            showerror('Bad attribute name',
                                   ' Attribute not in mesh.')


class  Courant_ThresholdDialog(Dialog):
    """
    Dialog box for thresholding a set by entering minimum 
    and maximum values
    """
    def __init__(self,
                 parent,
                 attribute_titles):
        self.attribute_titles=attribute_titles
        Dialog.__init__(self, parent)


    def body(self, master):
        """
        GUI description
        """
        self.title("Courant_Threshold selected set")
        
        Label(master, text='minimum attribute:').grid(row=0, sticky=W)
        Label(master, text='maximum attribute:').grid(row=1, sticky=W)
        Label(master, text='attribute name').grid(row=2, sticky=W)


        nameVar = StringVar()
        nameVar.set('elevation')

        self.minstr   = Entry(master, width = 16, name ="entry")
        self.maxstr   = Entry(master, width = 16)
        self.attstr   = Entry(master, width = 16,textvariable = nameVar)
       
        self.minstr.grid(row=0, column=1, sticky=W)
        self.maxstr.grid(row=1, column=1, sticky=W)
        self.attstr.grid(row=2, column=1, sticky=W)
        self.minstr.focus_force()
        self.min  = 0
        self.max  = 0
        self.attribute_name = 'elevation'
        self.minmaxValuesOk = False
        
    def apply(self):
        self.minmaxValuesOk = True
        try:
            self.min = float(self.minstr.get())
            self.max = float(self.maxstr.get())
        except ValueError:
            self.minmaxValuesOk = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')
        try:
            self.attribute_titles.index(self.attstr.get())#dodgey.
            self.attribute_name = self.attstr.get()
        except ValueError:
            self.minmaxValuesOk = False
            showerror('Bad attribute name',
                                   ' Attribute not in mesh.')


class SmoothDialog(Dialog):
    """
    Dialog box for setting the number of triangles
    used to make up dilation or erosion
    """
    def body(self, master):
        """
        GUI description
        """
        self.title("Enter the minimum radius to remove")
        
        Label(master, text='radius:').grid(row=0, sticky=W)

        self.min_radius = Entry(master, width = 16, name ="entry")
       
        self.min_radius.grid(row=0, column=1, sticky=W)
        self.min_radius.focus_force()
        self.min = 2.
        self.valueOK = False
        
    def apply(self):
        self.valueOK = True
        try:
            self.min = float(self.min_radius.get())
            self.min_radius = self.min
        except ValueError:
            self.valueOK = False
            showerror('Bad Number',
                                   ' Value not a number')

class  setStructureNumberDialog(Dialog):
    """
    Dialog box for setting the number of triangles
    used to make up dilation or erosion
    """
    def body(self, master):
        """
        GUI description
        """
        self.title("Set number of elements effected by morphing sets")
        
        Label(master, text='number:').grid(row=0, sticky=W)

        self.number = Entry(master, width = 16, name ="entry")
       
        self.number.grid(row=0, column=1, sticky=W)
        self.number.focus_force()
        self.number = 0
        self.numberOk = False
        
    def apply(self):
        self.numberOk = True
        try:
            self.number = int(self.number.get())
        except ValueError:
            self.numberOK = False
            showerror('Bad mesh generation values',
                                   ' Values are not numbers.')


      
