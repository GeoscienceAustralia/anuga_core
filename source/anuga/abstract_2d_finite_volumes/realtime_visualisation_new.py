
from visual import *

elevation_color = (0.3,0.4,0.5)
stage_color     = (0.1,0.4,0.99)
friction_color  = (0.1,0.99,0.1)
other_color     = (0.99,0.4,0.1)


class Visualiser:

    def __init__(self,domain,scale_z=1.0,rect=None,title='Test'):
        """Create visualisation of domain
        """

        self.scene = display(title=title)
        self.scene.center = (0.5,0.5,0.0)
        self.scene.forward = vector(-0.0, 0.5, -0.8)
        self.scene.ambient = 0.4
        self.scene.lights = [(0.6, 0.3, 0.2), (0.1, -0.5, 0.4)]

        self.frame  = frame()

        self.domain          = domain
        self.default_scale_z = scale_z
        self.vertices        = domain.vertex_coordinates


        # models for each quantity
        self.vpython_z_models = {}
        #print keys
        for key in self.domain.quantities:
            self.vpython_z_models[key] = faces(frame=self.frame)

        #print self.z_models

        if rect is None:
            self.max_x = max(max(self.vertices[:,0]),max(self.vertices[:,2]),max(self.vertices[:,4]))
            self.min_x = min(min(self.vertices[:,0]),min(self.vertices[:,2]),min(self.vertices[:,4]))
            self.max_y = max(max(self.vertices[:,1]),max(self.vertices[:,3]),max(self.vertices[:,5]))
            self.min_y = min(min(self.vertices[:,1]),min(self.vertices[:,3]),min(self.vertices[:,5]))
        else:
            self.max_x = rect[2]
            self.min_x = rect[0]
            self.max_y = rect[3]
            self.min_y = rect[1]


        self.range_x = self.max_x - self.min_x
        self.range_y = self.max_y - self.min_y
        self.range_xy = max(self.range_x, self.range_y)


#        print 'min_x=',self.min_x
#        print 'max_x=',self.max_x
#        print 'range_x=',self.range_x
#        print 'min_y=',self.min_y
#        print 'max_y=',self.max_y
#        print 'range_y',self.range_y

        self.stage_color     = stage_color
        self.elevation_color = elevation_color
        self.friction_color  = friction_color
        self.other_color     = other_color

        self.qcolor = {}
        self.scale_z = {}
        self.coloring = {}
        self.updating = {}
        self.setup = {}

        self.qcolor['stage']  = self.stage_color
        self.scale_z['stage']  = 1.0
        self.coloring['stage'] = False
        self.updating['stage'] = True
        self.setup['stage'] = True

        self.qcolor['elevation']  = self.elevation_color
        self.scale_z['elevation']  = 1.0
        self.coloring['elevation'] = False
        self.updating['elevation'] = False
        self.setup['elevation'] = False

        #print self.max_x, self.min_x, self.max_y, self.min_y,
        #self.min_bed print self.max_bed, self.range_x, self.range_y,
        #self.range_xy, self.range_bed
        #print 'Calculating triangles'

        self.pos     = zeros( (6*len(domain),3), Float)
        self.colour  = zeros( (6*len(domain),3), Float)
        self.normals = zeros( (6*len(domain),3), Float)

        #print 'keys',self.domain.quantities.keys()
        #print 'shape of stage',shape(self.stage)

        self.border_model = curve(frame = self.frame,
                                  pos=[(0,0),(0,1),(1,1),(1,0),(0,0)])
        self.timer=label(pos=(0.75,0.75,0.5),text='Time=%10.5e'%self.domain.time,
                         visible=False)


        scene.autoscale=0
        #self.update_all()
        #scene.uniform=0
        self.range_z  = 1.0
        if domain.visualise_range_z == None:
            self.setup_range_z()
        else:
            self.range_z = domain.visualise_range_z

        self.max_z = self.range_z/2.0
        self.min_z = -self.range_z/2.0
        self.range_z = max(self.max_z - self.min_z, 1.0e-10)

#        print self.range_z
#        print self.max_z
#        print self.min_z

    def setup_range_z(self,qname1='elevation',qname2='stage'):

        #print qname1
        #print qname2
        range_z = 1.0e-10
        try:
            q1 = self.domain.quantities[qname1].vertex_values
            max_z = max(max(q1))
            min_z = min(min(q1))
            print max_z, min_z
            range_z = max(range_z, max_z - min_z)
        except:
            #print 'Visualisation: could not find range of '+qname1
            pass

        try:
            q2 = self.domain.quantities[qname2].vertex_values
            max_z = max(max(q2))
            min_z = min(min(q2))
            print max_z, min_z
            range_z = max(range_z, max_z - min_z)
        except:
            #print 'Visualisation: could not find range of '+qname2
            pass

        self.range_z = max(range_z,self.range_z)


    def setup_all(self):

        self.scene.select()

        for key in ['elevation','stage']:
            if self.setup[key] | self.updating[key]:
                self.update_quantity(key)



    def update_all(self):

        self.scene.select()

        for key in ['elevation','stage']:
            if self.updating[key]:
                self.update_quantity(key)


    def update_timer(self):

        self.scene.select()

        if self.domain.visualise_timer == True:
            self.timer.visible = True
            self.timer.text='Time=%10.5e'%self.domain.time


    def update_quantity(self,qname):

        #print 'update '+qname+' arrays'
        try:
            qcolor   = self.qcolor[qname]
            scale_z  = self.scale_z[qname]
            coloring = self.coloring[qname]
            updating = self.updating[qname]
        except:
            qcolor   = None
            scale_z  = None
            coloring = False
            updating = True


        if qcolor is None:
            qcolor = self.other_color

            if qname=='elevation':
                qcolor = self.elevation_color

            if qname=='stage':
                qcolor = self.stage_color

            if qname=='friction':
                qcolor = self.friction_color

        if scale_z is None:
            scale_z = self.default_scale_z

        #try:
        if coloring:
            self.update_arrays_color(self.domain.quantities[qname].vertex_values, qcolor, scale_z)
        else:
            self.update_arrays(self.domain.quantities[qname].vertex_values, qcolor, scale_z)

        #print 'update bed image'
        if qname=='stage':
            self.pos[:,2] = self.pos[:,2] - self.domain.visualise_wet_dry_cutoff

        self.vpython_z_models[qname].pos    = self.pos
        self.vpython_z_models[qname].color  = self.colour
        self.vpython_z_models[qname].normal = self.normals
        #except:
        #    print 'Visualisation: Could not update quantity '+qname
        #    pass

    def update_quantity_color(self,qname,qcolor=None,scale_z=None):

        #print 'update '+qname+' arrays'

        if qcolor is None:
            if qname=='elevation':
                qcolor = self.elevation_color

            if qname=='stage':
                qcolor = self.stage_color

            if qname=='friction':
                qcolor = self.friction_color

        if scale_z is None:
            scale_z = self.scale_z

        #try:
        self.update_arrays_color(self.domain.quantities[qname].vertex_values, qcolor, scale_z)

            #print 'update bed image'
        self.vpython_z_models[qname].pos    = self.pos
        self.vpython_z_models[qname].color  = self.colour
        self.vpython_z_models[qname].normal = self.normals
        #except:
        #    print 'Visualisation: Could not update quantity '+qname
        #    pass


    def update_arrays(self,quantity,qcolor,scale_z):

        col = asarray(qcolor)

        N = len(self.domain)

        min_x = self.min_x
        min_y = self.min_y
        range_xy = self.range_xy
        range_z = self.range_z

        vertices = self.vertices
        pos      = self.pos
        normals  = self.normals
        colour   = self.colour

        try:
            update_arrays_weave(vertices,quantity,col,pos,normals,colour,N,
                                min_x,min_y,range_xy,range_z,scale_z)
        except:
            update_arrays_python(vertices,quantity,col,pos,normals,colour,N,
                                 min_x,min_y,range_xy,range_z,scale_z)



    def update_arrays_color(self,quantity,qcolor,scale_z):

        col = asarray(qcolor)

        N = len(self.domain)

        min_x = self.min_x
        min_y = self.min_y
        range_xy = self.range_xy
        range_z = self.range_z

        vertices = self.vertices
        pos      = self.pos
        normals  = self.normals
        colour   = self.colour

        #try:
        update_arrays_color_python(vertices,quantity,col,pos,normals,colour,N,
                                min_x,min_y,range_xy,range_z,scale_z)
        #except:
        #    update_arrays_color_python(vertices,quantity,col,pos,normals,colour,N,
        #                         min_x,min_y,range_xy,range_z,scale_z)


#==================================================================================

def  update_arrays_python(vertices,quantity,col,pos,normals,colour,N,
                          min_x,min_y,range_xy,range_z,scale_z):

    from math import sqrt
    normal = zeros( 3, Float)
    v = zeros( (3,3), Float)
    s = 1.0

    for i in range( N ):

        for j in range(3):
            v[j,0] = (vertices[i,2*j  ]-min_x)/range_xy
            v[j,1] = (vertices[i,2*j+1]-min_y)/range_xy
            v[j,2] = quantity[i,j]/range_z*scale_z*0.5

        v10 = v[1,:]-v[0,:]
        v20 = v[2,:]-v[0,:]

        normal[0] = v10[1]*v20[2] - v20[1]*v10[2]
        normal[1] = v10[2]*v20[0] - v20[2]*v10[0]
        normal[2] = v10[0]*v20[1] - v20[0]*v10[1]

        norm = sqrt( normal[0]**2 + normal[1]**2 + normal[2]**2)

        normal[0] = normal[0]/norm
        normal[1] = normal[1]/norm
        normal[2] = normal[2]/norm

        pos[6*i  ,:] = v[0,:]
        pos[6*i+1,:] = v[1,:]
        pos[6*i+2,:] = v[2,:]
        pos[6*i+3,:] = v[0,:]
        pos[6*i+4,:] = v[2,:]
        pos[6*i+5,:] = v[1,:]



        colour[6*i  ,:] = s*col
        colour[6*i+1,:] = s*col
        colour[6*i+2,:] = s*col
        colour[6*i+3,:] = s*col
        colour[6*i+4,:] = s*col
        colour[6*i+5,:] = s*col

        s =  15.0/8.0 - s

        normals[6*i  ,:] = normal
        normals[6*i+1,:] = normal
        normals[6*i+2,:] = normal
        normals[6*i+3,:] = -normal
        normals[6*i+4,:] = -normal
        normals[6*i+5,:] = -normal



def  update_arrays_weave(vertices,quantity,col,pos,normals,colour,
                         N,min_x,min_y,range_xy,range_z,scale_z):

    import weave
    from weave import converters

    from math import sqrt
    normal = zeros( 3, Float)
    v = zeros( (3,3), Float)
    v10 = zeros( 3, Float)
    v20 = zeros( 3, Float)

    code1 = """

        double s = 1.0;

        for (int i=0; i<N ; i++ ) {
            for (int j=0; j<3 ; j++) {
                v(j,0) = (vertices(i,2*j  )-min_x)/range_xy;
                v(j,1) = (vertices(i,2*j+1)-min_y)/range_xy;
                v(j,2) = quantity(i,j)/range_z*scale_z*0.5;
            }


            for (int j=0; j<3; j++) {
                v10(j) = v(1,j)-v(0,j);
                v20(j) = v(2,j)-v(0,j);
            }

            normal(0) = v10(1)*v20(2) - v20(1)*v10(2);
            normal(1) = v10(2)*v20(0) - v20(2)*v10(0);
            normal(2) = v10(0)*v20(1) - v20(0)*v10(1);

            double norm =  sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));

            normal(0) = normal(0)/norm;
            normal(1) = normal(1)/norm;
            normal(2) = normal(2)/norm;


            for (int j=0; j<3; j++) {
                pos(6*i  ,j) = v(0,j);
                pos(6*i+1,j) = v(1,j);
                pos(6*i+2,j) = v(2,j);
                pos(6*i+3,j) = v(0,j);
                pos(6*i+4,j) = v(2,j);
                pos(6*i+5,j) = v(1,j);

                colour(6*i  ,j) = s*col(j);
                colour(6*i+1,j) = s*col(j);
                colour(6*i+2,j) = s*col(j);
                colour(6*i+3,j) = s*col(j);
                colour(6*i+4,j) = s*col(j);
                colour(6*i+5,j) = s*col(j);

                normals(6*i  ,j) = normal(j);
                normals(6*i+1,j) = normal(j);
                normals(6*i+2,j) = normal(j);
                normals(6*i+3,j) = -normal(j);
                normals(6*i+4,j) = -normal(j);
                normals(6*i+5,j) = -normal(j);
                }

            s =  15.0/8.0 - s;
        }

        """

    weave.inline(code1, ['vertices','quantity','col','pos','normals','colour','N',
                         'min_x','min_y','range_xy','range_z','scale_z','v','normal','v10','v20'],
                 type_converters = converters.blitz, compiler='gcc');

def  update_arrays_color_python(vertices,quantity,col,pos,normals,colour,N,
                          min_x,min_y,range_xy,range_z,scale_z):

    from math import sqrt
    normal = zeros( 3, Float)
    v = zeros( (3,3), Float)
    s = 1.0

    for i in range( N ):

        for j in range(3):
            v[j,0] = (vertices[i,2*j  ]-min_x)/range_xy
            v[j,1] = (vertices[i,2*j+1]-min_y)/range_xy
            v[j,2] = quantity[i,j]/range_z*scale_z*0.5

        v10 = v[1,:]-v[0,:]
        v20 = v[2,:]-v[0,:]

        normal[0] = v10[1]*v20[2] - v20[1]*v10[2]
        normal[1] = v10[2]*v20[0] - v20[2]*v10[0]
        normal[2] = v10[0]*v20[1] - v20[0]*v10[1]

        norm = sqrt( normal[0]**2 + normal[1]**2 + normal[2]**2)

        normal[0] = normal[0]/norm
        normal[1] = normal[1]/norm
        normal[2] = normal[2]/norm

        pos[6*i  ,:] = v[0,:]
        pos[6*i+1,:] = v[1,:]
        pos[6*i+2,:] = v[2,:]
        pos[6*i+3,:] = v[0,:]
        pos[6*i+4,:] = v[2,:]
        pos[6*i+5,:] = v[1,:]

        q0 = quantity[i,0]/0.4
        q1 = quantity[i,1]/0.4
        q2 = quantity[i,2]/0.4

        q0r = min(q0,0.0)
        q1r = min(q1,0.0)
        q2r = min(q2,0.0)
        q0b = max(q0,0.0)
        q1b = max(q1,0.0)
        q2b = max(q2,0.0)


#        colour[6*i  ,:] = s*array([0.3 - q0r, 0.3 - q0, 0.3 + q0b])
#        colour[6*i+1,:] = s*array([0.3 - q1r, 0.3 - q1, 0.3 + q1b])
#        colour[6*i+2,:] = s*array([0.3 - q2r, 0.3 - q2, 0.3 + q2b])
#        colour[6*i+3,:] = s*array([0.3 - q0r, 0.3 - q0, 0.3 + q0b])
#        colour[6*i+4,:] = s*array([0.3 - q2r, 0.3 - q2, 0.3 + q2b])
#        colour[6*i+5,:] = s*array([0.3 - q1r, 0.3 - q1, 0.3 + q1b])

        s = 0.1
        colour[6*i  ,:] = s*array([0.3 - q0r, 0.3, 0.3 + q0b])
        colour[6*i+1,:] = s*array([0.3 - q1r, 0.3, 0.3 + q1b])
        colour[6*i+2,:] = s*array([0.3 - q2r, 0.3, 0.3 + q2b])
        colour[6*i+3,:] = s*array([0.3 - q0r, 0.3, 0.3 + q0b])
        colour[6*i+4,:] = s*array([0.3 - q2r, 0.3, 0.3 + q2b])
        colour[6*i+5,:] = s*array([0.3 - q1r, 0.3, 0.3 + q1b])


        s =  15.0/8.0 - s

        normals[6*i  ,:] = normal
        normals[6*i+1,:] = normal
        normals[6*i+2,:] = normal
        normals[6*i+3,:] = -normal
        normals[6*i+4,:] = -normal
        normals[6*i+5,:] = -normal



def  update_arrays_color_weave(vertices,quantity,col,pos,normals,colour,
                         N,min_x,min_y,range_xy,range_z,scale_z):

    import weave
    from weave import converters

    from math import sqrt
    normal = zeros( 3, Float)
    v = zeros( (3,3), Float)
    v10 = zeros( 3, Float)
    v20 = zeros( 3, Float)

    code1 = """
        double s = 1.0;

        for (int i=0; i<N ; i++ ) {
            for (int j=0; j<3 ; j++) {
                v(j,0) = (vertices(i,2*j  )-min_x)/range_xy;
                v(j,1) = (vertices(i,2*j+1)-min_y)/range_xy;
                v(j,2) = quantity(i,j)/range_z*scale_z*0.5;
            }


            for (int j=0; j<3; j++) {
                v10(j) = v(1,j)-v(0,j);
                v20(j) = v(2,j)-v(0,j);
            }

            normal(0) = v10(1)*v20(2) - v20(1)*v10(2);
            normal(1) = v10(2)*v20(0) - v20(2)*v10(0);
            normal(2) = v10(0)*v20(1) - v20(0)*v10(1);

            double norm =  sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));

            normal(0) = normal(0)/norm;
            normal(1) = normal(1)/norm;
            normal(2) = normal(2)/norm;



            for (int j=0; j<3; j++) {
                pos(6*i  ,j) = v(0,j);
                pos(6*i+1,j) = v(1,j);
                pos(6*i+2,j) = v(2,j);
                pos(6*i+3,j) = v(0,j);
                pos(6*i+4,j) = v(2,j);
                pos(6*i+5,j) = v(1,j);

                normals(6*i  ,j) = normal(j);
                normals(6*i+1,j) = normal(j);
                normals(6*i+2,j) = normal(j);
                normals(6*i+3,j) = -normal(j);
                normals(6*i+4,j) = -normal(j);
                normals(6*i+5,j) = -normal(j);
                }

            s = 0.2+fabs((v(0,2)+v(1,2)+v(2,2))/3.0);
            for (int j=0; j<3; j++) {
                colour(6*i  ,j) = s*col(j);
                colour(6*i+1,j) = s*col(j);
                colour(6*i+2,j) = s*col(j);
                colour(6*i+3,j) = s*col(j);
                colour(6*i+4,j) = s*col(j);
                colour(6*i+5,j) = s*col(j);
                }
            s =  15.0/8.0 - s;
        }

        """

    weave.inline(code1, ['vertices','quantity','col','pos','normals','colour','N',
                         'min_x','min_y','range_xy','range_z','scale_z','v','normal','v10','v20'],
                 type_converters = converters.blitz, compiler='gcc');


def setup_scene():

    #scene.width = 1000
    #scene.height = 800

    #Original
    scene.center = (0.5,0.5,0.0)
    scene.forward = vector(-0.0, 0.5, -0.8)

    #Temporary (for bedslope)
    #scene.forward = vector(0.0006, 0.7, -0.03)


    #Temporary for hackett - begin
    #scene.autoscale = 0
    #scene.scale = (0.002, 0.002, 0.01) #Scale z so that countours stand out more
    #scene.center = (300.0,500.0,-10)
    #Temporary for hackett - end


    scene.ambient = 0.4
    #scene.lights = [(0.6, 0.3, 0.2), (0.1, -0.5, 0.4), (-0.1, 0.1, -0.4),
    #               (-0.2, 0.2, 0.1)]

    scene.lights = [(0.6, 0.3, 0.2), (0.1, -0.5, 0.4)]



def create_surface(domain,scale_z=1.0):

    domain.initialise_visualiser(scale_z)

    domain.visualiser.update_quantity('elevation')
    domain.visualiser.update_quantity('stage')
    domain.visualiser.update_timer()

def update(domain):

    if domain.visualise_color_stage:
        domain.visualiser.update_quantity_color('stage')
    else:
        domain.visualiser.update_quantity('stage')

    domain.visualiser.update_timer()



if __name__ == '__main__':

    from advection import Domain as A_Domain
    from shallow_water import Domain as SW_Domain
    from shallow_water import Constant_height, Constant_stage
    from mesh_factory import rectangular

    points, vertices, boundary = rectangular(60, 60, len1=200, len2=200)
    a_domain  = A_Domain(points, vertices, boundary)
    print 'No of Elements' , a_domain.number_of_elements
    #sw_domain = SW_Domain(points, vertices, boundary)
    #sw_domain.set_quantity('elevation', Constant_stage(0.0))
    #sw_domain.set_quantity('stage', Constant_stage(0.0))

    a_domain.initialise_visualiser()
    a_domain.visualiser.setup_all()


    #sw_domain.initialise_visualiser()
    #sw_domain.visualiser.setup['elevation']=True
    #sw_domain.visualiser.setup_all()

    time = 0.0
    while (True):
        time = time+0.001
        #print 'Time = ',time
        #sw_domain.set_quantity('stage', time)
        a_domain.set_quantity('stage',time)

        #sw_domain.visualiser.update_quantity('stage')
        a_domain.visualiser.update_quantity('stage')

        if time>1.0 :
            break

    #a_v = Visualiser(domain=a_domain, title='advection')
    #a_v.update_all()


    #sw_v = Visualiser(domain=sw_domain, title='shallow water')
    #sw_v.update_all()
