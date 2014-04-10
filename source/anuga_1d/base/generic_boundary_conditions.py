"""boundary.py - Classes for implementing boundary conditions
"""

import numpy

class Boundary:
    """Base class for boundary conditions.
       Specific boundary conditions must provide values for
       the conserved_quantities

       A boundary object has one neighbour; the one it
       serves.
    """

    def __init__(self):
        pass

    def evaluate(self, vol_id=None, edge_id=None):
        msg = 'Generic class Boundary must be subclassed'
        raise msg


class Transmissive_boundary(Boundary):
    """Transmissive boundary returns same conserved quantities as
    those present in its neighbour volume.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for transmissive boundary'
            raise msg

        self.domain = domain

    def __repr__(self):
        return 'Transmissive_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries return the edge values
        of the volume they serve.
        """

        q = self.domain.get_conserved_quantities(vol_id, edge = edge_id)
        return q


class Dirichlet_boundary(Boundary):
    """Dirichlet boundary returns constant values for the
    conserved quantities
    """


    def __init__(self, conserved_quantities=None):
        Boundary.__init__(self)

        if conserved_quantities is None:
            msg = 'Must specify one value for each conserved quantity'
            raise msg

        self.conserved_quantities=numpy.array(conserved_quantities, numpy.float)

    def __repr__(self):
        return 'Dirichlet boundary (%s)' %self.conserved_quantities

    def evaluate(self, vol_id=None, edge_id=None):
        return self.conserved_quantities



class Time_boundary(Boundary):
    """Time dependent boundary returns values for the
    conserved quantities as a function of time.
    Must specify domain to get access to model time and a function f(t)
    which must return conserved quantities as a function time
    """


    def __init__(self, domain = None, f = None):
        Boundary.__init__(self)

        try:
            q = f(0.0)
        except Exception, e:
            msg = 'Function for time boundary could not be executed:\n%s' %e
            raise msg

        try:
            q = numpy.array(q, numpy.float)
        except:
            msg = 'Return value from time boundary function could '
            msg += 'not be converted into a numpy array of numpy.floats.\n'
            msg += 'Specified function should return either list or array.'
            raise msg

        msg = 'ERROR: Time boundary function must return a 1d list or array '
        assert len(q.shape) == 1, msg

        d = len(domain.conserved_quantities)
        msg = 'Return value for f must be a list or an array of length %d' %d
        assert len(q) == d, msg

        self.f = f
        self.domain = domain

    def __repr__(self):
        return 'Time boundary'

    def evaluate(self, vol_id=None, edge_id=None):
        return self.f(self.domain.time)


class File_boundary_time(Boundary):
    """Boundary values obtained from file and interpolated.
    conserved quantities as a function of time.

    Assumes that file contains a time series.

    No spatial info assumed.
    """

    #FIXME: Is this still necessary

    def __init__(self, filename, domain):
        import time
        from config import time_format
        from util import File_function

        Boundary.__init__(self)

        self.F = File_function(filename, domain)
        self.domain = domain

        #Test
        q = self.F(0)

        d = len(domain.conserved_quantities)
        msg = 'Values specified in file must be a list or an array of length %d' %d
        assert len(q) == d, msg


    def __repr__(self):
        return 'File boundary'

    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain.time

        vol_id and edge_id are ignored
        """

        t = self.domain.time
        return self.F(t)




class File_boundary(Boundary):
    """Boundary values obtained from file and interpolated.
    conserved quantities as a function of time.

    Assumes that file contains a time series and possibly
    also spatial info.
    See docstring for File_function in util.py for details about
    admissible file formats

    The full solution history is not exactly the same as if
    the models were coupled:
    File boundary must read and interpolate from *smoothed* version
    as stored in sww and cannot work with the discontinuos triangles.

    """

    def __init__(self, filename, domain, verbose = False):
        import time
        from config import time_format
        from util import file_function

        Boundary.__init__(self)

        #Get x,y vertex coordinates for all triangles
        V = domain.vertex_coordinates

        #Compute midpoint coordinates for all boundary elements
        #Only a subset may be invoked when boundary is evaluated but
        #we don't know which ones at this stage since this object can be attached to
        #any tagged boundary later on.

        if verbose: print 'Find midpoint coordinates of entire boundary'
        self.midpoint_coordinates = numpy.zeros( (len(domain.boundary), 2), numpy.float)
        boundary_keys = domain.boundary.keys()


        xllcorner = domain.geo_reference.get_xllcorner()
        yllcorner = domain.geo_reference.get_yllcorner()        
        

        #Make ordering unique #FIXME: should this happen in domain.py?
        boundary_keys.sort()


        #Record ordering #FIXME: should this also happen in domain.py?
        self.boundary_indices = {}
        for i, (vol_id, edge_id) in enumerate(boundary_keys):

            x0 = V[vol_id, 0]; y0 = V[vol_id, 1]
            x1 = V[vol_id, 2]; y1 = V[vol_id, 3]
            x2 = V[vol_id, 4]; y2 = V[vol_id, 5]

            #Compute midpoints
            if edge_id == 0: m = array([(x1 + x2)/2, (y1 + y2)/2])
            if edge_id == 1: m = array([(x0 + x2)/2, (y0 + y2)/2])
            if edge_id == 2: m = array([(x1 + x0)/2, (y1 + y0)/2])

            #Convert to absolute UTM coordinates
            m[0] += xllcorner
            m[1] += yllcorner
            
            #Register point and index
            self.midpoint_coordinates[i,:] = m
            self.boundary_indices[(vol_id, edge_id)] = i


        if verbose: print 'Initialise file_function'
        self.F = file_function(filename, domain,
	                       interpolation_points=self.midpoint_coordinates,
                               verbose = verbose)
        self.domain = domain

        #Test
        q = self.F(0, point_id=0)

        d = len(domain.conserved_quantities)
        msg = 'Values specified in file %s must be ' %filename
        msg += ' a list or an array of length %d' %d
        assert len(q) == d, msg


    def __repr__(self):
        return 'File boundary'


    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain.time

        vol_id and edge_id are ignored
        """

        t = self.domain.time

        if vol_id is not None and edge_id is not None:
            i = self.boundary_indices[ vol_id, edge_id ]
            return self.F(t, point_id = i)
        else:
            #raise 'Boundary call without point_id not implemented'
            #FIXME: What should the semantics be?
            return self.F(t)








#THIS FAR (10/8/4)
class Connective_boundary(Boundary):
    """Connective boundary returns values for the
    conserved quantities from a volume as defined by a connection table
    mapping between tuples of (volume id, face id) for volumes that
    have their boundaries connected.

    FIXME: Perhaps include possibility of mapping between
    different domains as well

    FIXME: In case of shallow water we may want to have a
    special version that casts this in terms of height rather than stage
    """


    def __init__(self, table):
        from domain import Volume

        Boundary.__init__(self)

        self.connection_table = table
        self.Volume = Volume

    def __repr__(self):
        return 'Connective boundary'

    #FIXME: IF we ever need to get field_values from connected volume,
    #that method could be overridden here (using same idea as in
    #get_conserved_quantities
    #def get_field_values()

    def get_conserved_quantities(self, volume, face=0):

        id = volume.id
        if self.connection_table.has_key((id, face)):
            other_id, other_face = self.connection_table[(id, face)]

            other_volume = self.Volume.instances[other_id]
            cmd = 'q = other_volume.conserved_quantities_face%d' %face;
            exec(cmd)
            return q
        else:
            msg = 'Volume, face tuple ($d, %d) has not been mapped'\
                  %(id, face)
            raise msg





#FIXME: Add a boundary with a general function of x,y, and t

#FIXME: Add periodic boundaries e.g.:
# Attempt at periodic conditions from advection_spik. Remember this
#
#first = 2*(N-1)*N
#for i in range(1,2*N+1,2):
#    k = first + i-1#
#
#    print i,k
#
#    domain[i].faces[2].neighbour = domain[k].faces[1]
#    domain[k].faces[1].neighbour = domain[i].faces[2]



class General_boundary(Boundary):
    """General boundary which can compute conserved quantities based on
    their previous value, conserved quantities of its neighbour and model time.

    Must specify initial conserved quantities,
    neighbour,
    domain to get access to model time
    a function f(q_old, neighbours_q, t) which must return
    new conserved quantities q as a function time

    FIXME: COMPLETE UNTESTED - JUST AN IDEA
    """

    def __init__(self, neighbour=None, conserved_quantities=None, domain=None, f=None):
        Boundary.__init__(self, neighbour=neighbour, conserved_quantities=conserved_quantities)

        self.f = f
        self.domain = domain


    def get_conserved_quantities(self, volume=None, face=0):
        return self.f(self.conserved_quantities,
                      neighbour.conserved_quantities,
                      self.domain.time)





class File_boundary_old(Boundary):
    """Boundary values obtained from file and interpolated.
    conserved quantities as a function of time.

    Assumes that file contains a time series.

    No spatial info assumed.
    """


    def __init__(self, domain = None, filename = None):
        import time
        from config import time_format

        Boundary.__init__(self)


        try:
            fid = open(filename)
        except Exception, e:
            msg = 'Boundary file %s could not be opened: %s\n' %(filename, e)
            raise msg


        line = fid.readline()
        fid.close()
        fields = line.split(',')
        msg = 'File %s must have the format date, values'
        assert len(fields) == 2, msg

        try:
            starttime = time.mktime(time.strptime(fields[0], time_format))
        except ValueError:
            msg = 'First field in file %s must be' %filename
            msg += ' date-time with format %s.\n' %time_format
            msg += 'I got %s instead.' %fields[0]
            raise msg

        #Split values
        values = []
        for value in fields[1].split():
            values.append(numpy.float(value))

        q = array(values)

        msg = 'ERROR: File boundary function must return a 1d list or array '
        assert len(q.shape) == 1, msg

        d = len(domain.conserved_quantities)
        msg = 'Values specified in file must be a list or an array of length %d' %d
        assert len(q) == d, msg

        self.filename = filename
        self.domain = domain
        self.starttime = starttime

        if domain.starttime is None:
            domain.starttime = starttime
        else:
            msg = 'Start time as specified in domain (%s) is earlier '
            'than the starttime of file %s: %s'\
                  %(domain.starttime, self.filename, self.starttime)
            if self.starttime > domain.starttime:
                raise msg

        self.read_time_boundary() #Now read all times in.


    def read_time_boundary(self):
        from config import time_format
        import time

        fid = open(self.filename)
        lines = fid.readlines()
        fid.close()

        N = len(lines)
        d = len(self.domain.conserved_quantities)

        T = numpy.zeros(N, numpy.float)
        Q = numpy.zeros((N, d), numpy.float)


        for i, line in enumerate(lines):
            fields = line.split(',')
            real_time = time.mktime(time.strptime(fields[0], time_format))

            T[i] = real_time - self.starttime


            for j, value in enumerate(fields[1].split()):
                Q[i, j] = numpy.float(value)

        msg = 'Time boundary must list time as a monotonuosly '
        msg += 'increasing sequence'

        assert alltrue( T[1:] - T[:-1] > 0 ), msg

        self.T = T     #Time
        self.Q = Q     #Boundary values
        self.index = 0 #Initial index


    def __repr__(self):
        return 'File boundary'

    def evaluate(self, vol_id=None, edge_id=None):
        """Return linearly interpolated values based on domain.time
        """

        t = self.domain.time

        msg = 'Time given in File boundary does not match model time'
        if t < self.T[0]: raise msg
        if t > self.T[-1]: raise msg

        oldindex = self.index

        #Find slot
        while t > self.T[self.index]: self.index += 1
        while t < self.T[self.index]: self.index -= 1

        #if oldindex != self.index:
        #    print 'Changing from %d to %d' %(oldindex, self.index)


        #t is now between index and index+1
        ratio = (t - self.T[self.index])/\
                (self.T[self.index+1] - self.T[self.index])

        return self.Q[self.index,:] +\
               ratio*(self.Q[self.index+1,:] - self.Q[self.index,:])




