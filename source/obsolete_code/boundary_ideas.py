
#MOVED from generic boundary conditions 27 Nov 2008





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
    
        # FIXME (Ole): I think this should be get_time(), see ticket:306    
        return self.f(self.conserved_quantities,
                      neighbour.conserved_quantities,
                      self.domain.time)




