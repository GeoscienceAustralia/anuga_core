#! /usr/bin/python

"""
Boundaries - specific to the shallow water wave equation
"""
__author__="Stephen Roberts"
__date__ ="$05/06/2010 5:44:05 PM$"

from anuga_1d.base.generic_domain import *

class Dirichlet_boundary(Boundary):
    """Dirichlet boundary returns constant values for the
    conserved quantities
    """

    def __init__(self, quantities=None):
        Boundary.__init__(self)

        
        if quantities is None:
            msg = 'Must specify one value for each evolved quantities' + \
                  ' area, discharge, elevation, height, velocity, width, top, stage'
            raise Exception(msg)
        
        msg = 'Must specify one value for each evolved quantities' + \
             ' area, discharge, elevation, height, velocity, width, top, stage'
        assert len(quantities)==8, msg

        self.quantities=numpy.array(quantities, numpy.float)

    def __repr__(self):
        return 'Dirichlet boundary (%s)' %self.quantities

    def evaluate(self, vol_id=None, edge_id=None):
        return self.quantities


class Reflective_boundary(Boundary):
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    conserved quantities.
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for Reflective boundary'
            raise Exception(msg)

        #Handy shorthands
        self.normals    = domain.normals
        self.area       = domain.quantities['area'].vertex_values
        self.discharge  = domain.quantities['discharge'].vertex_values
        self.bed        = domain.quantities['elevation'].vertex_values
        self.height     = domain.quantities['height'].vertex_values
        self.velocity   = domain.quantities['velocity'].vertex_values
        self.width      = domain.quantities['width'].vertex_values
        self.top        = domain.quantities['top'].vertex_values
        self.stage      = domain.quantities['stage'].vertex_values

        self.quantities = numpy.zeros(8, numpy.float)

    def __repr__(self):
        return 'Reflective_boundary'


    def evaluate(self, vol_id, edge_id):
        """Reflective boundaries reverses the outward momentum
        of the volume they serve.
        """

        q = self.quantities
        q[0] =  self.area[vol_id, edge_id]
        q[1] = -self.discharge[vol_id, edge_id]
        q[2] =  self.bed[vol_id, edge_id]
        q[3] =  self.height[vol_id, edge_id]
        q[4] = -self.velocity[vol_id, edge_id]
        q[5] =  self.width[vol_id, edge_id]
        q[6] =  self.top[vol_id, edge_id]
        q[7] =  self.stage[vol_id, edge_id]

        #normal = self.normals[vol_id,edge_id]

        return q


class Transmissive_boundary(Boundary):
    """Transmissive boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    quantities.
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for Transmissive boundary'
            raise Exception(msg)

        #Handy shorthands
        self.normals    = domain.normals
        self.area       = domain.quantities['area'].vertex_values
        self.discharge  = domain.quantities['discharge'].vertex_values
        self.bed        = domain.quantities['elevation'].vertex_values
        self.height     = domain.quantities['height'].vertex_values
        self.velocity   = domain.quantities['velocity'].vertex_values
        self.width      = domain.quantities['width'].vertex_values
        self.top        = domain.quantities['top'].vertex_values
        self.stage      = domain.quantities['stage'].vertex_values
        
        self.quantities = numpy.zeros(8, numpy.float)

    def __repr__(self):
        return 'Transmissive_boundary'


    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries mimics the outward momentum
        of the volume they serve.
        """

        q = self.quantities
        q[0] =  self.area[vol_id, edge_id]
        q[1] =  self.discharge[vol_id, edge_id]
        q[2] =  self.bed[vol_id, edge_id]
        q[3] =  self.height[vol_id, edge_id]
        q[4] =  self.velocity[vol_id, edge_id]
        q[5] =  self.width[vol_id, edge_id]
        q[6] =  self.top[vol_id, edge_id]
        q[7] =  self.stage[vol_id, edge_id]

        #normal = self.normals[vol_id,edge_id]

        return q



if __name__ == "__main__":
    print "Hello World";
