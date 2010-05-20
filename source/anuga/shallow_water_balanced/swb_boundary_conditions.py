
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Boundary

class Transmissive_boundary(Boundary):
    """Transmissive boundary returns same conserved quantities as
    those present in its neighbour volume.

    Underlying domain must be specified when boundary is instantiated
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for transmissive boundary'
            raise Exception, msg

        self.domain = domain

    def __repr__(self):
        return 'Transmissive_boundary(%s)' %self.domain

    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries return the edge values
        of the volume they serve.
        """

        if self.domain.get_centroid_transmissive_bc() :
            q = self.domain.get_evolved_quantities(vol_id)
        else:
            q = self.domain.get_evolved_quantities(vol_id, edge = edge_id)
        return q
