from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Boundary

import numpy as num

from anuga.shallow_water.shallow_water_ext import rotate


################################################################################
# Boundary conditions - specific to the shallow water wave equation
################################################################################

##
# @brief Class for a reflective boundary.
# @note Inherits from Boundary.
class Reflective_boundary(Boundary):
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    """

    ##
    # @brief Instantiate a Reflective_boundary.
    # @param domain 
    def __init__(self, model=None):
        Boundary.__init__(self)

        if model is None:
            msg = 'Model must be specified for reflective boundary'
            raise Exception, msg

        self.conserved_quantities = num.zeros(3, num.float)
        
        self.model = model

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Reflective_boundary'

    ##
    # @brief Calculate reflections (reverse outward momentum).
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Reflective boundaries reverses the outward momentum
        of the volume they serve.
        """

        # FIXME - arbitrary model construction time means we have to do these in the
        # RUN state.
        stage = self.model.domain.quantities['stage'].edge_values
        xmom = self.model.domain.quantities['xmomentum'].edge_values
        ymom = self.model.domain.quantities['ymomentum'].edge_values
        normals = self.model.domain.normals

        q = self.conserved_quantities
        q[0] = stage[vol_id, edge_id]
        q[1] = xmom[vol_id, edge_id]
        q[2] = ymom[vol_id, edge_id]

        normal = normals[vol_id, 2*edge_id:2*edge_id+2]

        r = rotate(q, normal, direction = 1)
        r[1] = -r[1]
        q = rotate(r, normal, direction = -1)

        return q


##
# @brief A transmissive boundary, momentum set to zero.
# @note Inherits from Bouondary.
class Transmissive_stage_zero_momentum_boundary(Boundary):
    """Return same stage as those present in its neighbour volume.
    Set momentum to zero.

    Underlying domain must be specified when boundary is instantiated
    """

    ##
    # @brief Instantiate a Transmissive (zero momentum) boundary.
    # @param domain
    def __init__(self, model=None):
        Boundary.__init__(self)

        if model is None:
            msg = ('Domain must be specified for '
                   'Transmissive_stage_zero_momentum boundary')
            raise Exception, msg

        self.model = model

    ##
    # @brief Return a representation of this instance.
    def __repr__(self):
        return 'Transmissive_stage_zero_momentum_boundary(%s)' % self.domain

    ##
    # @brief Calculate transmissive (zero momentum) results.
    # @param vol_id 
    # @param edge_id 
    def evaluate(self, vol_id, edge_id):
        """Transmissive boundaries return the edge values
        of the volume they serve.
        """

        q = self.model.domain.get_conserved_quantities(vol_id, edge=edge_id)

        q[1] = q[2] = 0.0
        return q
