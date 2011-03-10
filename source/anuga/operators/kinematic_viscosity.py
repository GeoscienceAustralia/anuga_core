
from anuga import Domain
from anuga import Quantity
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient
import anuga.abstract_2d_finite_volumes.neighbour_mesh as neighbour_mesh
from anuga import Dirichlet_boundary
import numpy as num
import kinematic_viscosity_ext
import anuga.utilities.log as log

class Kinematic_Viscosity:
    """
    Class for setting up structures and matrices for kinematic viscosity differential
    operator using centroid values.

    div ( diffusivity grad )

    where diffusvity is scalar quantity (defaults to quantity with values = 1)
    boundary values of f are used to setup entries associated with cells with boundaries

    There are procedures to apply this operator, ie

    (1) Calculate div( diffusivity grad u )
    using boundary values stored in u

    (2) Calculate ( u + dt div( diffusivity grad u )
    using boundary values stored in u

    (3) Solve div( diffusivity grad u ) = f
    for quantity f and using boundary values stored in u

    (4) Solve ( u + dt div( diffusivity grad u ) = f
    for quantity f using boundary values stored in u

    """

    def __init__(self, domain, diffusivity = None, triangle_areas=True, verbose=False):
        if verbose: log.critical('Kinematic Viscosity: Beginning Initialisation')
        #Expose the domain attributes

        self.domain = domain
        self.mesh = domain.mesh
        self.boundary = domain.boundary

        # Setup diffusivity quantity
        if diffusivity is None:
            diffusivity = Quantity(domain)
            diffusivity.set_values(1.0)
        #check that diffusivity is a quantity associated with domain
        assert diffusivity.domain == domain

        self.diffusivity = diffusivity
        self.diffusivity_bdry_data = self.diffusivity.boundary_values
        self.diffusivity_data = self.diffusivity.centroid_values


        self.n = len(self.domain)
        self.dt = 1e-6 #Need to set to domain.timestep
        self.boundary_len = len(domain.boundary)
        self.tot_len = self.n + self.boundary_len

        # FIXME SR: maybe this should already exist in mesh!
        self.boundary_enum = self.enumerate_boundary()
        
        self.verbose = verbose

        #Geometric Information
        if verbose: log.critical('Kinematic Viscosity: Building geometric structure')

        self.geo_structure_indices = num.zeros((self.n, 3), num.int)
        self.geo_structure_values = num.zeros((self.n, 3), num.float)

        # Only needs to built once, doesn't change
        kinematic_viscosity_ext.build_geo_structure(self, self.n, self.tot_len)

        # Setup type of scaling
        self.apply_triangle_areas = triangle_areas

        # FIXME SR: should this really be a matrix?
        temp  = Sparse(self.n, self.n)
        for i in range(self.n):
            temp[i, i] = 1.0 / self.mesh.areas[i]
            
        self.triangle_areas = Sparse_CSR(temp)

        # FIXME SR: More to do with solving equation
        self.qty_considered = 1 #1 or 2 (uh or vh respectively)

        #FIXME SR: Do we need Sparse or should we just go with Sparse_CSR
        self.operator_matrix = Sparse(self.n, self.tot_len)

        #Sparse_CSR.data
        self.operator_data = num.zeros((4 * self.n, ), num.float)
        #Sparse_CSR.colind
        self.operator_colind = num.zeros((4 * self.n, ), num.int)
        #Sparse_CSR.rowptr (4 entries in every row, we know this already) = [0,4,8,...,4*n]
        self.operator_rowptr = 4 * num.arange(self.n + 1)

        # Build matrix self.elliptic_operator_matrix [A B]
        self.build_elliptic_operator()

        # Build self.boundary_term
        self.build_operator_boundary_term()

        self.parabolic_solve = False #Are we doing a parabolic solve at the moment?

        if verbose: log.critical('Kinematic Viscosity: Initialisation Done')

        
    def enumerate_boundary(self):
        #Enumerate by boundary index
        # FIXME SR: Probably should be part of mesh
        enumeration = {}
        for i in range(self.n):
            for edge in range(3):
                j = self.mesh.neighbours[i, edge]
                if j < 0:
                    enumeration[(i, edge)] = -j-1
        return enumeration

    def set_qty_considered(self, qty):
        # FIXME SR: Probably should just be set by quantity to which operation is applied
        if qty == 1 or qty == 'u':
            self.qty_considered = 1
        elif qty == 2 or qty == 'v':
            self.qty_considered = 2
        else: #Raise an exception
            msg = "Incorrect input qty"
            assert 0 == 1, msg

    def build_elliptic_operator(self):
        """
        Builds matrix representing

        div ( diffusivity grad )

        which has the form [ A B ]
        """

        #Arrays self.operator_data, self.operator_colind, self.operator_rowptr
        # are setup via this call
        kinematic_viscosity_ext.build_operator_matrix(self,self.diffusivity_data, self.diffusivity_bdry_data)

        self.elliptic_operator_matrix = Sparse_CSR(None, \
                self.operator_data, self.operator_colind, self.operator_rowptr, \
                self.n, self.tot_len)

#        #Set up the scaling matrix
#        data = h
#        num.putmask(data, data != 0, 1 / data) #take the reciprocal of each entry unless it is zero
#        self.stage_heights_scaling = \
#            Sparse_CSR(None, self.diffusivity_data, num.arange(self.n), num.arange(self.n +1), self.n, self.n)

    def update_elliptic_operator(self):
        """
        Updates the data values of matrix representing

        div ( diffusivity grad )

        (as when diffusivitiy has changed)
        """

        #Array self.operator_data is changed by this call, which should flow
        # through to the Sparse_CSR matrix.

        kinematic_viscosity_ext.update_operator_matrix(self, self.n, self.tot_len, \
          self.diffusivity_data, self.diffusivity_bdry_data)




    def build_operator_boundary_term(self):
        """
        Operator has form [A B] and U = [ u ; b]

        This procedure calculates B b which can be calculated as

        [A B] [ 0 ; b]
        """

        n = self.n
        tot_len = self.tot_len

        self.boundary_term = num.zeros((n, ), num.float)
        
        X = num.zeros((tot_len,), num.float)

        X[n:] = self.diffusivity_bdry_data
        self.boundary_term[:] = self.operator_matrix * X

        #Tidy up
        if self.apply_triangle_areas:
            self.boundary_term[:] = self.triangle_areas * self.boundary_term


    def elliptic_multiply(self, V, qty_considered=None, include_boundary=True):
        msg = "(KV_Operator.apply_vector) V vector has incorrect dimensions"
        assert V.shape == (self.n, ) or V.shape == (self.n, 1), msg

        if qty_considered != None:
            print "Quantity Considered changed!"
            self.set_qty_considered(qty_considered)

        D = num.zeros((self.n, ), num.float)
        #Change (uh) to (u), (vh) to (v)
        X = self.stage_heights_scaling * V

        if self.apply_triangle_areas:
            X = self.triangle_areas * X

        X = num.append(X, num.zeros((self.boundary_len, )), axis=0)
        D = self.operator_matrix * X

        if include_boundary:
            D += self.boundary_vector[:, self.qty_considered-1]
        #make sure we return (uh) not (u), for the solver's benefit
        D = self.stage_heights * D

        return num.array(D).reshape(self.n, )

    #parabolic_multiply(V) = identity - dt*elliptic_multiply(V)
    #We do not need include boundary, as we will only use this
    #method for solving (include_boundary = False)
    #Here V is already either (u) or (v), not (uh) or (vh)
    def parabolic_multiply(self, V):
        msg = "(KV_Operator.parabolic_multiply) V vector has incorrect dimensions"
        assert V.shape == (self.n, ) or V.shape == (self.n, 1), msg

        D = V #The return value
        X = V #a temporary array

        #Apply triangle areas
        if self.apply_triangle_areas:
            X = self.triangle_areas * X

        #Multiply out
        X = num.append(X, num.zeros((self.boundary_len, )), axis=0)
        D = D - self.dt * (self.operator_matrix * X)

        return num.array(D).reshape(self.n, )

    def __mul__(self, other):
        try:
            B = num.array(other)
        except:
            msg = "Trying to multiply the Kinematic Viscosity Operator against a non-numeric type"
            raise msg

        if len(B.shape) == 0:
            #Scalar
            R = B * self
        elif len(B.shape) == 1 or (len(B.shape) == 2 and B.shape[1] == 1):
            #Vector
            if self.parabolic_solve:
                R = self.parabolic_multiply(other)
            else:
                #include_boundary=False is this is *only* used for cg_solve()
                R = self.elliptic_multiply(other, include_boundary=False)
        else:
            raise ValueError, 'Dimension too high: d=%d' % len(B.shape)
        return R
    
    def __rmul__(self, other):
        #Right multiply with scalar
        try:
            other = float(other)
        except:
            msg = 'Sparse matrix can only "right-multiply" onto a scalar'
            raise TypeError, msg
        else:
            new = self.operator_matrix * new
        return new

    def cg_solve(self, B, qty_to_consider=None):
        if len(B.shape) == 1:
            return self.cg_solve_vector(B, qty_to_consider)
        elif len(B.shape) == 2:
            return self.cg_solve_matrix(B)
        else:
            raise ValueError, 'Dimension too high: d=%d' % len(B.shape)

    def cg_solve_matrix(self, B):
        assert B.shape[1] < 3, "Input matrix has too many columns (max 2)"

        X = num.zeros(B.shape, num.float)
        for i in range(B.shape[1]):
            X[:, i] = self.cg_solve_vector(B[:, i], i + 1) #assuming B columns are (uh) and (vh)
        return X
    
    def cg_solve_vector(self, b, qty_to_consider=None):
        if not qty_to_consider == None:
            self.set_qty_considered(qty_to_consider)
            #Call the ANUGA conjugate gradient utility
        x = conjugate_gradient(self, b - self.boundary_vector[:, self.qty_considered-1])
        return x

    #Solve the parabolic equation to find u^(n+1), where u = u^n
    #Here B = [uh vh] where uh,vh = n*1 column vectors
    def parabolic_solver(self, B):
        assert B.shape == (self.n, 2), "Input matrix has incorrect dimensions"

        next_B = num.zeros((self.n, 2), num.float)
        #The equation is self.parabolic_multiply*u^(n+1) = u^n + (dt)*self.boundary_vector
        self.parabolic_solve = True
        #Change (uh) and (vh) to (u) and (v)
        b = self.stage_heights_scaling * B

        next_B = conjugate_gradient(self, b + (self.dt * self.boundary_vector), iprint=1)

        #Return (u) and (v) back to (uh) and (vh)
        next_B = self.stage_heights_scaling * next_B

        self.parabolic_solve = False

        return next_B
