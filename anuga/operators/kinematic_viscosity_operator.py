
from anuga import Quantity
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient
import anuga.abstract_2d_finite_volumes.neighbour_mesh as neighbour_mesh
from anuga import Dirichlet_boundary
import numpy as num
from . import kinematic_viscosity_operator_ext
import anuga.utilities.log as log

from anuga.operators.base_operator import Operator



class Kinematic_viscosity_operator(Operator):
    """
    Class for setting up structures and matrices for kinematic viscosity differential
    operator using centroid values.


    As an anuga operator, when the __call__ method is called one step of the parabolic
    step is applied. In particular the x and y velocities are updated using

    du/dt = div( h grad u )
    dv/dt = div( h grad v )

    """

    def __init__(self,
                 domain, diffusivity='height',
                 use_triangle_areas=True,
                 add_safety = False,
                 verbose=False):

        if verbose: log.critical('Kinematic Viscosity: Beginning Initialisation')
        

        Operator.__init__(self,domain)

        #Expose the domain attributes
        self.mesh = self.domain.mesh
        self.boundary = domain.boundary
        self.boundary_enumeration = domain.boundary_enumeration

        self.diffusivity = diffusivity
        # Setup a quantity as diffusivity
        if diffusivity is None:
            self.diffusivity = Quantity(self.domain)
            self.diffusivity.set_values(1.0)
            self.diffusivity.set_boundary_values(1.0)

        if isinstance(diffusivity, (int, float)):
            self.diffusivity = Quantity(self.domain)
            self.diffusivity.set_values(diffusivity)
            self.diffusivity.set_boundary_values(diffusivity)

        if isinstance(diffusivity, str):
            self.diffusivity = self.domain.get_quantity(diffusivity)
            

        self.add_safety = add_safety
        self.smooth = 0.1
        
        assert isinstance(self.diffusivity, Quantity)
        

        self.n = len(self.domain)

        self.dt = 0.0 #Need to set to domain.timestep
        self.dt_apply = 0.0

        self.boundary_len = len(self.domain.boundary)
        self.tot_len = self.n + self.boundary_len

        self.verbose = verbose

        #Geometric Information
        if verbose: log.critical('Kinematic Viscosity: Building geometric structure')

        self.geo_structure_indices = num.zeros((self.n, 3), int)
        self.geo_structure_values = num.zeros((self.n, 3), float)

        # Only needs to built once, doesn't change
        kinematic_viscosity_operator_ext.build_geo_structure(self)

        # Setup type of scaling
        self.set_triangle_areas(use_triangle_areas)        

        # FIXME SR: should this really be a matrix?
        temp  = Sparse(self.n, self.n)
        for i in range(self.n):
            temp[i, i] = 1.0 / self.mesh.areas[i]
            
        self.triangle_areas = Sparse_CSR(temp)
        #self.triangle_areas

        # FIXME SR: More to do with solving equation
        self.qty_considered = 1 #1 or 2 (uh or vh respectively)

        #Sparse_CSR.data
        self.operator_data = num.zeros((4 * self.n, ), float)
        #Sparse_CSR.colind
        self.operator_colind = num.zeros((4 * self.n, ), int)
        #Sparse_CSR.rowptr (4 entries in every row, we know this already) = [0,4,8,...,4*n]
        self.operator_rowptr = 4 * num.arange(self.n + 1)

        # Build matrix self.elliptic_matrix [A B]
        self.build_elliptic_matrix(self.diffusivity)

        self.boundary_term = num.zeros((self.n, ), float)

        self.parabolic = False #Are we doing a parabolic solve at the moment?

        self.u_stats = None
        self.v_stats = None

        if verbose: log.critical('Kinematic Viscosity: Initialisation Done')


    def __call__(self):
        """ Parabolic update of x and y velocity

        (I + dt (div d grad) ) U^{n+1} = U^{n}

        """

        domain = self.domain

        self.dt = self.dt + domain.get_timestep()


        if self.dt < self.dt_apply:
            return


        # Setup initial values of velocity
        domain.update_centroids_of_velocities_and_height()

        # diffusivity
        if self.add_safety:
            d = self.diffusivity + 0.1
        else:
            d = self.diffusivity


        assert num.all(d.centroid_values >= 0.0)

        #d.set_values(num.where(h.centroid_values<1.0, 0.0, 1.0), location= 'centroids')
        #d.set_boundary_values(num.where(h.boundary_values<1.0, 0.0, 1.0))

        # Quantities to solve
        # Boundary values are implied from BC set in advection part of code
        u = domain.quantities['xvelocity']
        #u.set_boundary_values(0.0)

        v = domain.quantities['yvelocity']
        #v.set_boundary_values(0.0)

        #Update operator using current height
        self.update_elliptic_matrix(d)

        (u, self.u_stats) = self.parabolic_solve(u, u, d, u_out=u, update_matrix=False, output_stats=True)

        (v, self.v_stats) = self.parabolic_solve(v, v, d, u_out=v, update_matrix=False, output_stats=True)

        # Update the conserved quantities
        domain.update_centroids_of_momentum_from_velocity()


        self.dt = 0.0



    def statistics(self):

        message = 'Kinematic_viscosity_operator '
        return message


    def timestepping_statistics(self):

        from anuga import indent

        message = indent+'Kinematic Viscosity Operator: \n'
        if self.u_stats is not None:
            message  += indent + indent + 'u: ' + self.u_stats.__str__() +'\n'

        if self.v_stats is not None:
            message  += indent + indent + 'v: ' + self.v_stats.__str__()

        return message


    def set_triangle_areas(self,flag=True):

        self.apply_triangle_areas = flag
        

    def set_parabolic_solve(self,flag):

        self.parabolic = flag


    def build_elliptic_matrix(self, a):
        """
        Builds matrix representing

        div ( a grad )

        which has the form [ A B ]
        """

        #Arrays self.operator_data, self.operator_colind, self.operator_rowptr
        # are setup via this call
        kinematic_viscosity_operator_ext.build_elliptic_matrix(self, \
                a.centroid_values, \
                a.boundary_values)

        self.elliptic_matrix = Sparse_CSR(None, \
                self.operator_data, self.operator_colind, self.operator_rowptr, \
                self.n, self.tot_len)


    def update_elliptic_matrix(self, a=None):
        """
        Updates the data values of matrix representing

        div ( a grad )

        If a is None then we set a = quantity which is set to 1
        """

        #Array self.operator_data is changed by this call, which should flow
        # through to the Sparse_CSR matrix.

        if a is None:
            a = Quantity(self.domain)
            a.set_values(1.0)
            a.set_boundary_values(1.0)
            
        kinematic_viscosity_operator_ext.update_elliptic_matrix(self, \
                a.centroid_values, \
                a.boundary_values)
        




    def update_elliptic_boundary_term(self, boundary):


        if isinstance(boundary, Quantity):

            self._update_elliptic_boundary_term(boundary.boundary_values)

        elif isinstance(boundary, num.ndarray):

            self._update_elliptic_boundary_term(boundary)

        else:

            raise  TypeError('expecting quantity or numpy array')


    def _update_elliptic_boundary_term(self, b):
        """
        Operator has form [A B] and u = [ u_1 ; b]

        u_1 associated with centroid values of u
        u_2 associated with boundary_values of u

        This procedure calculates B u_2 which can be calculated as

        [A B] [ 0 ; b]

        Assumes that update_elliptic_matrix has just been run.
        """

        n = self.n
        tot_len = self.tot_len

        X = num.zeros((tot_len,), float)

        X[n:] = b
        self.boundary_term[:] = self.elliptic_matrix * X

        #Tidy up
        if self.apply_triangle_areas:
            self.boundary_term[:] = self.triangle_areas * self.boundary_term



    def elliptic_multiply(self, input, output=None):


        if isinstance(input, Quantity):

            assert isinstance(output, Quantity) or output is None

            output = self._elliptic_multiply_quantity(input, output)

        elif isinstance(input, num.ndarray):

            assert isinstance(output, num.ndarray) or output is None

            output = self._elliptic_multiply_array(input, output)

        else:

            raise TypeError('expecting quantity or numpy array')
        
        return output


    def _elliptic_multiply_quantity(self, quantity_in, quantity_out):
        """
        Assumes that update_elliptic_matrix has been run
        """

        if quantity_out is None:
            quantity_out = Quantity(self.domain)

        array_in = quantity_in.centroid_values
        array_out = quantity_out.centroid_values

        X = self._elliptic_multiply_array(array_in, array_out)

        quantity_out.set_values(X, location = 'centroids')
        
        return quantity_out

    def _elliptic_multiply_array(self, array_in, array_out):
        """
        calculates [A B] [array_in ; 0]
        """

        n = self.n
        tot_len = self.tot_len

        V = num.zeros((tot_len,), float)

        assert len(array_in) == n

        if array_out is None:
            array_out = num.zeros_like(array_in)

        V[0:n] = array_in
        V[n:] = 0.0


        if self.apply_triangle_areas:
            V[0:n] = self.triangle_areas * V[0:n]


        array_out[:] = self.elliptic_matrix * V


        return array_out





    def parabolic_multiply(self, input, output=None):


        if isinstance(input, Quantity):

            assert isinstance(output, Quantity) or output is None

            output = self._parabolic_multiply_quantity(input, output)

        elif isinstance(input, num.ndarray):

            assert isinstance(output, num.ndarray) or output is None

            output = self._parabolic_multiply_array(input, output)

        else:

            raise TypeError('expecting quantity or numpy array')

        return output


    def _parabolic_multiply_quantity(self, quantity_in, quantity_out):
        """
        Assumes that update_elliptic_matrix has been run
        """

        if quantity_out is None:
            quantity_out = Quantity(self.domain)

        array_in = quantity_in.centroid_values
        array_out = quantity_out.centroid_values

        X = self._parabolic_multiply_array(array_in, array_out)

        quantity_out.set_values(X, location = 'centroids')

        return quantity_out

    def _parabolic_multiply_array(self, array_in, array_out):
        """
        calculates ( [ I 0 ; 0  0] + dt [A B] ) [array_in ; 0]
        """

        n = self.n
        tot_len = self.tot_len

        V = num.zeros((tot_len,), float)

        assert len(array_in) == n

        if array_out is None:
            array_out = num.zeros_like(array_in)

        V[0:n] = array_in
        V[n:] = 0.0


        if self.apply_triangle_areas:
            V[0:n] = self.triangle_areas * V[0:n]


        array_out[:] = array_in - self.dt * (self.elliptic_matrix * V)

        return array_out




    def __mul__(self, vector):
        
        #Vector
        if self.parabolic:
            R = self.parabolic_multiply(vector)
        else:
            #include_boundary=False is this is *only* used for cg_solve()
            R = self.elliptic_multiply(vector)

        return R
    
    def __rmul__(self, other):
        #Right multiply with scalar
        try:
            other = float(other)
        except:
            msg = 'Sparse matrix can only "right-multiply" onto a scalar'
            raise TypeError(msg)
        else:
            new = self.elliptic_matrix * new
        return new

    
    def elliptic_solve(self, u_in, b, a = None, u_out = None, update_matrix=True, \
                       imax=10000, tol=1.0e-8, atol=1.0e-8,
                       iprint=None, output_stats=False):
        """ Solving div ( a grad u ) = b
        u | boundary = g

        u_in, u_out, f anf g are Quantity objects

        Dirichlet BC g encoded into u_in boundary_values

        Initial guess for iterative scheme is given by
        centroid values of u_in

        Centroid values of a and b provide diffusivity and rhs

        Solution u is retruned in u_out
        """

        if u_out is None:
            u_out = Quantity(self.domain)

        if update_matrix :
            self.update_elliptic_matrix(a) 

        self.update_elliptic_boundary_term(u_in)

        # Pull out arrays and a matrix operator
        A = self
        rhs = b.centroid_values - self.boundary_term
        x0 = u_in.centroid_values

        x, stats = conjugate_gradient(A,rhs,x0,imax=imax, tol=tol, atol=atol,
                               iprint=iprint, output_stats=True)

        u_out.set_values(x, location='centroids')
        u_out.set_boundary_values(u_in.boundary_values)

        if output_stats:
            return u_out, stats
        else:
            return u_out

    

    def parabolic_solve(self, u_in, b, a = None, u_out = None, update_matrix=True, \
                       output_stats=False, use_dt_tol=True, iprint=None, imax=10000):
        """
        Solve for u in the equation

        ( I + dt div a grad ) u = b

        u | boundary = g


        u_in, u_out, f anf g are Quantity objects

        Dirichlet BC g encoded into u_in boundary_values

        Initial guess for iterative scheme is given by
        centroid values of u_in

        Centroid values of a and b provide diffusivity and rhs

        Solution u is retruned in u_out

        """


        if use_dt_tol:
            tol  = min(self.dt,1.0e-5)
            atol = min(self.dt,1.0e-5)
        else:
            tol  = 1.0e-5
            atol = 1.0e-5



        
        if u_out is None:
            u_out = Quantity(self.domain)

        if update_matrix :
            self.update_elliptic_matrix(a)

        self.update_elliptic_boundary_term(u_in)

        self.set_parabolic_solve(True)


        # Pull out arrays and a matrix operator
        IdtA = self
        rhs = b.centroid_values + (self.dt * self.boundary_term)
        x0 = u_in.centroid_values

        x, stats = conjugate_gradient(IdtA,rhs,x0,imax=imax, tol=tol, atol=atol,
                                      iprint=iprint, output_stats=True)

        self.set_parabolic_solve(False)

        u_out.set_values(x, location='centroids')
        u_out.set_boundary_values(u_in.boundary_values)

        if output_stats:
            return u_out, stats
        else:
            return u_out

 
