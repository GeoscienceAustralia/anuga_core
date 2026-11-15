
import warnings

from anuga import Quantity
from anuga.utilities.sparse import Sparse, Sparse_CSR
from anuga.utilities.cg_solve import conjugate_gradient, Stats, ConvergenceError
from anuga.utilities.cg_ext import cg_solve_c_precon, jacobi_precon_c
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

    du/dt = div( d grad u )
    dv/dt = div( d grad v )

    where d is the diffusivity: the water depth h by default
    (diffusivity='height'), or a constant or quantity given by the
    diffusivity argument. The velocities, not the momenta, are diffused, so
    where the depth varies this form does not conserve momentum.

    With form='momentum' the operator solves the momentum-conserving form

    d(hu)/dt = div( nu h grad u )
    d(hv)/dt = div( nu h grad v )

    where nu is the eddy viscosity given by the diffusivity argument (a
    number, a Quantity or a quantity name, not 'height'). An edge between a
    wet and a dry cell carries no stress, so total momentum changes only
    through the domain boundary. In either form a Reflective boundary is a
    free-slip wall and carries no stress either; other boundaries prescribe
    the velocity at the edge.

    Either form is solved implicitly, one step per call:

    A_i m_i u_i - dt sum_edges a_e (u_j - u_i) l_e / d_e = A_i m_i u_i^old

    with A_i the triangle area, m_i = 1 (velocity form) or h_i (momentum
    form), a_e the mean of the diffusivity (or of nu h) across the edge, l_e
    its length and d_e the distance between the centroids. The matrix is
    symmetric and positive definite. Each momentum then changes by
    h (u - u^old).

    With scheme='explicit' (form='momentum' and a constant nu only) the
    momentum form is instead advanced by forward Euler, with the edge depth
    taken as the harmonic mean 2 h_i h_j / (h_i + h_j). There is no linear
    solve: one pass over the edges, which runs on the GPU in unified mode.
    The step is stable while 2 nu dt sum_e (l_e / d_e) / A_i <= 1 for every
    cell; the operator splits a longer timestep into equal sub-steps that
    meet it. At the timesteps the flow solver takes nu dt / dx^2 is usually
    far below this, so one step is enough.

    For a physical eddy viscosity use form='momentum', scheme='explicit'
    with nu in m^2/s: it is the form the depth-averaged equations give,
    conserves momentum, and runs on the GPU. The velocity form, the
    default, is a smoother: it diffuses velocity rather than momentum,
    which creates momentum where the depth varies. Not passing form gives
    a FutureWarning, as the default may change; pass form='velocity' to
    keep it.

    """

    def __init__(self,
                 domain, diffusivity='height',
                 use_triangle_areas=True,
                 add_safety = False,
                 form=None,
                 scheme='implicit',
                 verbose=False):

        if verbose: log.info('Kinematic Viscosity: Beginning Initialisation')


        Operator.__init__(self,domain)

        if form is None:
            warnings.warn(
                "Kinematic_viscosity_operator: no form given, so using the "
                "velocity form, du/dt = div(d grad u), which does not conserve "
                "momentum where the depth varies. Pass form='velocity' to keep "
                "it, or form='momentum' with an eddy viscosity nu (and "
                "scheme='explicit') for the momentum-conserving form. The "
                "default may change in a future release.",
                FutureWarning, stacklevel=2)
            form = 'velocity'
        if form not in ('velocity', 'momentum'):
            raise ValueError("form must be 'velocity' or 'momentum', not %r" % (form,))
        if form == 'momentum' and isinstance(diffusivity, str) and diffusivity == 'height':
            raise ValueError("form='momentum' needs the eddy viscosity nu as the "
                             "diffusivity (a number, Quantity or quantity name); "
                             "the depth h is already part of the momentum form")
        if scheme not in ('implicit', 'explicit'):
            raise ValueError("scheme must be 'implicit' or 'explicit', not %r" % (scheme,))
        if scheme == 'explicit':
            if form != 'momentum':
                raise ValueError("scheme='explicit' needs form='momentum'")
            if isinstance(diffusivity, bool) or not isinstance(diffusivity, (int, float)):
                raise ValueError("scheme='explicit' needs a constant eddy viscosity "
                                 "nu (a number) as the diffusivity")
            if diffusivity < 0:
                raise ValueError('the eddy viscosity must not be negative')
        self.form = form
        self.scheme = scheme

        if scheme == 'explicit':
            self._init_explicit(float(diffusivity), verbose)
            return

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

        # n: total local triangles (full + ghost in parallel; all in serial)
        # n_full: full (owned) triangles only
        self.n = len(self.domain)
        self.n_full = domain.number_of_full_triangles

        self.dt = 0.0 #Need to set to domain.timestep
        self.dt_apply = 0.0

        self.boundary_len = len(self.domain.boundary)
        self.tot_len = self.n + self.boundary_len

        self.verbose = verbose

        #Geometric Information
        if verbose: log.info('Kinematic Viscosity: Building geometric structure')

        self.geo_structure_indices = num.zeros((self.n, 3), int)
        self.geo_structure_values = num.zeros((self.n, 3), float)

        # Only needs to built once, doesn't change
        kinematic_viscosity_operator_ext.build_geo_structure(self)
        self._geo_structure_base = self.geo_structure_values.copy()
        self._geo_boundary_objects = None
        self._geo_walls_set = False

        # Setup type of scaling
        self.set_triangle_areas(use_triangle_areas)

        # FIXME SR: should this really be a matrix?
        temp  = Sparse(self.n, self.n)
        for i in range(self.n):
            temp[i, i] = 1.0 / self.mesh.areas[i]

        self.triangle_areas = Sparse_CSR(temp)

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

        # Pre-allocate workspace for the distributed CG matvec
        self._kv_V = num.zeros(self.tot_len, dtype=float)

        if verbose: log.info('Kinematic Viscosity: Initialisation Done')


    def _init_explicit(self, nu, verbose):
        """Set up the explicit scheme: no matrix, just the stability bound."""

        domain = self.domain
        self.mesh = domain.mesh
        self.nu = nu
        self.verbose = verbose
        self.n = len(domain)
        self.u_stats = None
        self.v_stats = None
        self.substeps = 0          # sub-steps taken by the last call
        self._scratch = {}         # 'host' / 'gpu' -> scratch handle
        self._scratch_walls = {}   # 'host' / 'gpu' -> boundary_objects last set

        # G = max over cells of sum_e (l_e / d_e) / A_i, the geometric part
        # of the stability bound 2 nu dt G <= 1. The same on every rank, so
        # that all ranks take the same number of sub-steps.
        mesh = self.mesh
        nbr = mesh.neighbours
        cc = mesh.centroid_coordinates
        em = mesh.edge_midpoint_coordinates.reshape(-1, 3, 2)
        other = num.where((nbr >= 0)[:, :, None], cc[num.maximum(nbr, 0)], em)
        dist = num.sqrt(((other - cc[:, None, :])**2).sum(axis=2))
        g = (mesh.edgelengths / dist).sum(axis=1) / mesh.areas
        G = num.array([g.max() if len(g) else 0.0])
        if getattr(domain, 'numproc', 1) > 1:
            import anuga.utilities.parallel_abstraction as pypar
            from mpi4py.MPI import MAX
            G_all = num.zeros(1)
            pypar.comm.Allreduce(G, G_all, op=MAX)
            G = G_all
        self.stability_G = float(G[0])

        if verbose:
            log.info('Kinematic Viscosity: explicit, nu = %g, stable for dt <= %g'
                     % (nu, self.max_stable_timestep()))

    def max_stable_timestep(self):
        """Longest timestep the explicit scheme takes in a single step."""

        if self.nu <= 0.0 or self.stability_G <= 0.0:
            return float('inf')
        return 1.0 / (2.0 * self.nu * self.stability_G)

    def _explicit_step(self):
        """One explicit step of the momentum form, on the device in mode 2."""

        domain = self.domain
        dt = domain.get_timestep()
        if dt <= 0.0 or self.nu <= 0.0:
            self.substeps = 0
            return

        # As for the sediment operator: inside the host-coherent window that
        # a CPU-only operator opens (_gpu_host_writes_suppressed) the host
        # copy is the one that gets pushed back, so work on the host there.
        on_gpu = (domain.multiprocessor_mode == 2
                  and getattr(domain, 'gpu_interface', None) is not None
                  and not getattr(domain, '_gpu_host_writes_suppressed', False))
        key = 'gpu' if on_gpu else 'host'
        if on_gpu:
            from anuga.shallow_water.sw_domain_gpu_ext import (
                kv_scratch_alloc_gpu as alloc, kv_scratch_set_walls_gpu as set_walls,
                kinematic_viscosity_explicit_gpu as step)
            target = domain.gpu_interface.gpu_dom
        else:
            from anuga.shallow_water.sw_domain_openmp_ext import (
                kv_scratch_alloc as alloc, kv_scratch_set_walls as set_walls,
                kinematic_viscosity_explicit as step)
            target = domain
        if key not in self._scratch:
            self._scratch[key] = alloc(self.n, len(domain.boundary))
        scratch = self._scratch[key]
        objects = getattr(domain, 'boundary_objects', None)
        if self._scratch_walls.get(key) is not objects or key not in self._scratch_walls:
            flags = num.zeros(len(domain.boundary), dtype=float)
            flags[[k for k, _vol, _edge in self._wall_edges()]] = 1.0
            set_walls(scratch, self.n, flags)
            self._scratch_walls[key] = objects

        nsub = max(1, int(num.ceil(dt / self.max_stable_timestep())))
        self.substeps = nsub
        for _ in range(nsub):
            # Fractional steps run before the step's ghost exchange, and the
            # stress on an edge shared with another rank needs the ghost's
            # current momentum.
            if domain.parallel:
                if on_gpu:
                    domain.gpu_interface.exchange_ghosts()
                else:
                    domain.update_ghosts()
            step(target, self.nu, dt / nsub, scratch)

    def _wall_edges(self):
        """(boundary index, triangle, edge) of each Reflective boundary edge.

        A Reflective wall is free-slip, so it exerts no viscous stress: it
        keeps the flow off the wall through the pressure and the fluxes alone.
        Its boundary values mirror the flow rather than prescribe a velocity,
        and using them as one, a sixth of a cell from the centroid on a
        rectangular cross mesh, swamps the stress between the cells.
        """
        from anuga.shallow_water.boundaries import Reflective_boundary

        objects = getattr(self.domain, 'boundary_objects', None) or []
        return [(k, vol, edge) for k, ((vol, edge), B) in enumerate(objects)
                if isinstance(B, Reflective_boundary)]

    def _update_wall_edges(self):
        """Take Reflective edges out of the geometry when the boundaries change."""
        objects = getattr(self.domain, 'boundary_objects', None)
        if objects is self._geo_boundary_objects and self._geo_walls_set:
            return
        self.geo_structure_values[:] = self._geo_structure_base
        for _k, vol, edge in self._wall_edges():
            self.geo_structure_values[vol, edge] = 0.0
        self._geo_boundary_objects = objects
        self._geo_walls_set = True

    def __del__(self):
        # Host scratch only; device scratch is released with the device
        # context, as the GPU runtime may already be gone at interpreter exit.
        handle = getattr(self, '_scratch', {}).get('host')
        if handle is not None:
            try:
                from anuga.shallow_water.sw_domain_openmp_ext import kv_scratch_free
                kv_scratch_free(handle)
            except Exception:
                pass

    def parallel_safe(self):
        return True


    def __call__(self):
        """ Parabolic update of x and y velocity

        (I + dt (div d grad) ) U^{n+1} = U^{n}

        """

        if self.scheme == 'explicit':
            self._explicit_step()
            return

        domain = self.domain

        self.dt = self.dt + domain.get_timestep()


        if self.dt < self.dt_apply:
            return


        self._update_wall_edges()

        # The solve runs on the host and reads boundary values, which mode 2
        # keeps on the device. Without them every open edge looks dry and still.
        if domain.multiprocessor_mode == 2 and getattr(domain, 'gpu_interface', None) is not None:
            domain.gpu_interface.sync_boundary_values_from_device()

        # Setup initial values of velocity
        domain.update_centroids_of_velocities_and_height()

        # diffusivity
        if self.add_safety:
            d = self.diffusivity + 0.1
        else:
            d = self.diffusivity


        assert num.all(d.centroid_values >= 0.0)

        # Quantities to solve
        # Boundary values are implied from BC set in advection part of code
        u = domain.quantities['xvelocity']

        v = domain.quantities['yvelocity']

        h = domain.quantities['height']

        if self.form == 'momentum':
            # Edge coefficient nu h; a cell with no water takes no part.
            wet = h.centroid_values > domain.minimum_allowed_height
            a_c = num.where(wet, d.centroid_values * h.centroid_values, 0.0)
            a_b = num.ascontiguousarray(d.boundary_values * h.boundary_values)
            kinematic_viscosity_operator_ext.update_elliptic_matrix(self, a_c, a_b)
            mass = h.centroid_values
            close_dry = True
        else:
            #Update operator using current height
            self.update_elliptic_matrix(d)
            mass = None
            close_dry = False

        h_c = h.centroid_values
        for vel, mom, name in ((u, 'xmomentum', 'u_stats'), (v, 'ymomentum', 'v_stats')):
            vel_old = vel.centroid_values.copy()
            vel, stats = self.parabolic_solve(vel, vel, d, u_out=vel, update_matrix=False,
                                              output_stats=True, mass=mass,
                                              close_dry=close_dry)
            setattr(self, name, stats)
            if self.form == 'momentum':
                # The change in momentum is the stress the solved velocity
                # gives across each edge, so the stresses on the two sides
                # of an edge cancel and only the boundary changes the total,
                # whatever the solver tolerance. The stresses on full cells
                # need the solved velocities of their ghost neighbours.
                self._exchange_ghost_vector(vel.centroid_values)
                stress = (self._parabolic_rhs(vel.centroid_values, vel.boundary_values)
                          - self._P * vel.centroid_values)
                domain.quantities[mom].centroid_values[:] += stress / self._areas()
            else:
                # By the change in velocity, so that cells the solve leaves
                # alone keep their momentum exactly.
                domain.quantities[mom].centroid_values[:] += h_c * (vel.centroid_values - vel_old)
        domain.distribute_to_vertices_and_edges()


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
        calculates ( [ I 0 ; 0  0] - dt T^{-1} [A B] ) [array_in ; 0]

        with T the diagonal of triangle areas (if apply_triangle_areas).
        """

        n = self.n
        tot_len = self.tot_len

        V = num.zeros((tot_len,), float)

        assert len(array_in) == n

        if array_out is None:
            array_out = num.zeros_like(array_in)

        V[0:n] = array_in
        V[n:] = 0.0

        AV = self.elliptic_matrix * V

        if self.apply_triangle_areas:
            AV = AV / self.mesh.areas

        array_out[:] = array_in - self.dt * AV

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
        except (ValueError, TypeError):
            msg = 'Sparse matrix can only "right-multiply" onto a scalar'
            raise TypeError(msg)
        else:
            new = self.elliptic_matrix * other
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


    # ------------------------------------------------------------------
    # Distributed (MPI) CG helpers
    # ------------------------------------------------------------------

    def _exchange_ghost_vector(self, v):
        """In-place ghost exchange of vector v (length n = n_full + n_ghost).

        Full-triangle values at v[Idf] are sent to neighbouring ranks and
        placed into v[Idg] (ghost indices) on each receiving rank.
        Uses MPI tag 198 to avoid collisions with the existing tag-123 path.
        """
        domain = self.domain
        if domain.numproc == 1:
            return

        import mpi4py
        import anuga.utilities.parallel_abstraction as pypar

        recv_bufs = {}
        recv_reqs  = []

        for recv_proc, recv_entry in domain.ghost_recv_dict.items():
            Idg = recv_entry[0]
            buf = num.empty(len(Idg), dtype=float)
            recv_bufs[recv_proc] = (Idg, buf)
            recv_reqs.append(pypar.comm.Irecv(buf, recv_proc, tag=198))

        send_bufs = []
        send_reqs = []

        for send_proc, send_entry in domain.full_send_dict.items():
            Idf = send_entry[0]
            buf = v[Idf].copy()   # contiguous copy — MPI must not alias with recv
            send_bufs.append(buf)
            send_reqs.append(pypar.comm.Isend(buf, send_proc, tag=198))

        mpi4py.MPI.Request.Waitall(recv_reqs + send_reqs)

        for Idg, buf in recv_bufs.values():
            v[Idg] = buf


    def _distributed_dot(self, a, b):
        """Global dot product: sum of local dot(a, b) across all ranks."""
        local = num.array([num.dot(a, b)])
        if self.domain.numproc == 1:
            return local[0]
        import anuga.utilities.parallel_abstraction as pypar
        from mpi4py.MPI import SUM
        result = num.zeros(1, dtype=float)
        pypar.comm.Allreduce(local, result, op=SUM)
        return result[0]


    def _parabolic_matvec_distributed(self, d_full):
        """Compute P * d_full for the n_full owned triangles.

        d_full : 1-D array of length n_full (search direction on full cells).

        Returns array of length n_full:
            result = (P * d_extended)[:n_full]

        where P is the matrix from _build_parabolic_csr and
        d_extended = [d_full; d_ghost] after ghost exchange.
        """
        n      = self.n        # full + ghost
        n_full = self.n_full

        # Reuse pre-allocated workspace
        V = self._kv_V
        V[:n_full]  = d_full
        V[n_full:n] = 0.0       # ghost portion — will be filled by exchange

        # Fill ghost values of the search direction
        self._exchange_ghost_vector(V)  # operates on V[0:n]

        PV = self._P * V[:n]

        return PV[:n_full]


    def _parabolic_solve_distributed(self, rhs_full, x0_full, imax, tol, atol):
        """Distributed CG: solve P x = rhs for x on n_full owned triangles.

        Dot products are reduced globally via Allreduce.
        Each SpMV exchanges ghost values of the search direction so that
        off-diagonal entries coupling full triangles to ghost triangles are
        correct.

        Returns x of length n_full.
        """
        n_full = self.n_full

        x = x0_full.copy()
        r = rhs_full - self._parabolic_matvec_distributed(x)
        d = r.copy()

        rTr  = self._distributed_dot(r, r)
        rTr0 = rTr

        if rTr0 == 0.0:
            return x

        for i in range(1, imax + 1):
            q     = self._parabolic_matvec_distributed(d)
            dTq   = self._distributed_dot(d, q)
            alpha = rTr / dTq

            x = x + alpha * d
            r = r - alpha * q

            rTr_new = self._distributed_dot(r, r)

            if rTr_new <= tol**2 * rTr0 or rTr_new <= atol**2:
                break

            bt = rTr_new / rTr
            d  = r + bt * d
            rTr = rTr_new
        else:
            raise ConvergenceError(
                'parabolic_solve (distributed): CG did not converge '
                f'after {imax} iterations; rTr={rTr_new:.3e} rTr0={rTr0:.3e}')

        return x


    # ------------------------------------------------------------------
    # _build_parabolic_csr — used only in the serial (single-rank) path
    # ------------------------------------------------------------------

    def _build_parabolic_csr(self, mass=None, close_dry=False):
        """Return P = T M - dt * A_{nn} as a Sparse_CSR (n x n).

        T is the diagonal of triangle areas (the identity if
        apply_triangle_areas is False) and M = diag(mass) (the identity if
        mass is None). The equation for row i is

            T_i M_i u_i - dt sum_j A_ij u_j = T_i M_i b_i + (boundary terms)

        A is the elliptic matrix, which is symmetric apart from the rows it
        leaves empty (cells whose diffusivity is zero): such a cell keeps
        its value, u_i = b_i. Its column is taken out of the other rows,
        either moved to the right-hand side with the cell's value
        (close_dry=False) or dropped together with its share of the diagonal,
        so that the edge carries nothing (close_dry=True). P is then
        symmetric and positive definite.

        Sets self._P, and the edges for the right-hand side: self._rhs_bdry
        (rows, boundary indices, values) and self._rhs_fixed (rows, columns,
        values).
        """
        n   = self.n
        dt  = self.dt

        data   = self.operator_data.reshape(n, 4)
        colind = self.operator_colind.reshape(n, 4)

        empty = num.all(data == 0.0, axis=1)

        rows = num.repeat(num.arange(n, dtype=num.int64), 4)
        cols = colind.ravel().astype(num.int64)
        vals = data.ravel().copy()

        diagonal = cols == rows
        interior = cols < n
        to_empty = interior & ~diagonal & empty[num.minimum(cols, n - 1)]

        if close_dry:
            # The diagonal is minus the sum of the off-diagonals: take the
            # dropped edges out of it too.
            vals[diagonal] += num.bincount(rows[to_empty], vals[to_empty], minlength=n)
            self._rhs_fixed = (rows[:0], cols[:0], vals[:0])
        else:
            self._rhs_fixed = (rows[to_empty], cols[to_empty], vals[to_empty])

        bdry = ~interior
        self._rhs_bdry = (rows[bdry], cols[bdry] - n, vals[bdry])

        keep = interior & ~to_empty
        i_rows = rows[keep]
        i_cols = cols[keep]
        i_vals = -dt * vals[keep]

        self._row_mass = self._mass(mass, empty)
        i_vals[i_rows == i_cols] += self._row_mass

        counts = num.bincount(i_rows, minlength=n).astype(num.int64)
        rowptr = num.zeros(n + 1, dtype=num.int64)
        num.cumsum(counts, out=rowptr[1:])

        self._P = Sparse_CSR(None, i_vals, i_cols, rowptr, n, n)
        return self._P


    def _areas(self):
        """The diagonal T: triangle areas, or ones if apply_triangle_areas is False."""
        if self.apply_triangle_areas:
            return self.mesh.areas
        return num.ones(self.n, dtype=float)


    def _mass(self, mass, empty):
        """The diagonal T M of the parabolic matrix. A cell the elliptic
        matrix leaves empty keeps T, so that its row stays u_i = b_i."""
        t = self._areas().copy()
        if mass is not None:
            t = num.where(empty, t, t * mass)
        return t


    def _parabolic_rhs(self, b, g):
        """T M b + dt (boundary terms and fixed empty cells) for the matrix
        from _build_parabolic_csr; b on the cells, g on the boundary."""
        n = self.n
        rhs = self._row_mass * b
        r, k, v = self._rhs_bdry
        rhs += self.dt * num.bincount(r, v * g[k], minlength=n)
        r, k, v = self._rhs_fixed
        if len(r):
            rhs += self.dt * num.bincount(r, v * b[k], minlength=n)
        return rhs


    def parabolic_solve(self, u_in, b, a = None, u_out = None, update_matrix=True, \
                       output_stats=False, use_dt_tol=True, iprint=None, imax=10000,
                       mass=None, close_dry=False):
        """
        Solve for u in the equation

        ( M - dt div a grad ) u = M b

        u | boundary = g

        M is the cell weight mass (1 if None; the depth h for the momentum
        form). See _build_parabolic_csr for cells with zero diffusivity and
        close_dry.


        u_in, u_out, f anf g are Quantity objects

        Dirichlet BC g encoded into u_in boundary_values

        Initial guess for iterative scheme is given by
        centroid values of u_in

        Centroid values of a and b provide diffusivity and rhs

        Solution u is retruned in u_out

        """

        if use_dt_tol:
            tol  = min(self.dt, 1.0e-5)
            atol = min(self.dt, 1.0e-5)
        else:
            tol  = 1.0e-5
            atol = 1.0e-5

        if u_out is None:
            u_out = Quantity(self.domain)

        if update_matrix:
            self.update_elliptic_matrix(a)

        P = self._build_parabolic_csr(mass, close_dry)

        rhs = num.ascontiguousarray(self._parabolic_rhs(b.centroid_values, u_in.boundary_values))
        x0  = u_in.centroid_values.copy()

        # The system is scaled by the triangle areas (and mass): scale the
        # absolute tolerance with it. The scale is global, so that every
        # rank of a parallel run stops on the same iteration.
        m = self._row_mass[:self.n_full]
        ones = num.ones(self.n_full, dtype=float)
        atol = atol * num.sqrt(self._distributed_dot(m, m) / self._distributed_dot(ones, ones))

        if self.domain.numproc > 1:
            # ----------------------------------------------------------
            # Distributed path: CG with ghost exchange + global Allreduce
            # ----------------------------------------------------------
            n_full = self.n_full

            x = self._parabolic_solve_distributed(
                    rhs[:n_full], x0[:n_full], imax, tol, atol)

            # Only update owned (full) triangles; ghosts are refreshed by
            # the domain's normal ghost exchange in the next timestep.
            u_out.centroid_values[:n_full] = x

        else:
            # ----------------------------------------------------------
            # Serial path: call the OpenMP C CG on P
            # ----------------------------------------------------------
            precon = num.empty(self.n, dtype=float)
            jacobi_precon_c(P, precon)

            err = cg_solve_c_precon(P, x0, rhs, imax, tol, atol, 1, precon)

            if err == -1:
                raise ConvergenceError('parabolic_solve: conjugate gradient did not converge')
            elif err != 0:
                raise ConvergenceError(
                    'parabolic_solve: conjugate gradient broke down '
                    '(C error code %d): the operator matrix is singular, '
                    'indefinite, or contains NaN' % err)

            u_out.set_values(x0, location='centroids')

        u_out.set_boundary_values(u_in.boundary_values)

        if output_stats:
            stats = Stats()
            stats.iter = 0
            stats.rTr  = 0.0
            stats.rTr0 = 0.0
            stats.dx   = 0.0
            stats.x    = num.linalg.norm(u_out.centroid_values[:self.n_full])
            stats.x0   = num.linalg.norm(x0[:self.n_full])
            return u_out, stats
        else:
            return u_out
