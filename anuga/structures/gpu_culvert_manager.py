"""
GPU Culvert Manager - registers Boyd box/pipe operators with the GPU
culvert system for batched execution.

Instead of each culvert doing its own GPU↔CPU sync (2 syncs × 20 culverts = 40),
this manager batches all culverts into a single gather → compute → scatter cycle
(exactly 2 GPU sync points total).

Usage:
    # After all Boyd operators are created and domain is in GPU mode:
    manager = GPUCulvertManager(domain)
    manager.register_all()

    # Then in the evolve loop, instead of per-culvert __call__:
    manager.apply_all()  # replaces domain.apply_fractional_steps() for culverts
"""

import numpy as np


CULVERT_TYPE_BOX = 0
CULVERT_TYPE_PIPE = 1
CULVERT_TYPE_WEIR_TRAPEZOID = 2


class GPUCulvertManager:
    """Manages batched GPU execution of all Boyd culvert operators."""

    def __init__(self, domain):
        self.domain = domain
        self.culvert_ids = []       # GPU-side IDs returned by init
        self.operators = []          # Python operator references (for stats)
        self._initialized = False
        self._mapped = False

    @staticmethod
    def _is_boyd_box(op):
        """Check if op is a Boyd box operator (serial or parallel)."""
        from anuga.structures.boyd_box_operator import Boyd_box_operator
        if isinstance(op, Boyd_box_operator):
            return True
        try:
            from anuga.parallel.parallel_boyd_box_operator import Parallel_Boyd_box_operator
            return isinstance(op, Parallel_Boyd_box_operator)
        except ImportError:
            return False

    @staticmethod
    def _is_boyd_pipe(op):
        """Check if op is a Boyd pipe operator (serial or parallel)."""
        from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
        if isinstance(op, Boyd_pipe_operator):
            return True
        try:
            from anuga.parallel.parallel_boyd_pipe_operator import Parallel_Boyd_pipe_operator
            return isinstance(op, Parallel_Boyd_pipe_operator)
        except ImportError:
            return False

    @staticmethod
    def _is_weir_trapezoid(op):
        """Check if op is a Weir_orifice_trapezoid_operator (serial or parallel)."""
        from anuga.structures.weir_orifice_trapezoid_operator import Weir_orifice_trapezoid_operator
        return isinstance(op, Weir_orifice_trapezoid_operator)

    @staticmethod
    def is_boyd_operator(op):
        """Check if op is any supported GPU culvert operator (box, pipe, or weir trapezoid)."""
        return (GPUCulvertManager._is_boyd_box(op)
                or GPUCulvertManager._is_boyd_pipe(op)
                or GPUCulvertManager._is_weir_trapezoid(op))

    @staticmethod
    def _is_fully_local(op):
        """Check if a parallel Boyd operator has all data local (no MPI needed).

        For serial operators, always returns True.
        For parallel operators, checks that this process owns both enquiry
        points and both inlet regions.
        """
        if not hasattr(op, 'myid'):
            return True  # Serial operator, always local

        myid = op.myid
        # Must be master proc (computes discharge)
        if myid != op.master_proc:
            return False
        # Both enquiry points must be local
        if myid != op.enquiry_proc[0] or myid != op.enquiry_proc[1]:
            return False
        # Both inlet regions must be local
        if op.inlet_master_proc[0] != myid or op.inlet_master_proc[1] != myid:
            return False
        return True

    def register_operator(self, op):
        """Register a single Boyd box/pipe operator (serial or parallel).

        Must be called before map_to_gpu().
        Cross-boundary parallel operators are registered with MPI topology
        so the C code handles the inter-rank communication.
        """
        gpu_interface = self.domain.gpu_interface
        if gpu_interface is None:
            raise RuntimeError("Domain has no GPU interface - set_multiprocessor_mode(2) first")

        # Import the Cython extension
        from anuga.shallow_water import sw_domain_gpu_ext as ext

        # Determine culvert type and geometry
        if self._is_boyd_box(op):
            culvert_type = CULVERT_TYPE_BOX
            width = op.culvert_width
            height = op.culvert_height
            diameter = 0.0
            z1 = 0.0
            z2 = 0.0
        elif self._is_boyd_pipe(op):
            culvert_type = CULVERT_TYPE_PIPE
            width = 0.0
            height = 0.0
            diameter = op.culvert_diameter
            z1 = 0.0
            z2 = 0.0
        elif self._is_weir_trapezoid(op):
            culvert_type = CULVERT_TYPE_WEIR_TRAPEZOID
            width = op.culvert_width
            height = op.culvert_height
            diameter = 0.0
            z1 = op.culvert_z1
            z2 = op.culvert_z2
        else:
            raise TypeError(f"Unsupported operator type: {type(op)}")

        is_local = self._is_fully_local(op)
        is_parallel = hasattr(op, 'myid')

        # Extract inlet indices and areas (may be empty for remote inlets)
        inlet0 = op.inlets[0]
        inlet1 = op.inlets[1]

        if inlet0 is not None and hasattr(inlet0, 'triangle_indices') and len(inlet0.triangle_indices) > 0:
            inlet0_indices = np.ascontiguousarray(inlet0.triangle_indices, dtype=np.intc)
            inlet0_areas = np.ascontiguousarray(
                self.domain.areas[inlet0.triangle_indices], dtype=np.float64)
        else:
            inlet0_indices = np.array([], dtype=np.intc)
            inlet0_areas = np.array([], dtype=np.float64)

        if inlet1 is not None and hasattr(inlet1, 'triangle_indices') and len(inlet1.triangle_indices) > 0:
            inlet1_indices = np.ascontiguousarray(inlet1.triangle_indices, dtype=np.intc)
            inlet1_areas = np.ascontiguousarray(
                self.domain.areas[inlet1.triangle_indices], dtype=np.float64)
        else:
            inlet1_indices = np.array([], dtype=np.intc)
            inlet1_areas = np.array([], dtype=np.float64)

        # Enquiry point triangle indices (-1 if not on this rank)
        if inlet0 is not None and hasattr(inlet0, 'enquiry_index') and inlet0.enquiry_index >= 0:
            enquiry_index_0 = int(inlet0.enquiry_index)
        else:
            enquiry_index_0 = -1

        if inlet1 is not None and hasattr(inlet1, 'enquiry_index') and inlet1.enquiry_index >= 0:
            enquiry_index_1 = int(inlet1.enquiry_index)
        else:
            enquiry_index_1 = -1

        # Outward culvert vectors (use operator's if inlet is None)
        if inlet0 is not None and hasattr(inlet0, 'outward_culvert_vector'):
            outward_vector_0 = inlet0.outward_culvert_vector
        else:
            outward_vector_0 = getattr(op, 'outward_vector_0', np.array([1.0, 0.0]))

        if inlet1 is not None and hasattr(inlet1, 'outward_culvert_vector'):
            outward_vector_1 = inlet1.outward_culvert_vector
        else:
            outward_vector_1 = getattr(op, 'outward_vector_1', np.array([-1.0, 0.0]))

        # Invert elevations
        if inlet0 is not None:
            has_invert_0 = 1 if inlet0.invert_elevation is not None else 0
            invert_elev_0 = float(inlet0.invert_elevation) if has_invert_0 else 0.0
        else:
            has_invert_0 = 0
            invert_elev_0 = 0.0

        if inlet1 is not None:
            has_invert_1 = 1 if inlet1.invert_elevation is not None else 0
            invert_elev_1 = float(inlet1.invert_elevation) if has_invert_1 else 0.0
        else:
            has_invert_1 = 0
            invert_elev_1 = 0.0

        # MPI topology
        if is_parallel:
            master_proc = op.master_proc
            enquiry_proc_0 = op.enquiry_proc[0]
            enquiry_proc_1 = op.enquiry_proc[1]
            inlet_master_proc_0 = op.inlet_master_proc[0]
            inlet_master_proc_1 = op.inlet_master_proc[1]
        else:
            rank = getattr(self.domain, 'processor', 0)
            master_proc = rank
            enquiry_proc_0 = rank
            enquiry_proc_1 = rank
            inlet_master_proc_0 = rank
            inlet_master_proc_1 = rank

        # Assign unique MPI tag base from operator label (deterministic across ranks)
        import hashlib
        label = getattr(op, 'label', str(len(self.culvert_ids)))
        label_hash = int(hashlib.md5(label.encode()).hexdigest()[:8], 16)
        mpi_tag_base = (label_hash % 100000) * 10 + 1000

        # Register with GPU
        culvert_id = ext.init_culvert_operator(
            gpu_interface.gpu_dom,
            culvert_type=culvert_type,
            width=width,
            height=height,
            diameter=diameter,
            z1=z1,
            z2=z2,
            length=op.culvert_length,
            manning=op.manning,
            sum_loss=op.sum_loss,
            blockage=op.culvert_blockage if hasattr(op, 'culvert_blockage') else op.blockage,
            barrels=op.culvert_barrels if hasattr(op, 'culvert_barrels') else op.barrels,
            use_velocity_head=1 if op.use_velocity_head else 0,
            use_momentum_jet=1 if op.use_momentum_jet else 0,
            use_old_momentum_method=1 if getattr(op, 'use_old_momentum_method', True) else 0,
            always_use_Q_wetdry_adjustment=1 if getattr(op, 'always_use_Q_wetdry_adjustment', True) else 0,
            max_velocity=op.max_velocity,
            smoothing_timescale=op.smoothing_timescale,
            outward_vector_0=outward_vector_0,
            outward_vector_1=outward_vector_1,
            invert_elevation_0=invert_elev_0,
            invert_elevation_1=invert_elev_1,
            has_invert_elevation_0=has_invert_0,
            has_invert_elevation_1=has_invert_1,
            enquiry_index_0=enquiry_index_0,
            enquiry_index_1=enquiry_index_1,
            inlet0_indices=inlet0_indices,
            inlet0_areas=inlet0_areas,
            inlet1_indices=inlet1_indices,
            inlet1_areas=inlet1_areas,
            master_proc=master_proc,
            enquiry_proc_0=enquiry_proc_0,
            enquiry_proc_1=enquiry_proc_1,
            inlet_master_proc_0=inlet_master_proc_0,
            inlet_master_proc_1=inlet_master_proc_1,
            is_local=1 if is_local else 0,
            mpi_tag_base=mpi_tag_base,
            init_smooth_Q=float(getattr(op, 'smooth_Q', 0.0)),
            init_smooth_delta_total_energy=float(
                getattr(op, 'smooth_delta_total_energy', 0.0)),
        )

        if culvert_id < 0:
            raise RuntimeError(f"Failed to register culvert '{op.label}' on GPU")

        self.culvert_ids.append(culvert_id)
        self.operators.append(op)

        return culvert_id

    def register_all(self):
        """Find and register all Boyd operators in domain.fractional_step_operators.

        Both fully-local and cross-boundary culverts are registered.
        Cross-boundary culverts use MPI in the C layer for data exchange.
        """
        n_local = 0
        n_parallel = 0
        # With tracers, a culvert split across ranks is left to the host MPI
        # path (Parallel_Structure_operator), which sums the tracer it carries
        # over the structure's ranks; the batched kernel only carries tracers
        # for culverts local to one rank. Unregistered, it is synced like any
        # other CPU operator.
        tracers = getattr(self.domain, 'number_of_tracers', 0) > 0
        for op in self.domain.fractional_step_operators:
            if self.is_boyd_operator(op):
                if tracers and not self._is_fully_local(op):
                    continue
                self.register_operator(op)
                if self._is_fully_local(op):
                    n_local += 1
                else:
                    n_parallel += 1

        self._initialized = True

        rank = getattr(self.domain, 'processor', 0)
        if self.culvert_ids:
            # The unified culvert kernels target a GPU only when offload is active;
            # otherwise they run as CPU-multicore OpenMP. Report where they run so
            # the log is not misleading on a CPU build / with offload disabled.
            try:
                from anuga import gpu_offload_enabled
                backend = 'GPU offload' if gpu_offload_enabled() else 'CPU multicore'
            except Exception:
                backend = 'CPU multicore'
            msg = (f"[Rank {rank}] CulvertManager: registered "
                   f"{len(self.culvert_ids)} culverts (unified C kernels, {backend})")
            if n_parallel:
                msg += f" ({n_local} local, {n_parallel} cross-boundary MPI)"
            print(msg)

    def map_to_gpu(self):
        """Map scratch buffers to GPU. Call after register_all() and GPU domain init."""
        if not self._initialized:
            raise RuntimeError("Call register_all() before map_to_gpu()")
        if self._mapped:
            return

        from anuga.shallow_water import sw_domain_gpu_ext as ext
        ext.map_culvert_operators(self.domain.gpu_interface.gpu_dom)
        self._mapped = True

    def apply_all(self):
        """Execute all culverts for the current timestep.

        This replaces the per-culvert __call__() that would normally happen
        in apply_fractional_steps(). Only 2 GPU sync points total.
        """
        if not self._mapped:
            self.map_to_gpu()

        from anuga.shallow_water import sw_domain_gpu_ext as ext
        timestep = self.domain.get_timestep()
        ext.apply_all_culvert_operators(
            self.domain.gpu_interface.gpu_dom, timestep)
        self._update_operator_stats(ext)

    def _update_operator_stats(self, ext):
        """Copy the per-culvert stats the C layer just computed onto the Python
        operator objects, so their ``.log`` files carry the same discharge /
        velocity / energy / accumulated-flow columns as mode-1.

        The C values are non-zero only on each culvert's master proc (the one
        that computes the discharge), which is also the only proc that logs, so
        writing them on every proc is harmless. We mirror the mode-1 discharge
        routine exactly: ``accumulated_flow`` and ``discharge_abs_timemean``
        accumulate every timestep (the logger resets the latter each yieldstep),
        while the rest are instantaneous.
        """
        gpu_dom = self.domain.gpu_interface.gpu_dom
        yieldstep = getattr(self.domain, 'yieldstep', None)
        for op, culvert_id in zip(self.operators, self.culvert_ids):
            gain, discharge, velocity, driving_energy, delta_total_energy = \
                ext.get_culvert_report_stats(gpu_dom, culvert_id)
            op.accumulated_flow += gain
            if yieldstep:
                op.discharge_abs_timemean += gain / yieldstep
            op.discharge = discharge
            op.velocity = velocity
            op.driving_energy = driving_energy
            op.delta_total_energy = delta_total_energy

    def finalize(self):
        """Clean up all GPU culvert resources."""
        if self._initialized:
            from anuga.shallow_water import sw_domain_gpu_ext as ext
            ext.finalize_all_culvert_operators(self.domain.gpu_interface.gpu_dom)
            self._initialized = False
            self._mapped = False
            self.culvert_ids.clear()
            self.operators.clear()

    @property
    def num_culverts(self):
        return len(self.culvert_ids)

    @property
    def is_gpu_safe(self):
        """Returns True — this manager handles culverts on GPU."""
        return True
