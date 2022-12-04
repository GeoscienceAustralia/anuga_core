

# To change this template, choose Tools | Templates
# and open the template in the editor.



from builtins import str
from builtins import range
import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function

import anuga
from math import pi, pow, sqrt
import numpy as num
from .parallel_inlet_operator import Parallel_Inlet_operator
from .parallel_structure_operator import Parallel_Structure_operator
from .parallel_boyd_box_operator import Parallel_Boyd_box_operator
from .parallel_boyd_pipe_operator import Parallel_Boyd_pipe_operator
from .parallel_weir_orifice_trapezoid_operator import Parallel_Weir_orifice_trapezoid_operator
from .parallel_internal_boundary_operator import Parallel_Internal_boundary_operator

from . import distribute, myid, numprocs, finalize
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect

import anuga.structures.boyd_box_operator
import anuga.structures.boyd_pipe_operator
import anuga.structures.inlet_operator
import anuga.structures.internal_boundary_operator

import anuga.structures.weir_orifice_trapezoid_operator

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.parallel.parallel_shallow_water import Parallel_domain

import math

"""
Factory method for Parallel Inlet_operator. All parameters are the same
as normal Inlet_Operators master_proc coordinates the allocation process,
procs contains the potential list of processors to allocate the inlet to.

Returns None for calling processors not associated with inlet. Otherwise
return an instance of Parallel_Inlet_Operator
"""

def Inlet_operator(domain,
                   region,
                   Q,
                   velocity = None,
                   zero_velocity = False,
                   default = 0.0,
                   description = None,
                   label = None,
                   logging = False,
                   master_proc = 0,
                   procs = None,
                   verbose = False):

    """Inlet Operator - add water to a domain via an inlet.
    
    :param domain: Specify domain
    :param region: Apply Inlet flow over a region (which can be a Region, Polygon or line)
    :param Q: function(t) or scalar discharge (m^3/s)
    :param velocity: Optional [u,v] to set velocity of applied discharge
    :param zero_velocity: If set to True, velocity of inlet region set to 0
    :param default: If outside time domain of the Q function, use this default discharge
    :param description: Describe the Inlet_operator
    :param label: Give Inlet_operator a label (name)
    :param verbose: Provide verbose output
    
    

    Example:

    >>> inflow_region  = anuga.Region(domain, center=[0.0,0.0], radius=1.0)
    >>> inflow = anuga.Inlet_operator(domain, inflow_region, Q = lambda t : 1 + 0.5*math.sin(t/60))

    """


    from anuga.utilities import parallel_abstraction as pypar
    # if procs is None:
    #     procs = list(range(0,pypar.size()))

    myid = pypar.rank()

    if isinstance(domain, Parallel_domain) is False and procs == None:
        procs = [myid]
    elif isinstance(domain, Parallel_domain) and procs == None:
        procs = list(range(0, pypar.size()))

    # poly can be either a line, polygon or a regions
    if isinstance(region, anuga.Region):
        pass
    else:
        region = anuga.Region(domain, poly=region, expand_polygon=True)

    alloc, inlet_master_proc, inlet_procs, enquiry_proc = allocate_inlet_procs(domain,
                                                                               region,
                                                                               master_proc = master_proc,
                                                                               procs = procs,
                                                                               verbose = verbose)



    if alloc:
        if verbose and myid == inlet_master_proc:
            print("Parallel Inlet Operator =================")
            print("Poly = " + str(region.get_type()))
            print("Master Processor is P%d" %(inlet_master_proc))
            print("Processors are P%s" %(inlet_procs))
            print("=========================================")

        return Parallel_Inlet_operator(domain,
                                       region,
                                       Q,
                                       velocity = velocity,
                                       zero_velocity = zero_velocity,
                                       default = default,
                                       description = description,
                                       label = label,
                                       logging = logging,
                                       master_proc = inlet_master_proc,
                                       procs = inlet_procs,
                                       verbose = verbose)
    else:
        return None

"""
Factory method for Parallel Boyd_box_operator. All parameters are the same
as normal Boyd_box_operator master_proc coordinates the allocation process,
procs contains the potential list of processors to allocate the inlet to.

Returns None for calling processors not associated with structure. Otherwise
return an instance of Parallel Boyd_box_operator
"""

def Boyd_box_operator(domain,
                       losses,
                       width,
                       height=None,
                       blockage=0.0,
                       barrels=1.0,
                       z1=0.0,
                       z2=0.0,
                       end_points=None,
                       exchange_lines=None,
                       enquiry_points=None,
                       invert_elevations=None,
                       apron=0.1,
                       manning=0.013,
                       enquiry_gap=0.0,
                       smoothing_timescale=0.0,
                       use_momentum_jet=True,
                       use_velocity_head=True,
                       description=None,
                       label=None,
                       structure_type='boyd_box',
                       logging=False,
                       verbose=False,
                       master_proc=0,
                       procs=None):

 
    from anuga.utilities import parallel_abstraction as pypar
    # if procs is None:
    #     procs = list(range(0,pypar.size()))

    myid = pypar.rank()

    if isinstance(domain, Parallel_domain) is False and procs == None:
        procs = [myid]
    elif isinstance(domain, Parallel_domain) and procs == None:
        procs = list(range(0, pypar.size()))

    end_points = ensure_numeric(end_points)
    exchange_lines = ensure_numeric(exchange_lines)
    enquiry_points = ensure_numeric(enquiry_points)

    if height is None:
        height = width

    diameter = None

    if apron is None:
        apron = width






    # Calculate location of inlet enquiry points and exchange lines
    if myid == master_proc:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = __process_skew_culvert(exchange_lines, end_points, enquiry_points, apron, enquiry_gap)

            for i in procs:
                if i == master_proc: continue
                pypar.send(enquiry_points_tmp, i)

        elif end_points is not None:
            exchange_lines_tmp, enquiry_points_tmp = __process_non_skew_culvert(end_points, width,
                                                                                enquiry_points, apron, enquiry_gap)
            for i in procs:
                if i == master_proc: continue
                pypar.send(exchange_lines_tmp, i)
                pypar.send(enquiry_points_tmp, i)
        else:
            raise Exception('Define either exchange_lines or end_points')

    else:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = pypar.receive(master_proc)
        elif end_points is not None:
            exchange_lines_tmp = pypar.receive(master_proc)
            enquiry_points_tmp = pypar.receive(master_proc)

    # Determine processors associated with first inlet
    line0 = exchange_lines_tmp[0]
    enquiry_point0 = enquiry_points_tmp[0]

    alloc0, inlet0_master_proc, inlet0_procs, enquiry0_proc = allocate_inlet_procs(domain, line0, enquiry_point =  enquiry_point0,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    # Determine processors associated with second inlet
    line1 = exchange_lines_tmp[1]
    enquiry_point1 = enquiry_points_tmp[1]

    alloc1, inlet1_master_proc, inlet1_procs, enquiry1_proc = allocate_inlet_procs(domain, line1, enquiry_point =  enquiry_point1,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    structure_procs = list(set(inlet0_procs + inlet1_procs))
    inlet_master_proc = [inlet0_master_proc, inlet1_master_proc]
    inlet_procs = [inlet0_procs, inlet1_procs]
    enquiry_proc = [enquiry0_proc, enquiry1_proc]

    if myid == master_proc and verbose:
        print("Parallel Boyd Box Operator =============================")
        print("Structure Master Proc is P" + str(inlet0_master_proc))
        print("Structure Procs are P" + str(structure_procs))
        print("Inlet Master Procs are P" + str(inlet_master_proc))
        print("Inlet Procs are P" + str(inlet_procs[0]) + " and " + str(inlet_procs[1]))
        print("Inlet Enquiry Procs are P" + str(enquiry_proc))
        print("Enquiry Points are " + str(enquiry_point0) + " and " + str(enquiry_point1))
        print("Inlet Exchange Lines are " + str(line0) + " and " + str(line1))
        print("========================================================")

    if alloc0 or alloc1:
        return Parallel_Boyd_box_operator(domain=domain,
                                         losses=losses,
                                         width=width,
                                         height=height,
                                         blockage=blockage,
                                         barrels=barrels,
                                         end_points=end_points,
                                         exchange_lines=exchange_lines,
                                         enquiry_points=enquiry_points,
                                         invert_elevations=invert_elevations,
                                         apron=apron,
                                         manning=manning,
                                         enquiry_gap=enquiry_gap,
                                         smoothing_timescale=smoothing_timescale,
                                         use_momentum_jet=use_momentum_jet,
                                         use_velocity_head=use_velocity_head,
                                         description=description,
                                         label=label,
                                         structure_type=structure_type,
                                         logging=logging,
                                         verbose=verbose,
                                         master_proc = inlet0_master_proc,
                                         procs = structure_procs,
                                         inlet_master_proc = inlet_master_proc,
                                         inlet_procs = inlet_procs,
                                         enquiry_proc = enquiry_proc)
    else:
        return None



"""
Factory method for Parallel Boyd_pipe_operator. All parameters are the same
as normal Boyd_pipe_operator master_proc coordinates the allocation process,
procs contains the potential list of processors to allocate the inlet to.

Returns None for calling processors not associated with structure. Otherwise
return an instance of Parallel Boyd_pipe_operator
"""

def Boyd_pipe_operator(domain,
                       losses,
                       diameter,
                       blockage=0.0,
                       barrels=1.0,
                       end_points=None,
                       exchange_lines=None,
                       enquiry_points=None,
                       invert_elevations=None,
                       apron=0.1,
                       manning=0.013,
                       enquiry_gap=0.0,
                       smoothing_timescale=0.0,
                       use_momentum_jet=True,
                       use_velocity_head=True,
                       description=None,
                       label=None,
                       structure_type='boyd_pipe',
                       logging=False,
                       verbose=False,
                       master_proc=0,
                       procs=None):

    # If not parallel domain then allocate serial Boyd box operator
    if isinstance(domain, Parallel_domain) is False:
        if verbose: print("Allocating non parallel boyd pipe operator .....")
        return anuga.structures.boyd_pipe_operator.Boyd_pipe_operator(domain=domain,
                                                                    losses=losses,
                                                                    diameter=diameter,
                                                                    blockage=blockage,
                                                                    barrels=barrels,
                                                                    end_points=end_points,
                                                                    exchange_lines=exchange_lines,
                                                                    enquiry_points=enquiry_points,
                                                                    invert_elevations=invert_elevations,
                                                                    apron=apron,
                                                                    manning=manning,
                                                                    enquiry_gap=enquiry_gap,
                                                                    smoothing_timescale=smoothing_timescale,
                                                                    use_momentum_jet=use_momentum_jet,
                                                                    use_velocity_head=use_velocity_head,
                                                                    description=description,
                                                                    label=label,
                                                                    structure_type=structure_type,
                                                                    logging=logging,
                                                                    verbose=verbose)

    from anuga.utilities import parallel_abstraction as pypar
    if procs is None:
        procs = list(range(0,pypar.size()))

    myid = pypar.rank()

    end_points = ensure_numeric(end_points)
    exchange_lines = ensure_numeric(exchange_lines)
    enquiry_points = ensure_numeric(enquiry_points)

    width = diameter

    assert diameter is not None

    if apron is None:
        apron = width

    # Calculate location of inlet enquiry points and exchange lines
    if myid == master_proc:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = __process_skew_culvert(exchange_lines, end_points, enquiry_points, apron, enquiry_gap)

            for i in procs:
                if i == master_proc: continue
                pypar.send(enquiry_points_tmp, i)

        elif end_points is not None:
            exchange_lines_tmp, enquiry_points_tmp = __process_non_skew_culvert(end_points, width,
                                                                                enquiry_points, apron, enquiry_gap)
            for i in procs:
                if i == master_proc: continue
                pypar.send(exchange_lines_tmp, i)
                pypar.send(enquiry_points_tmp, i)
        else:
            raise Exception('Define either exchange_lines or end_points')

    else:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = pypar.receive(master_proc)
        elif end_points is not None:
            exchange_lines_tmp = pypar.receive(master_proc)
            enquiry_points_tmp = pypar.receive(master_proc)

    # Determine processors associated with first inlet
    line0 = exchange_lines_tmp[0]
    enquiry_point0 = enquiry_points_tmp[0]

    alloc0, inlet0_master_proc, inlet0_procs, enquiry0_proc = allocate_inlet_procs(domain, line0, enquiry_point =  enquiry_point0,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    # Determine processors associated with second inlet
    line1 = exchange_lines_tmp[1]
    enquiry_point1 = enquiry_points_tmp[1]

    alloc1, inlet1_master_proc, inlet1_procs, enquiry1_proc = allocate_inlet_procs(domain, line1, enquiry_point =  enquiry_point1,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    structure_procs = list(set(inlet0_procs + inlet1_procs))
    inlet_master_proc = [inlet0_master_proc, inlet1_master_proc]
    inlet_procs = [inlet0_procs, inlet1_procs]
    enquiry_proc = [enquiry0_proc, enquiry1_proc]

    if myid == master_proc and verbose:
        print("Parallel Boyd Pipe Operator =============================")
        print("Structure Master Proc is P" + str(inlet0_master_proc))
        print("Structure Procs are P" + str(structure_procs))
        print("Inlet Master Procs are P" + str(inlet_master_proc))
        print("Inlet Procs are P" + str(inlet_procs[0]) + " and " + str(inlet_procs[1]))
        print("Inlet Enquiry Procs are P" + str(enquiry_proc))
        print("Enquiry Points are " + str(enquiry_point0) + " and " + str(enquiry_point1))
        print("Inlet Exchange Lines are " + str(line0) + " and " + str(line1))
        print("========================================================")

    if alloc0 or alloc1:
       return Parallel_Boyd_pipe_operator(domain=domain,
                                         losses=losses,
                                         diameter=diameter,
                                         blockage=blockage,
                                         barrels=barrels,
                                         end_points=end_points,
                                         exchange_lines=exchange_lines,
                                         enquiry_points=enquiry_points,
                                         invert_elevations=invert_elevations,
                                         apron=apron,
                                         manning=manning,
                                         enquiry_gap=enquiry_gap,
                                         smoothing_timescale=smoothing_timescale,
                                         use_momentum_jet=use_momentum_jet,
                                         use_velocity_head=use_velocity_head,
                                         description=description,
                                         label=label,
                                         structure_type=structure_type,
                                         logging=logging,
                                         verbose=verbose,
                                         master_proc = inlet0_master_proc,
                                         procs = structure_procs,
                                         inlet_master_proc = inlet_master_proc,
                                         inlet_procs = inlet_procs,
                                         enquiry_proc = enquiry_proc)
    else:
        return None




"""
Factory method for Parallel Weir_orifice_trapezoid_operator. All parameters are the same
as normal Weir_orifice_trapezoid_operator master_proc coordinates the allocation process,
procs contains the potential list of processors to allocate the inlet to.

Returns None for calling processors not associated with structure. Otherwise
return an instance of Parallel Weir_orifice_trapezoid_operator

This was addede by PM 22/10/2013
"""

def Weir_orifice_trapezoid_operator(domain,
                       losses,
                       width,
                       blockage=0.0,
                       barrels=1.0,
                       z1=None,
                       z2=None,
                       height=None,
                       end_points=None,
                       exchange_lines=None,
                       enquiry_points=None,
                       invert_elevations=None,
                       #culvert_slope=None,
                       apron=0.1,
                       manning=0.013,
                       enquiry_gap=0.0,
                       smoothing_timescale=0.0,
                       use_momentum_jet=True,
                       use_velocity_head=True,
                       description=None,
                       label=None,
                       structure_type='weir_orifice_trapezoid',
                       logging=False,
                       verbose=False,
                       master_proc=0,
                       procs=None):

    # If not parallel domain then allocate serial Weir orifice trapezoid operator
    if isinstance(domain, Parallel_domain) is False:
        if verbose: print("Allocating non parallel weir orifice trapzezoid operator .....")
        return anuga.structures.weir_orifice_trapezoid_operator.Weir_orifice_trapezoid_operator(domain=domain,
                                                                    losses=losses,
                                                                    width=width,
                                                                    height=height,
                                                                    blockage=blockage,
                                                                    barrels=barrels,
                                                                    z1=z1,
                                                                    z2=z2,
                                                                    #culvert_slope=culvert_slope,
                                                                    end_points=end_points,
                                                                    exchange_lines=exchange_lines,
                                                                    enquiry_points=enquiry_points,
                                                                    invert_elevations=invert_elevations,
                                                                    apron=apron,
                                                                    manning=manning,
                                                                    enquiry_gap=enquiry_gap,
                                                                    smoothing_timescale=smoothing_timescale,
                                                                    use_momentum_jet=use_momentum_jet,
                                                                    use_velocity_head=use_velocity_head,
                                                                    description=description,
                                                                    label=label,
                                                                    structure_type=structure_type,
                                                                    logging=logging,
                                                                    verbose=verbose)

    from anuga.utilities import parallel_abstraction as pypar
    if procs is None:
        procs = list(range(0,pypar.size()))

    myid = pypar.rank()

    end_points = ensure_numeric(end_points)
    exchange_lines = ensure_numeric(exchange_lines)
    enquiry_points = ensure_numeric(enquiry_points)

    if height is None:
        height = width

    diameter = None

    if apron is None:
        apron = width






    # Calculate location of inlet enquiry points and exchange lines
    if myid == master_proc:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = __process_skew_culvert(exchange_lines, end_points, enquiry_points, apron, enquiry_gap)

            for i in procs:
                if i == master_proc: continue
                pypar.send(enquiry_points_tmp, i)

        elif end_points is not None:
            exchange_lines_tmp, enquiry_points_tmp = __process_non_skew_culvert(end_points, width,
                                                                                enquiry_points, apron, enquiry_gap)
            for i in procs:
                if i == master_proc: continue
                pypar.send(exchange_lines_tmp, i)
                pypar.send(enquiry_points_tmp, i)
        else:
            raise Exception('Define either exchange_lines or end_points')

    else:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = pypar.receive(master_proc)
        elif end_points is not None:
            exchange_lines_tmp = pypar.receive(master_proc)
            enquiry_points_tmp = pypar.receive(master_proc)

    # Determine processors associated with first inlet
    line0 = exchange_lines_tmp[0]
    enquiry_point0 = enquiry_points_tmp[0]

    alloc0, inlet0_master_proc, inlet0_procs, enquiry0_proc = allocate_inlet_procs(domain, line0, enquiry_point =  enquiry_point0,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    # Determine processors associated with second inlet
    line1 = exchange_lines_tmp[1]
    enquiry_point1 = enquiry_points_tmp[1]

    alloc1, inlet1_master_proc, inlet1_procs, enquiry1_proc = allocate_inlet_procs(domain, line1, enquiry_point =  enquiry_point1,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    structure_procs = list(set(inlet0_procs + inlet1_procs))
    inlet_master_proc = [inlet0_master_proc, inlet1_master_proc]
    inlet_procs = [inlet0_procs, inlet1_procs]
    enquiry_proc = [enquiry0_proc, enquiry1_proc]

    if myid == master_proc and verbose:
        print("Parallel Weir Orifice Trapezoid Operator =============================")
        print("Structure Master Proc is P" + str(inlet0_master_proc))
        print("Structure Procs are P" + str(structure_procs))
        print("Inlet Master Procs are P" + str(inlet_master_proc))
        print("Inlet Procs are P" + str(inlet_procs[0]) + " and " + str(inlet_procs[1]))
        print("Inlet Enquiry Procs are P" + str(enquiry_proc))
        print("Enquiry Points are " + str(enquiry_point0) + " and " + str(enquiry_point1))
        print("Inlet Exchange Lines are " + str(line0) + " and " + str(line1))
        print("========================================================")

    if alloc0 or alloc1:
        return Parallel_Weir_orifice_trapezoid_operator(domain=domain,
                                         losses=losses,
                                         width=width,
                                         height=height,
                                         blockage=blockage,
                                         barrels=barrels,
                                         z1=z1,
                                         z2=z2,
                                         #culvert_slope=culvert_slope,
                                         end_points=end_points,
                                         exchange_lines=exchange_lines,
                                         enquiry_points=enquiry_points,
                                         invert_elevations=invert_elevations,
                                         apron=apron,
                                         manning=manning,
                                         enquiry_gap=enquiry_gap,
                                         smoothing_timescale=smoothing_timescale,
                                         use_momentum_jet=use_momentum_jet,
                                         use_velocity_head=use_velocity_head,
                                         description=description,
                                         label=label,
                                         structure_type=structure_type,
                                         logging=logging,
                                         verbose=verbose,
                                         master_proc = inlet0_master_proc,
                                         procs = structure_procs,
                                         inlet_master_proc = inlet_master_proc,
                                         inlet_procs = inlet_procs,
                                         enquiry_proc = enquiry_proc)
    else:
        return None


"""
Factory method for Parallel Internal_boundary_operator. All parameters are the same
as normal Internal_boundary_operator. master_proc coordinates the allocation process,
procs contains the potential list of processors to allocate the inlet to.

Returns None for calling processors not associated with structure. Otherwise
return an instance of Parallel Internal_boundary_operator
"""

def Internal_boundary_operator(domain,
                               internal_boundary_function,
                               width=1.,
                               height=1.,
                               end_points=None,
                               exchange_lines=None,
                               enquiry_points=None,
                               invert_elevation=None,
                               apron=0.0,
                               enquiry_gap=0.0,
                               use_velocity_head=False,
                               zero_outflow_momentum=False,
                               force_constant_inlet_elevations=True,
                               smoothing_timescale=0.0,
                               compute_discharge_implicitly=True,
                               description=None,
                               label=None,
                               structure_type='internal_boundary',
                               logging=False,
                               verbose=True,
                               master_proc = 0,
                               procs = None,
                               inlet_master_proc = [0,0],
                               inlet_procs = None,
                               enquiry_proc = [0,0]):

    # If not parallel domain then allocate serial Internal boundary operator
    if isinstance(domain, Parallel_domain) is False:
        if verbose: print("Allocating non parallel internal_boundary operator .....")
        return anuga.structures.internal_boundary_operator.Internal_boundary_operator(domain=domain,
                                                                    internal_boundary_function=internal_boundary_function,
                                                                    width=width,
                                                                    height=height,
                                                                    end_points=end_points,
                                                                    exchange_lines=exchange_lines,
                                                                    enquiry_points=enquiry_points,
                                                                    invert_elevation=invert_elevation,
                                                                    apron=apron,
                                                                    enquiry_gap=enquiry_gap,
                                                                    use_velocity_head=use_velocity_head,
                                                                    zero_outflow_momentum=zero_outflow_momentum,
                                                                    force_constant_inlet_elevations=force_constant_inlet_elevations,
                                                                    smoothing_timescale=smoothing_timescale,
                                                                    compute_discharge_implicitly=compute_discharge_implicitly,
                                                                    description=description,
                                                                    label=label,
                                                                    structure_type=structure_type,
                                                                    logging=logging,
                                                                    verbose=verbose)

    from anuga.utilities import parallel_abstraction as pypar
    if procs is None:
        procs = list(range(0,pypar.size()))

    myid = pypar.rank()

    end_points = ensure_numeric(end_points)
    exchange_lines = ensure_numeric(exchange_lines)
    enquiry_points = ensure_numeric(enquiry_points)

    if height is None:
        height = width

    diameter = None

    if apron is None:
        apron = width


    # Calculate location of inlet enquiry points and exchange lines
    if myid == master_proc:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = __process_skew_culvert(exchange_lines, end_points, enquiry_points, apron, enquiry_gap)

            for i in procs:
                if i == master_proc: continue
                pypar.send(enquiry_points_tmp, i)

        elif end_points is not None:
            exchange_lines_tmp, enquiry_points_tmp = __process_non_skew_culvert(end_points, width,
                                                                                enquiry_points, apron, enquiry_gap)
            for i in procs:
                if i == master_proc: continue
                pypar.send(exchange_lines_tmp, i)
                pypar.send(enquiry_points_tmp, i)
        else:
            raise Exception('Define either exchange_lines or end_points')

    else:
        if exchange_lines is not None:
            exchange_lines_tmp = exchange_lines
            enquiry_points_tmp = pypar.receive(master_proc)
        elif end_points is not None:
            exchange_lines_tmp = pypar.receive(master_proc)
            enquiry_points_tmp = pypar.receive(master_proc)

    # Determine processors associated with first inlet
    line0 = exchange_lines_tmp[0]
    enquiry_point0 = enquiry_points_tmp[0]

    alloc0, inlet0_master_proc, inlet0_procs, enquiry0_proc = allocate_inlet_procs(domain, line0, enquiry_point =  enquiry_point0,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    # Determine processors associated with second inlet
    line1 = exchange_lines_tmp[1]
    enquiry_point1 = enquiry_points_tmp[1]

    alloc1, inlet1_master_proc, inlet1_procs, enquiry1_proc = allocate_inlet_procs(domain, line1, enquiry_point =  enquiry_point1,
                                                                                   master_proc = master_proc,
                                                                                   procs = procs, verbose=verbose)

    structure_procs = list(set(inlet0_procs + inlet1_procs))
    inlet_master_proc = [inlet0_master_proc, inlet1_master_proc]
    inlet_procs = [inlet0_procs, inlet1_procs]
    enquiry_proc = [enquiry0_proc, enquiry1_proc]

    if myid == master_proc and verbose:
        print("Parallel Internal boundary Operator =============================")
        print("Structure Master Proc is P" + str(inlet0_master_proc))
        print("Structure Procs are P" + str(structure_procs))
        print("Inlet Master Procs are P" + str(inlet_master_proc))
        print("Inlet Procs are P" + str(inlet_procs[0]) + " and " + str(inlet_procs[1]))
        print("Inlet Enquiry Procs are P" + str(enquiry_proc))
        print("Enquiry Points are " + str(enquiry_point0) + " and " + str(enquiry_point1))
        print("Inlet Exchange Lines are " + str(line0) + " and " + str(line1))
        print("========================================================")

    if alloc0 or alloc1:
        return Parallel_Internal_boundary_operator(domain=domain,
                                         internal_boundary_function=internal_boundary_function,
                                         width=width,
                                         height=height,
                                         end_points=end_points,
                                         exchange_lines=exchange_lines,
                                         enquiry_points=enquiry_points,
                                         invert_elevation=invert_elevation,
                                         apron=apron,
                                         enquiry_gap=enquiry_gap,
                                         use_velocity_head=use_velocity_head,
                                         zero_outflow_momentum=zero_outflow_momentum,
                                         force_constant_inlet_elevations=force_constant_inlet_elevations,
                                         smoothing_timescale=smoothing_timescale,
                                         compute_discharge_implicitly=compute_discharge_implicitly,
                                         description=description,
                                         label=label,
                                         structure_type=structure_type,
                                         logging=logging,
                                         verbose=verbose,
                                         master_proc = inlet0_master_proc,
                                         procs = structure_procs,
                                         inlet_master_proc = inlet_master_proc,
                                         inlet_procs = inlet_procs,
                                         enquiry_proc = enquiry_proc)
    else:
        return None




def __process_non_skew_culvert(end_points, width, enquiry_points, apron, enquiry_gap):
    """Create lines at the end of a culvert inlet and outlet.
    At either end two lines will be created; one for the actual flow to pass through and one a little further away
    for enquiring the total energy at both ends of the culvert and transferring flow.
    """

    culvert_vector = end_points[1] - end_points[0]
    culvert_length = math.sqrt(num.sum(culvert_vector**2))
    assert culvert_length > 0.0, 'The length of culvert is less than 0'

    culvert_vector /= culvert_length

    culvert_normal = num.array([-culvert_vector[1], culvert_vector[0]])  # Normal vector
    w = 0.5*width*culvert_normal # Perpendicular vector of 1/2 width

    exchange_lines = []

    # Build exchange polyline and enquiry point
    if enquiry_points is None:

        gap = (apron + enquiry_gap)*culvert_vector
        enquiry_points = []

        for i in [0, 1]:
            p0 = end_points[i] + w
            p1 = end_points[i] - w
            exchange_lines.append(num.array([p0, p1]))
            ep = end_points[i] + (2*i - 1)*gap #(2*i - 1) determines the sign of the points
            enquiry_points.append(ep)

    else:
        for i in [0, 1]:
            p0 = end_points[i] + w
            p1 = end_points[i] - w
            exchange_lines.append(num.array([p0, p1]))

    return exchange_lines, enquiry_points

def __process_skew_culvert(exchange_lines, end_points, enquiry_points, apron, enquiry_gap):

    """Compute skew culvert.
    If exchange lines are given, the enquiry points are determined. This is for enquiring
    the total energy at both ends of the culvert and transferring flow.
    """

    centre_point0 = 0.5*(exchange_lines[0][0] + exchange_lines[0][1])
    centre_point1 = 0.5*(exchange_lines[1][0] + exchange_lines[1][1])

    if end_points is None:
        culvert_vector = centre_point1 - centre_point0
    else:
        culvert_vector = end_points[1] - end_points[0]

    culvert_length = math.sqrt(num.sum(culvert_vector**2))
    assert culvert_length > 0.0, 'The length of culvert is less than 0'

    if enquiry_points is None:

        culvert_vector /= culvert_length
        gap = (apron + enquiry_gap)*culvert_vector

        enquiry_points = []

        enquiry_points.append(centre_point0 - gap)
        enquiry_points.append(centre_point1 + gap)

    return enquiry_points


def allocate_inlet_procs(domain, region, enquiry_point = None, master_proc = 0, procs = None, verbose = False):


    from anuga.utilities import parallel_abstraction as pypar
    # if procs is None:
    #     procs = list(range(0, pypar.size()))

    myid = pypar.rank()
    vertex_coordinates = domain.get_full_vertex_coordinates(absolute=True)
    domain_centroids = domain.centroid_coordinates
    size = 0
    has_enq_point = False
    numprocs = pypar.size()

    inlet_procs = []
    max_size = -1
    inlet_master_proc = -1
    inlet_enq_proc = -1

    # Calculate the number of points of the line inside full polygon

    #tri_id = line_intersect(vertex_coordinates, poly)
    if isinstance(region, anuga.Region):
        tri_id = region.get_indices()
    elif len(region) == 2: # poly is a line
        if verbose : print("======================")
        tri_id = line_intersect(vertex_coordinates, region)
    else: # region is a polygon
        if verbose : print("+++++++++++++++++++++++")
        tris_0 = line_intersect(vertex_coordinates, [region[0],region[1]])
        tris_1 = inside_polygon(domain_centroids, region)
        tri_id = num.union1d(tris_0, tris_1)


    if verbose:
        if isinstance(region, anuga.Region):
            print("P%d has %d triangles in region %s" %(myid, len(tri_id), region.get_type()))
        else:
            print("P%d has %d triangles in region %s" %(myid, len(tri_id), region))
            

    size = len(tri_id)

    if enquiry_point is not None:
        try:
            k = domain.get_triangle_containing_point(enquiry_point)

            if domain.tri_full_flag[k] == 1:
                size = size + 1
                has_enq_point = True
                if verbose: print("P%d has enq point %s" %(myid, enquiry_point))
            else:
                if verbose: print("P%d contains ghost copy of enq point %s" %(myid, enquiry_point))
                has_enq_point = False
        except:
            if verbose: print("P%d does not contain enq point %s" %(myid, enquiry_point))
            has_enq_point = False

    if myid == master_proc:
        # Recieve size of overlap from each processor
        # Initialize line_master_proc and inlet_procs

        if size > 0:
            inlet_procs = [master_proc]
            max_size = size
            inlet_master_proc = master_proc
            if has_enq_point:
                inlet_enq_proc = master_proc

        # Recieve size of overlap
        for i in procs:
            if i == master_proc: continue
            x = pypar.receive(i)
            y = pypar.receive(i)

            if x > 0:
                inlet_procs.append(i)

                # Choose inlet_master_proc as the one with the most overlap
                if x > max_size:
                    max_size = x
                    inlet_master_proc = i

                if y is True:
                    assert inlet_enq_proc == -1, "Enquiry point correspond to more than one proc"
                    inlet_enq_proc = i

        assert len(inlet_procs) > 0, "Line does not intersect any domain"
        assert inlet_master_proc >= 0, "No master processor assigned"

        if enquiry_point is not None:
            msg = "Enquiry point %s doesn't intersect mesh, maybe inside a building, try reducing enquiry_gap" % str(enquiry_point)
            if inlet_enq_proc < 0:
                raise Exception(msg)

        # Send inlet_master_proc and inlet_procs to all processors in inlet_procs
        for i in procs:
            if i != master_proc:
                pypar.send(inlet_master_proc, i)
                pypar.send(inlet_procs, i)
                pypar.send(inlet_enq_proc, i)

    else:
        pypar.send(size, master_proc)
        pypar.send(has_enq_point, master_proc)

        inlet_master_proc = pypar.receive(master_proc)
        inlet_procs = pypar.receive(master_proc)
        inlet_enq_proc = pypar.receive(master_proc)
        if has_enq_point: assert inlet_enq_proc == myid, "Enquiry found in proc, but not declared globally"

    if size > 0:
        return True, inlet_master_proc, inlet_procs, inlet_enq_proc
    else:
        return False, inlet_master_proc, inlet_procs, inlet_enq_proc


__author__="pete"
__date__ ="$06/09/2011 1:17:57 PM$"

if __name__ == "__main__":
    print("Parallel operator factory")
