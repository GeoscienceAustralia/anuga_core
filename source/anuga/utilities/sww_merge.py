import numpy as num
from anuga.utilities.numerical_tools import ensure_numeric

from Scientific.IO.NetCDF import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, \
                            netcdf_mode_a, netcdf_float
from anuga.file.sww import SWW_file, Write_sww


def sww_merge(swwfiles, output, verbose = False):
    """
        Merge a list of sww files into a single file.
    """
    
    static_quantities = ['elevation']
    dynamic_quantities = ['stage', 'xmomentum', 'ymomentum']
    
    first_file = True
    tri_offset = 0
    for filename in swwfiles:
        if verbose:
            print 'Reading file ', filename, ':'    
    
        fid = NetCDFFile(filename, netcdf_mode_r)
        tris = fid.variables['volumes'][:]       
         
        if first_file:
            times = fid.variables['time'][:]
            x = []
            y = []
            out_tris = list(tris)  
            out_s_quantities = {}
            out_d_quantities = {}
            
            for quantity in static_quantities:
                out_s_quantities[quantity] = []

            for quantity in dynamic_quantities:
                out_d_quantities[quantity] = [ [] for _ in range(len(times))]
                 
            description = 'merged:' + getattr(fid, 'description')          
            first_file = False
        else:
            for tri in tris:
                verts = [vertex+tri_offset for vertex in tri]
                out_tris.append(verts)

        num_pts = fid.dimensions['number_of_points']
        tri_offset += num_pts
        
        if verbose:
            print '  new triangle index offset is ', tri_offset
            
        x.extend(list(fid.variables['x'][:]))
        y.extend(list(fid.variables['y'][:]))
        
        for quantity in static_quantities:
            out_s_quantities[quantity].extend(fid.variables[quantity][:])
            
        for quantity in dynamic_quantities:
            time_chunks = fid.variables[quantity][:]
            for i, time_chunk in enumerate(time_chunks):
                out_d_quantities[quantity][i].extend(time_chunk)            
        
    points = [[xx, yy] for xx, yy in zip(x, y)]
    fid.close()
    
    # NetCDF file definition
    fido = NetCDFFile(output, netcdf_mode_w)
    sww = Write_sww(static_quantities, dynamic_quantities)
    sww.store_header(fido, times,
                             len(out_tris),
                             len(points),
                             description=description,
                             sww_precision=netcdf_float)


    sww.store_triangulation(fido, points, out_tris)
       
    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)

    for q in dynamic_quantities:
        q_values = ensure_numeric(out_d_quantities[quantity])
        x = q_values.astype(netcdf_float)
        fido.variables[q] = x

        # This updates the _range values
        q_range = fido.variables[q + Write_sww.RANGE][:]
        q_values_min = num.min(q_values)
        if q_values_min < q_range[0]:
            fido.variables[q + Write_sww.RANGE][0] = q_values_min
        q_values_max = num.max(q_values)
        if q_values_max > q_range[1]:
            fido.variables[q + Write_sww.RANGE][1] = q_values_max        

                                        
    fido.close()
    

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular, \
                                                            rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file.sww import SWW_file
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import \
    Dirichlet_boundary

Bd = Dirichlet_boundary([0.5, 0., 0.])

# Create shallow water domain
domain = Domain(*rectangular_cross(2, 2))
domain.set_name('test1')
domain.set_quantity('elevation', 2)
domain.set_quantity('stage', 5)
domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
for t in domain.evolve(yieldstep=0.5, finaltime=1):
    pass
    
domain = Domain(*rectangular(3, 3))
domain.set_name('test2')
domain.set_quantity('elevation', 3)
domain.set_quantity('stage', 50)
domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
for t in domain.evolve(yieldstep=0.5, finaltime=1):
    pass
        
sww_merge(['test1.sww', 'test2.sww'], 'test_out.sww', verbose = True)
