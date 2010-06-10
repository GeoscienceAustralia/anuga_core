import numpy as num

from Scientific.IO.NetCDF import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file.sww import SWW_file, Write_sww

def sww_merge(swwfiles, output, verbose = False):
    """
        Merge a list of sww files into a single file.
    """
    
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
            description = 'merged:' + getattr(fid, 'description')          
            first_file = False
        else:
            for tri in tris:
                verts = [vertex+tri_offset for vertex in tri]
                out_tris.append(verts)

        tri_offset += fid.dimensions['number_of_points']
        
        if verbose:
            print '  new triangle index offset is ', tri_offset
            
        x.extend(list(fid.variables['x'][:]))
        y.extend(list(fid.variables['y'][:]))
      
    points = [[xx, yy] for xx, yy in zip(x, y)]    


    # NetCDF file definition
    fido = NetCDFFile(output, netcdf_mode_w)
    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
    sww.store_header(fido, times,
                             len(out_tris),
                             len(points),
                             description=description)

    sww.store_triangulation(fido,
                                    points,
                                    num.array(out_tris).astype(num.float32))


    # Get names of static quantities
    static_quantities = {}
    for name in ['elevation']:
        static_quantities[name] = fid.variables[name][:]
    
    # Store static quantities        
    self.writer.store_static_quantities(fido, **static_quantities)
                                        
    fid.close()



    domain = Domain(points, out_tris)
    domain.set_name(output)
    sww = SWW_file(domain)
    sww.store_connectivity()
    
    for _ in times:
        sww.store_timestep()
    

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular, \
                                                            rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file.sww import SWW_file

# Create shallow water domain
domain = Domain(*rectangular_cross(2, 2))
domain.set_name('test1')
sww = SWW_file(domain)
sww.store_connectivity()

domain = Domain(*rectangular(3, 3))
domain.set_name('test2')
sww = SWW_file(domain)
sww.store_connectivity()

        
sww_merge(['test1.sww', 'test2.sww'], 'test_out.sww', verbose = True)
