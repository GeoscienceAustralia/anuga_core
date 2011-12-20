"""
    Merge a list of .sww files together into a single file.
"""

import numpy as num
from anuga.utilities.numerical_tools import ensure_numeric

from Scientific.IO.NetCDF import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, netcdf_float32, netcdf_int
from anuga.file.sww import SWW_file, Write_sww

def sww_merge(domain_global_name, np, verbose=False):

    output = domain_global_name+".sww"
    swwfiles = [ domain_global_name+"_P"+str(np)+"_"+str(v)+".sww" for v in range(np)]

    _sww_merge(swwfiles, output, verbose)


def sww_merge_parallel(domain_global_name, np, verbose=False, delete_old=False):

    output = domain_global_name+".sww"
    swwfiles = [ domain_global_name+"_P"+str(np)+"_"+str(v)+".sww" for v in range(np)]

    _sww_merge_parallel(swwfiles, output, verbose, delete_old)


def _sww_merge(swwfiles, output, verbose):
    """
        Merge a list of sww files into a single file.
        
        May be useful for parallel runs. Note that colinear points and
        edges are not merged: there will essentially be multiple meshes within
        the one sww file.
        
        The sww files to be merged must have exactly the same timesteps. Note
        that some advanced information and custom quantities may not be
        exported.
        
        swwfiles is a list of .sww files to merge.
        output is the output filename, including .sww extension.
        verbose True to log output information
    """

    if verbose:
        print "MERGING SWW Files"
        
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


            xllcorner = fid.xllcorner
            yllcorner = fid.yllcorner

            order      = fid.order
            xllcorner  = fid.xllcorner;
            yllcorner  = fid.yllcorner ;
            zone       = fid.zone;
            false_easting  = fid.false_easting;
            false_northing = fid.false_northing;
            datum      = fid.datum;
            projection = fid.projection;

            
            for quantity in static_quantities:
                out_s_quantities[quantity] = []

            # Quantities are stored as a 2D array of timesteps x data.
            for quantity in dynamic_quantities:
                out_d_quantities[quantity] = [ [] for _ in range(len(times))]
                 
            description = 'merged:' + getattr(fid, 'description')          
            first_file = False
        else:
            for tri in tris:
                # Advance new tri indices to point at newly appended points.
                verts = [vertex+tri_offset for vertex in tri]
                out_tris.append(verts)

        num_pts = fid.dimensions['number_of_points']
        tri_offset += num_pts
        
        if verbose:
            print '  new triangle index offset is ', tri_offset
            
        x.extend(list(fid.variables['x'][:]))
        y.extend(list(fid.variables['y'][:]))
        
        # Grow the list of static quantities associated with the x,y points
        for quantity in static_quantities:
            out_s_quantities[quantity].extend(fid.variables[quantity][:])
            
        #Collate all dynamic quantities according to their timestep
        for quantity in dynamic_quantities:
            time_chunks = fid.variables[quantity][:]
            for i, time_chunk in enumerate(time_chunks):
                out_d_quantities[quantity][i].extend(time_chunk)            
    
    # Mash all points into a single big list    
    points = [[xx, yy] for xx, yy in zip(x, y)]

    points = num.asarray(points).astype(netcdf_float32)

    fid.close()

    #---------------------------
    # Write out the SWW file
    #---------------------------

    if verbose:
        print 'Writing file ', output, ':'
    fido = NetCDFFile(output, netcdf_mode_w)
    sww = Write_sww(static_quantities, dynamic_quantities)
    sww.store_header(fido, times,
                             len(out_tris),
                             len(points),
                             description=description,
                             sww_precision=netcdf_float32)



    from anuga.coordinate_transforms.geo_reference import Geo_reference
    geo_reference = Geo_reference()
    
    sww.store_triangulation(fido, points, out_tris, points_georeference=geo_reference)

    fido.order      = order
    fido.xllcorner  = xllcorner;
    fido.yllcorner  = yllcorner ;
    fido.zone       = zone;
    fido.false_easting  = false_easting;
    fido.false_northing = false_northing;
    fido.datum      = datum;
    fido.projection = projection;
       
    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)

    # Write out all the dynamic quantities for each timestep
    for q in dynamic_quantities:
        q_values = out_d_quantities[q]
        for i, time_slice in enumerate(q_values):
            fido.variables[q][i] = num.array(time_slice, netcdf_float32)
        
        # This updates the _range values
        q_range = fido.variables[q + Write_sww.RANGE][:]
        q_values_min = num.min(q_values)
        if q_values_min < q_range[0]:
            fido.variables[q + Write_sww.RANGE][0] = q_values_min
        q_values_max = num.max(q_values)
        if q_values_max > q_range[1]:
            fido.variables[q + Write_sww.RANGE][1] = q_values_max        

                                        
    fido.close()


def _sww_merge_parallel(swwfiles, output,  verbose, delete_old):
    """
        Merge a list of sww files into a single file.
        
        Use to merge files created by parallel runs.

        The sww files to be merged must have exactly the same timesteps.

        Note that some advanced information and custom quantities may not be
        exported.
        
        swwfiles is a list of .sww files to merge.
        output is the output filename, including .sww extension.
        verbose True to log output information
    """

    if verbose:
        print "MERGING SWW Files"
        
    static_quantities = ['elevation']
    dynamic_quantities = ['stage', 'xmomentum', 'ymomentum']
    
    first_file = True
    tri_offset = 0
    for filename in swwfiles:
        if verbose:
            print 'Reading file ', filename, ':'    
    
        fid = NetCDFFile(filename, netcdf_mode_r)
         
        if first_file:

            times    = fid.variables['time'][:]
            n_steps = len(times)
            number_of_timesteps = fid.dimensions['number_of_timesteps']
            starttime = int(fid.starttime)
            
            out_s_quantities = {}
            out_d_quantities = {}


            xllcorner = fid.xllcorner
            yllcorner = fid.yllcorner

            number_of_global_triangles = int(fid.number_of_global_triangles)
            number_of_global_nodes     = int(fid.number_of_global_nodes)

            order      = fid.order
            xllcorner  = fid.xllcorner;
            yllcorner  = fid.yllcorner ;
            zone       = fid.zone;
            false_easting  = fid.false_easting;
            false_northing = fid.false_northing;
            datum      = fid.datum;
            projection = fid.projection;

            g_volumes = num.zeros((number_of_global_triangles,3),num.int)
            g_x = num.zeros((number_of_global_nodes,),num.float32)
            g_y = num.zeros((number_of_global_nodes,),num.float32)

            g_points = num.zeros((number_of_global_nodes,2),num.float32)

            for quantity in static_quantities:
                out_s_quantities[quantity] = num.zeros((number_of_global_nodes,),num.float32)

            # Quantities are stored as a 2D array of timesteps x data.
            for quantity in dynamic_quantities:
                out_d_quantities[quantity] = \
                      num.zeros((n_steps,number_of_global_nodes),num.float32)
                 
            description = 'merged:' + getattr(fid, 'description')          
            first_file = False


        # Read in from files and add to global arrays

        tri_l2g  = fid.variables['tri_l2g'][:]
        node_l2g = fid.variables['node_l2g'][:]
        tri_full_flag = fid.variables['tri_full_flag'][:]
        volumes = num.array(fid.variables['volumes'][:],dtype=num.int)
        l_volumes = num.zeros_like(volumes)


        # Change the local node ids to global id in the
        # volume array
 
        for i in range(len(l_volumes)):
            g_n0 = node_l2g[volumes[i,0]]
            g_n1 = node_l2g[volumes[i,1]]
            g_n2 = node_l2g[volumes[i,2]]
        
            l_volumes[i,:] = [g_n0,g_n1,g_n2]

        # Just pick out the full triangles
        ftri_l2g = num.compress(tri_full_flag, tri_l2g)

        #print l_volumes
        #print tri_full_flag
        #print tri_l2g
        #print ftri_l2g
    
        g_volumes[ftri_l2g] = num.compress(tri_full_flag,l_volumes,axis=0)




        #g_x[node_l2g] = fid.variables['x']
        #g_y[node_l2g] = fid.variables['y']

        g_points[node_l2g,0] = fid.variables['x']
        g_points[node_l2g,1] = fid.variables['y']
        

        #print number_of_timesteps
        
        # Read in static quantities
        for quantity in static_quantities:
            out_s_quantities[quantity][node_l2g] = \
                         num.array(fid.variables[quantity],dtype=num.float32)

        
        #Collate all dynamic quantities according to their timestep
        for quantity in dynamic_quantities:
            q = fid.variables[quantity]
            #print q.shape
            for i in range(n_steps):
                out_d_quantities[quantity][i][node_l2g] = \
                           num.array(q[i],dtype=num.float32)




        fid.close()


    #---------------------------
    # Write out the SWW file
    #---------------------------
    #print g_points.shape

    #print number_of_global_triangles
    #print number_of_global_nodes


    if verbose:
            print 'Writing file ', output, ':'
    fido = NetCDFFile(output, netcdf_mode_w)
    sww = Write_sww(static_quantities, dynamic_quantities)
    sww.store_header(fido, starttime,
                             number_of_global_triangles,
                             number_of_global_nodes,
                             description=description,
                             sww_precision=netcdf_float32)



    from anuga.coordinate_transforms.geo_reference import Geo_reference
    geo_reference = Geo_reference()
    
    sww.store_triangulation(fido, g_points, g_volumes, points_georeference=geo_reference)

    fido.order      = order
    fido.xllcorner  = xllcorner;
    fido.yllcorner  = yllcorner ;
    fido.zone       = zone;
    fido.false_easting  = false_easting;
    fido.false_northing = false_northing;
    fido.datum      = datum;
    fido.projection = projection;
       
    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)

    # Write out all the dynamic quantities for each timestep

    for i in range(n_steps):
        fido.variables['time'][i] = times[i]
        
    for q in dynamic_quantities:
        q_values = out_d_quantities[q]
        for i in range(n_steps):
            fido.variables[q][i] = q_values[i]
        
        # This updates the _range values
        q_range = fido.variables[q + Write_sww.RANGE][:]
        q_values_min = num.min(q_values)
        if q_values_min < q_range[0]:
            fido.variables[q + Write_sww.RANGE][0] = q_values_min
        q_values_max = num.max(q_values)
        if q_values_max > q_range[1]:
            fido.variables[q + Write_sww.RANGE][1] = q_values_max        

                                        
    #print out_s_quantities
    #print out_d_quantities
    
    #print g_x
    #print g_y

    #print g_volumes
    
    fido.close()
    
    if delete_old:
        import os
        for filename in swwfiles:

            if verbose:
                print 'Deleting file ', filename, ':'
            os.remove(filename)

if __name__ == "__main__":

    import argparse
    from anuga.anuga_exceptions import ANUGAError


    parser = argparse.ArgumentParser(description='Merge sww files created from parallel run')
    parser.add_argument('-np', type=int, default = 4,
                   help='number of processors used to produce sww files')
    parser.add_argument('-f', type=str, default="domain",
                   help='domain global name')
    parser.add_argument('-v', nargs='?', type=bool, const=True, default=False,
                   help='verbosity')
    parser.add_argument('-delete_old', nargs='?', type=bool, const=True, default=False,
                   help='Flag to delete the input files')
    args = parser.parse_args()

    np = args.np
    domain_global_name = args.f
    verbose = args.v
    delete_old = args.delete_old


    try:
        sww_merge_parallel(domain_global_name, np, verbose, delete_old)
    except:
        msg = 'ERROR: When merging sww files %s '% domain_global_name
        print msg
        raise
