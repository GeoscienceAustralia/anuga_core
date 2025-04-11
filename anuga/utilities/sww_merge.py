
"""
    Merge a list of .sww files together into a single file.
"""


from builtins import zip
from builtins import str
from builtins import range
import numpy as num
from anuga.utilities.numerical_tools import ensure_numeric

from anuga.file.netcdf import NetCDFFile
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

    fid = NetCDFFile(swwfiles[0], netcdf_mode_r)

    try: # works with netcdf4
        number_of_volumes = len(fid.dimensions['number_of_volumes'])
        number_of_points = len(fid.dimensions['number_of_points'])
    except: # works with scientific.io.netcdf
        number_of_volumes = int(fid.dimensions['number_of_volumes'])
        number_of_points = int(fid.dimensions['number_of_points'])

    fid.close()

    if 3*number_of_volumes == number_of_points:
        _sww_merge_parallel_non_smooth(swwfiles, output, verbose, delete_old)
    else:
        _sww_merge_parallel_smooth(swwfiles, output, verbose, delete_old)
        

def _sww_merge(swwfiles, output, verbose=False):
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
        print("MERGING SWW Files")
        
    static_quantities = ['elevation']
    dynamic_quantities = ['stage', 'xmomentum', 'ymomentum']
    
    first_file = True
    tri_offset = 0
    for filename in swwfiles:
        if verbose:
            print('Reading file ', filename, ':')    
    
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



        try: # works with netcdf4
            num_pts = len(fid.dimensions['number_of_points'])
        except: # works with scientific.io.netcdf
            num_pts = int(fid.dimensions['number_of_points'])

        tri_offset += num_pts
        
        if verbose:
            print('  new triangle index offset is ', tri_offset)
            
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
        print('Writing file ', output, ':')
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


def _sww_merge_parallel_smooth(swwfiles, output,  verbose=False, delete_old=False):
    """
        Merge a list of sww files into a single file.
        
        Use to merge files created by parallel runs.

        The sww files to be merged must have exactly the same timesteps.

        It is assumed that the separate sww files have been stored in non_smooth
        format.

        Note that some advanced information and custom quantities may not be
        exported.
        
        swwfiles is a list of .sww files to merge.
        output is the output filename, including .sww extension.
        verbose True to log output information
    """

    if verbose:
        print("MERGING SWW Files")
        
    
    first_file = True
    tri_offset = 0
    for filename in swwfiles:
        if verbose:
            print('Reading file ', filename, ':')    
    
        fid = NetCDFFile(filename, netcdf_mode_r)
         
        if first_file:

            times    = fid.variables['time'][:]
            n_steps = len(times)
            #number_of_timesteps = fid.dimensions['number_of_timesteps']
            #print n_steps, number_of_timesteps
            starttime = int(fid.starttime)
            
            out_s_quantities = {}
            out_d_quantities = {}

            out_s_c_quantities = {}
            out_d_c_quantities = {}


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

            g_volumes = num.zeros((number_of_global_triangles,3),int)
            g_x = num.zeros((number_of_global_nodes,),num.float32)
            g_y = num.zeros((number_of_global_nodes,),num.float32)

            g_points = num.zeros((number_of_global_nodes,2),num.float32)

            #=====================================
            # Deal with the vertex based variables
            #=====================================
            quantities = set(['elevation', 'friction', 'stage', 'xmomentum',
                              'ymomentum', 'xvelocity', 'yvelocity', 'height'])
            variables = set(fid.variables.keys())

            quantities = list(quantities & variables)
            
            static_quantities = []
            dynamic_quantities = []

            for quantity in quantities:
                # Test if quantity is static
                if n_steps == fid.variables[quantity].shape[0]:
                    dynamic_quantities.append(quantity)
                else:
                    static_quantities.append(quantity)
                
            for quantity in static_quantities:
                out_s_quantities[quantity] = num.zeros((number_of_global_nodes,),num.float32)

            # Quantities are stored as a 2D array of timesteps x data.
            for quantity in dynamic_quantities:
                out_d_quantities[quantity] = \
                      num.zeros((n_steps,number_of_global_nodes),num.float32)

            #=======================================
            # Deal with the centroid based variables
            #=======================================
            quantities = set(['elevation_c', 'friction_c', 'stage_c', 'xmomentum_c',
                              'ymomentum_c', 'xvelocity_c', 'yvelocity_c', 'height_c'])
            variables = set(fid.variables.keys())

            quantities = list(quantities & variables)
            
            static_c_quantities = []
            dynamic_c_quantities = []

            for quantity in quantities:
                # Test if quantity is static
                if n_steps == fid.variables[quantity].shape[0]:
                    dynamic_c_quantities.append(quantity)
                else:
                    static_c_quantities.append(quantity)
                
            for quantity in static_c_quantities:
                out_s_c_quantities[quantity] = num.zeros((number_of_global_triangles,),num.float32)

            # Quantities are stored as a 2D array of timesteps x data.
            for quantity in dynamic_c_quantities:
                out_d_c_quantities[quantity] = \
                      num.zeros((n_steps,number_of_global_triangles),num.float32)
                 
            description = 'merged:' + getattr(fid, 'description')          
            first_file = False


        # Read in from files and add to global arrays

        tri_l2g  = fid.variables['tri_l2g'][:]
        node_l2g = fid.variables['node_l2g'][:]
        tri_full_flag = fid.variables['tri_full_flag'][:]
        volumes = num.array(fid.variables['volumes'][:],dtype=int)
        l_volumes = num.zeros_like(volumes)
        l_old_volumes = num.zeros_like(volumes)


        # Change the local node ids to global id in the
        # volume array

        # FIXME SR: Surely we can knock up a numpy way of doing this
        #for i in range(len(l_volumes)):
        #    g_n0 = node_l2g[volumes[i,0]]
        #    g_n1 = node_l2g[volumes[i,1]]
        #    g_n2 = node_l2g[volumes[i,2]]
        #
        #    l_old_volumes[i,:] = [g_n0,g_n1,g_n2]

        g_n0 = node_l2g[volumes[:,0]].reshape(-1,1)
        g_n1 = node_l2g[volumes[:,1]].reshape(-1,1)
        g_n2 = node_l2g[volumes[:,2]].reshape(-1,1)

        #print g_n0.shape
        l_volumes = num.hstack((g_n0,g_n1,g_n2))

        #assert num.allclose(l_volumes, l_old_volumes)

        # Just pick out the full triangles
        ftri_ids = num.where(tri_full_flag>0)
        ftri_l2g = num.compress(tri_full_flag, tri_l2g)
        
        #f_ids = num.argwhere(tri_full_flag==1).reshape(-1,)
        #f_gids = tri_l2g[f_ids]

        #print l_volumes
        #print tri_full_flag
        #print tri_l2g
        #print ftri_l2g
        
        f_volumes0 = num.compress(tri_full_flag,volumes[:,0])
        f_volumes1 = num.compress(tri_full_flag,volumes[:,1])
        f_volumes2 = num.compress(tri_full_flag,volumes[:,2])
        
        g_volumes[ftri_l2g,0] = node_l2g[f_volumes0]
        g_volumes[ftri_l2g,1] = node_l2g[f_volumes1]
        g_volumes[ftri_l2g,2] = node_l2g[f_volumes2]

        #fg_volumes = num.compress(tri_full_flag,l_volumes,axis=0)
        #g_volumes[ftri_l2g] = fg_volumes




        #g_x[node_l2g] = fid.variables['x']
        #g_y[node_l2g] = fid.variables['y']

        g_points[node_l2g,0] = fid.variables['x'][:]
        g_points[node_l2g,1] = fid.variables['y'][:]
        

        #print number_of_timesteps


        # FIXME SR: It seems that some of the "ghost" node quantity values
        # are being storded. We should only store those nodes which are associated with
        # full triangles. So we need an index array of "full" nodes, ie those in
        # full triangles

        #use numpy.compress and numpy.unique to get "full nodes

        f_volumes = num.compress(tri_full_flag,volumes,axis=0)
        fl_nodes = num.unique(f_volumes)
        f_node_l2g = node_l2g[fl_nodes]

        #print len(node_l2g)
        #print len(fl_nodes)

        # Read in static quantities
        for quantity in static_quantities:
            #out_s_quantities[quantity][node_l2g] = \
            #             num.array(fid.variables[quantity],dtype=num.float32)
            q = fid.variables[quantity]
            #print quantity, q.shape
            out_s_quantities[quantity][f_node_l2g] = \
                         num.array(q[:],dtype=num.float32)[fl_nodes]

        
        #Collate all dynamic quantities according to their timestep
        for quantity in dynamic_quantities:
            q = fid.variables[quantity]
            #print q.shape
            for i in range(n_steps):
                #out_d_quantities[quantity][i][node_l2g] = \
                #           num.array(q[i],dtype=num.float32)
                out_d_quantities[quantity][i][f_node_l2g] = \
                           num.array(q[i],dtype=num.float32)[fl_nodes]


        # Read in static c quantities
        for quantity in static_c_quantities:
            #out_s_quantities[quantity][node_l2g] = \
            #             num.array(fid.variables[quantity],dtype=num.float32)
            q = fid.variables[quantity]
            out_s_c_quantities[quantity][ftri_l2g] = \
                         num.array(q).astype(num.float32)[ftri_ids]

        
        #Collate all dynamic c quantities according to their timestep
        for quantity in dynamic_c_quantities:
            q = fid.variables[quantity]
            #print q.shape
            for i in range(n_steps):
                out_d_c_quantities[quantity][i][ftri_l2g] = \
                           num.array(q[i]).astype(num.float32)[ftri_ids]


        fid.close()


    #---------------------------
    # Write out the SWW file
    #---------------------------
    #print g_points.shape

    #print number_of_global_triangles
    #print number_of_global_nodes


    if verbose:
            print('Writing file ', output, ':')
    fido = NetCDFFile(output, netcdf_mode_w)

    sww = Write_sww(static_quantities, dynamic_quantities, static_c_quantities, dynamic_c_quantities)
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
    sww.store_static_quantities_centroid(fido, verbose=verbose, **out_s_c_quantities)

    # Write out all the dynamic quantities for each timestep

    for i in range(n_steps):
        fido.variables['time'][i] = times[i]

        
    for q in dynamic_quantities:
        q_values = out_d_quantities[q]
        if verbose:
            print('  Writing quantity: ',q)
            
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

    for q in dynamic_c_quantities:
        if verbose:
            print('  Writing quantity: ',q)
            
        q_values = out_d_c_quantities[q]
        for i in range(n_steps):
            fido.variables[q][i] = q_values[i]

                                        
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
                print('Deleting file ', filename, ':')
            os.remove(filename)


def _sww_merge_parallel_non_smooth(swwfiles, output,  verbose=False, delete_old=False):
    """
        Merge a list of sww files into a single file.

        Used to merge files created by parallel runs.

        The sww files to be merged must have exactly the same timesteps.

        It is assumed that the separate sww files have been stored in non_smooth
        format.

        Note that some advanced information and custom quantities may not be
        exported.

        swwfiles is a list of .sww files to merge.
        output is the output filename, including .sww extension.
        verbose True to log output information
    """

    if verbose:
        print("MERGING SWW Files")


    first_file = True
    tri_offset = 0
    for filename in swwfiles:
        if verbose:
            print('Reading file ', filename, ':')

        fid = NetCDFFile(filename, netcdf_mode_r)

        if first_file:

            times    = fid.variables['time'][:]
            n_steps = len(times)
            number_of_timesteps = fid.dimensions['number_of_timesteps']
            #print n_steps, number_of_timesteps
            starttime = int(fid.starttime)

            out_s_quantities = {}
            out_d_quantities = {}

            out_s_c_quantities = {}
            out_d_c_quantities = {}


            xllcorner = fid.xllcorner
            yllcorner = fid.yllcorner

            number_of_global_triangles = int(fid.number_of_global_triangles)
            number_of_global_nodes     = int(fid.number_of_global_nodes)
            number_of_global_triangle_vertices = 3*number_of_global_triangles


            order      = fid.order
            xllcorner  = fid.xllcorner;
            yllcorner  = fid.yllcorner ;
            zone       = fid.zone;
            false_easting  = fid.false_easting;
            false_northing = fid.false_northing;
            datum      = fid.datum;
            projection = fid.projection;

            g_volumes = num.arange(number_of_global_triangles*3).reshape(-1,3)



            g_x = num.zeros((number_of_global_triangle_vertices,),num.float32)
            g_y = num.zeros((number_of_global_triangle_vertices,),num.float32)

            g_points = num.zeros((number_of_global_triangle_vertices,2),num.float32)

            #=======================================
            # Deal with the vertex based variables
            #=======================================
            quantities = set(['elevation', 'friction', 'stage', 'xmomentum',
                              'ymomentum', 'xvelocity', 'yvelocity', 'height'])
            variables = set(fid.variables.keys())

            quantities = list(quantities & variables)

            static_quantities = []
            dynamic_quantities = []

            for quantity in quantities:
                # Test if elevation is static
                if n_steps == fid.variables[quantity].shape[0]:
                    dynamic_quantities.append(quantity)
                else:
                    static_quantities.append(quantity)

            # Static Quantities are stored as a 1D array
            for quantity in static_quantities:
                out_s_quantities[quantity] = num.zeros((3*number_of_global_triangles,),num.float32)

            #=======================================
            # Deal with the centroid based variables
            #=======================================
            quantities = set(['elevation_c', 'friction_c', 'stage_c', 'xmomentum_c',
                              'ymomentum_c', 'xvelocity_c', 'yvelocity_c', 'height_c'])
            variables = set(fid.variables.keys())

            quantities = list(quantities & variables)
            
            static_c_quantities = []
            dynamic_c_quantities = []

            for quantity in quantities:
                # Test if quantity is static
                if n_steps == fid.variables[quantity].shape[0]:
                    dynamic_c_quantities.append(quantity)
                else:
                    static_c_quantities.append(quantity)
                
            for quantity in static_c_quantities:
                out_s_c_quantities[quantity] = num.zeros((number_of_global_triangles,),num.float32)

            description = 'merged:' + getattr(fid, 'description')
            first_file = False


        # Read in from files and add to global arrays

        tri_l2g  = fid.variables['tri_l2g'][:]
        node_l2g = fid.variables['node_l2g'][:]
        tri_full_flag = fid.variables['tri_full_flag'][:]

        f_ids = num.argwhere(tri_full_flag==1).reshape(-1,)
        f_gids = tri_l2g[f_ids]

        g_vids = (3*f_gids.reshape(-1,1) + num.array([0,1,2])).reshape(-1,)
        l_vids = (3*f_ids.reshape(-1,1) + num.array([0,1,2])).reshape(-1,)


        l_x = num.array(fid.variables['x'][:],dtype=num.float32)
        l_y = num.array(fid.variables['y'][:],dtype=num.float32)

        
        g_x[g_vids] = l_x[l_vids]
        g_y[g_vids] = l_y[l_vids]

        g_points[g_vids,0] = g_x[g_vids]
        g_points[g_vids,1] = g_y[g_vids]


        ## Read in static quantities
        for quantity in static_quantities:
            q = fid.variables[quantity]
            out_s_quantities[quantity][g_vids] = \
                         num.array(q).astype(num.float32)[l_vids]
                         #num.array(q,dtype=num.float32)[l_vids]


        # Read in static c quantities
        for quantity in static_c_quantities:
            q = fid.variables[quantity]
            out_s_c_quantities[quantity][f_gids] = \
                         num.array(q).astype(num.float32)[f_ids]
                         #num.array(q,dtype=num.float32)[f_ids]

        
        fid.close()

    #---------------------------
    # Write out the SWW file
    #---------------------------

    if verbose:
            print('Writing file ', output, ':')

    fido = NetCDFFile(output, netcdf_mode_w)
    sww = Write_sww(static_quantities, dynamic_quantities, static_c_quantities, dynamic_c_quantities)
    sww.store_header(fido, starttime,
                             number_of_global_triangles,
                             number_of_global_triangles*3,
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
    sww.store_static_quantities_centroid(fido, verbose=verbose, **out_s_c_quantities)
    
    # Write out all the dynamic quantities for each timestep

    for i in range(n_steps):
        fido.variables['time'][i] = times[i]

    for q in (dynamic_quantities + dynamic_c_quantities):

        if verbose:
            print('  Writing quantity: ',q)
                    
        # Initialise q_values with zeros
        if q in dynamic_quantities:
            q_values = num.zeros((n_steps, 3*number_of_global_triangles), num.float32)
        elif q in dynamic_c_quantities:
            q_values = num.zeros((n_steps, number_of_global_triangles), num.float32)


        # Read the quantities one at a time, to reduce memory usage
        for filename in swwfiles:
            fid = NetCDFFile(filename, netcdf_mode_r)

            # Index information
            tri_l2g  = fid.variables['tri_l2g'][:]
            node_l2g = fid.variables['node_l2g'][:]
            tri_full_flag = fid.variables['tri_full_flag'][:]
            f_ids = num.argwhere(tri_full_flag==1).reshape(-1,)
            f_gids = tri_l2g[f_ids]
            g_vids = (3*f_gids.reshape(-1,1) + num.array([0,1,2])).reshape(-1,)
            l_vids = (3*f_ids.reshape(-1,1) + num.array([0,1,2])).reshape(-1,)
            for i in range(n_steps):
                # Different indices for vertex and centroid quantities
                if q in dynamic_quantities:
                    q_values[i][g_vids] = \
                    num.array(fid.variables[q][i], dtype=num.float32)[l_vids]
                elif q in dynamic_c_quantities:
                    q_values[i][f_gids] = \
                    num.array(fid.variables[q][i], dtype=num.float32)[f_ids]

            fid.close()

        # Write to the file
        for i in range(n_steps):
            fido.variables[q][i] = q_values[i]

        if q in dynamic_quantities:
            # This updates the _range values
            q_range = fido.variables[q + Write_sww.RANGE][:]
            q_values_min = num.min(q_values)
            if q_values_min < q_range[0]:
                fido.variables[q + Write_sww.RANGE][0] = q_values_min
            q_values_max = num.max(q_values)
            if q_values_max > q_range[1]:
                fido.variables[q + Write_sww.RANGE][1] = q_values_max

    fido.close()

    if delete_old:
        import os
        for filename in swwfiles:

            if verbose:
                print('Deleting file ', filename, ':')
            os.remove(filename)


