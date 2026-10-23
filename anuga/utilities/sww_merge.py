
"""
    Merge a list of .sww files together into a single file.
"""
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


def sww_merge_parallel(domain_global_name, np, verbose=False, delete_old=False,
                       chunk_size=None):
    """Merge parallel SWW files produced by an MPI run into a single file.

    Parameters
    ----------
    domain_global_name : str
        Base name of the domain (without .sww extension).  The per-process
        files are expected to be named
        ``<domain_global_name>_P<np>_<rank>.sww``.
    np : int
        Number of MPI processes (= number of per-process files to merge).
    verbose : bool, optional
        Print progress messages.
    delete_old : bool, optional
        Remove the per-process files after a successful merge.
    chunk_size : int or None, optional
        Number of timesteps to process in each pass through the input files
        during the dynamic-quantity merge.  Use this to bound peak memory
        when the full time series does not fit in RAM.

        * ``None`` (default) — read all timesteps at once (fastest, most RAM).
        * Positive integer — read at most *chunk_size* timesteps per pass
          (more I/O passes, but memory is bounded to
          ``chunk_size × n_global_nodes × 4 bytes`` per dynamic quantity).

        For example, ``chunk_size=100`` processes 100 timesteps at a time.
    """

    output = domain_global_name+".sww"
    swwfiles = [ domain_global_name+"_P"+str(np)+"_"+str(v)+".sww" for v in range(np)]

    fid = NetCDFFile(swwfiles[0], netcdf_mode_r)

    try: # works with netcdf4
        number_of_volumes = len(fid.dimensions['number_of_volumes'])
        number_of_points = len(fid.dimensions['number_of_points'])
    except (TypeError, AttributeError): # works with scientific.io.netcdf
        number_of_volumes = int(fid.dimensions['number_of_volumes'])
        number_of_points = int(fid.dimensions['number_of_points'])

    fid.close()

    if 3*number_of_volumes == number_of_points:
        _sww_merge_parallel_non_smooth(swwfiles, output, verbose, delete_old,
                                       chunk_size=chunk_size)
    else:
        _sww_merge_parallel_smooth(swwfiles, output, verbose, delete_old,
                                   chunk_size=chunk_size)


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

    # Accumulate per-file arrays for later concatenation
    x_parts = []
    y_parts = []
    tris_parts = []
    s_parts = {q: [] for q in static_quantities}
    d_parts = {q: [] for q in dynamic_quantities}

    for filename in swwfiles:
        if verbose:
            print('Reading file ', filename, ':')

        fid = NetCDFFile(filename, netcdf_mode_r)

        tris = num.array(fid.variables['volumes'][:], dtype=int)

        if first_file:
            times = fid.variables['time'][:]
            out_s_quantities = {}
            out_d_quantities = {}

            order          = fid.order
            xllcorner      = fid.xllcorner
            yllcorner      = fid.yllcorner
            timezone       = getattr(fid, 'timezone', 'UTC')

            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference(NetCDFObject=fid)

            description = 'merged:' + fid.description
            first_file = False

        try: # works with netcdf4
            num_pts = len(fid.dimensions['number_of_points'])
        except (TypeError, AttributeError): # works with scientific.io.netcdf
            num_pts = int(fid.dimensions['number_of_points'])

        # Offset triangle vertex indices so they point into the global array
        tris_parts.append(tris + tri_offset)
        tri_offset += num_pts

        if verbose:
            print('  new triangle index offset is ', tri_offset)

        x_parts.append(num.array(fid.variables['x'][:], dtype=num.float32))
        y_parts.append(num.array(fid.variables['y'][:], dtype=num.float32))

        for quantity in static_quantities:
            s_parts[quantity].append(num.array(fid.variables[quantity][:], dtype=num.float32))

        # Bulk read: shape (n_steps, n_pts_in_file)
        for quantity in dynamic_quantities:
            d_parts[quantity].append(num.array(fid.variables[quantity][:], dtype=num.float32))

        fid.close()

    # Concatenate all parts into global arrays
    x = num.concatenate(x_parts)
    y = num.concatenate(y_parts)
    points = num.column_stack((x, y)).astype(netcdf_float32)
    out_tris = num.concatenate(tris_parts, axis=0)

    for quantity in static_quantities:
        out_s_quantities[quantity] = num.concatenate(s_parts[quantity])

    # Concatenate along the points axis (axis=1): shape (n_steps, total_pts)
    for quantity in dynamic_quantities:
        out_d_quantities[quantity] = num.concatenate(d_parts[quantity], axis=1).astype(netcdf_float32)

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
                             sww_precision=netcdf_float32,
                             timezone=timezone)

    sww.store_triangulation(fido, points, out_tris,
                            points_georeference=Geo_reference())

    fido.order     = order
    fido.xllcorner = xllcorner
    fido.yllcorner = yllcorner
    geo_reference.write_NetCDF(fido)

    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)

    # Bulk write all dynamic quantities
    for q in dynamic_quantities:
        q_values = out_d_quantities[q]
        fido.variables[q][:] = q_values

        # Update _range values
        q_range = fido.variables[q + Write_sww.RANGE][:]
        q_values_min = num.min(q_values)
        if q_values_min < q_range[0]:
            fido.variables[q + Write_sww.RANGE][0] = q_values_min
        q_values_max = num.max(q_values)
        if q_values_max > q_range[1]:
            fido.variables[q + Write_sww.RANGE][1] = q_values_max

    fido.close()


def _sww_merge_parallel_smooth(swwfiles, output, verbose=False, delete_old=False,
                               chunk_size=None):
    """
    Merge a list of sww files into a single file.

    Used to merge files created by parallel runs stored in smooth format
    (shared nodes between triangles).

    The sww files to be merged must have exactly the same timesteps.

    swwfiles is a list of .sww files to merge.
    output is the output filename, including .sww extension.
    verbose True to log output information
    chunk_size maximum timesteps to hold in RAM at once; None means all.
    """

    if verbose:
        print("MERGING SWW Files")

    # ---------------------------------------------------------------
    # Static pass: read geometry, index arrays, and static quantities.
    # Cache per-file index arrays so the dynamic pass can reuse them
    # without re-opening files for every quantity / chunk.
    # ---------------------------------------------------------------

    first_file = True
    file_index_cache = {}   # filename -> (fl_nodes, f_node_l2g, ftri_ids, ftri_l2g)

    for filename in swwfiles:
        if verbose:
            print('Reading static data from ', filename, ':')

        fid = NetCDFFile(filename, netcdf_mode_r)

        if first_file:
            times     = fid.variables['time'][:]
            n_steps   = len(times)
            starttime = int(fid.starttime)

            out_s_quantities   = {}
            out_s_c_quantities = {}

            number_of_global_triangles = int(fid.number_of_global_triangles)
            number_of_global_nodes     = int(fid.number_of_global_nodes)

            order          = fid.order
            xllcorner      = fid.xllcorner
            yllcorner      = fid.yllcorner
            timezone       = getattr(fid, 'timezone', 'UTC')

            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference(NetCDFObject=fid)

            g_volumes = num.zeros((number_of_global_triangles, 3), int)
            g_points  = num.zeros((number_of_global_nodes, 2), num.float32)

            # Classify vertex-based quantities as static or dynamic
            candidates = set(['elevation', 'friction', 'stage', 'xmomentum',
                               'ymomentum', 'xvelocity', 'yvelocity', 'height'])
            present = set(fid.variables.keys())
            static_quantities  = []
            dynamic_quantities = []
            for q in candidates & present:
                if fid.variables[q].shape[0] == n_steps:
                    dynamic_quantities.append(q)
                else:
                    static_quantities.append(q)

            for q in static_quantities:
                out_s_quantities[q] = num.zeros((number_of_global_nodes,), num.float32)

            # Classify centroid-based quantities (any _c variable in the file)
            static_c_quantities  = []
            dynamic_c_quantities = []
            for q in sorted(v for v in present if v.endswith('_c')):
                shape = fid.variables[q].shape
                if len(shape) == 2 and shape[0] == n_steps:
                    dynamic_c_quantities.append(q)
                elif len(shape) == 1:
                    static_c_quantities.append(q)

            for q in static_c_quantities:
                out_s_c_quantities[q] = num.zeros((number_of_global_triangles,), num.float32)

            description = 'merged:' + fid.description
            first_file = False

        # --- Geometry and index arrays ---
        tri_l2g       = fid.variables['tri_l2g'][:]
        node_l2g      = fid.variables['node_l2g'][:]
        tri_full_flag = fid.variables['tri_full_flag'][:]
        volumes       = num.array(fid.variables['volumes'][:], dtype=int)

        # Full triangles only
        ftri_ids = num.where(tri_full_flag > 0)
        ftri_l2g = num.compress(tri_full_flag, tri_l2g)

        f_volumes0 = num.compress(tri_full_flag, volumes[:, 0])
        f_volumes1 = num.compress(tri_full_flag, volumes[:, 1])
        f_volumes2 = num.compress(tri_full_flag, volumes[:, 2])
        g_volumes[ftri_l2g, 0] = node_l2g[f_volumes0]
        g_volumes[ftri_l2g, 1] = node_l2g[f_volumes1]
        g_volumes[ftri_l2g, 2] = node_l2g[f_volumes2]

        g_points[node_l2g, 0] = fid.variables['x'][:]
        g_points[node_l2g, 1] = fid.variables['y'][:]

        # Full nodes: unique vertices belonging to full triangles only
        f_volumes_full = num.compress(tri_full_flag, volumes, axis=0)
        fl_nodes   = num.unique(f_volumes_full)
        f_node_l2g = node_l2g[fl_nodes]

        # Cache index arrays for the dynamic pass
        file_index_cache[filename] = (fl_nodes, f_node_l2g, ftri_ids, ftri_l2g)

        # --- Static vertex quantities ---
        for q in static_quantities:
            out_s_quantities[q][f_node_l2g] = \
                num.array(fid.variables[q][:], dtype=num.float32)[fl_nodes]

        # --- Static centroid quantities ---
        for q in static_c_quantities:
            out_s_c_quantities[q][ftri_l2g] = \
                num.array(fid.variables[q][:]).astype(num.float32)[ftri_ids]

        fid.close()

    # ---------------------------------------------------------------
    # Write static data to the output file
    # ---------------------------------------------------------------

    if verbose:
        print('Writing file ', output, ':')

    fido = NetCDFFile(output, netcdf_mode_w)
    sww  = Write_sww(static_quantities, dynamic_quantities,
                     static_c_quantities, dynamic_c_quantities)
    sww.store_header(fido, starttime,
                     number_of_global_triangles,
                     number_of_global_nodes,
                     description=description,
                     sww_precision=netcdf_float32,
                     timezone=timezone)

    sww.store_triangulation(fido, g_points, g_volumes,
                            points_georeference=Geo_reference())

    fido.order     = order
    fido.xllcorner = xllcorner
    fido.yllcorner = yllcorner
    geo_reference.write_NetCDF(fido)

    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)
    sww.store_static_quantities_centroid(fido, verbose=verbose, **out_s_c_quantities)

    fido.variables['time'][:] = times

    # ---------------------------------------------------------------
    # Chunked dynamic pass: process at most chunk_size timesteps at
    # a time so peak RAM is bounded regardless of n_steps.
    # ---------------------------------------------------------------

    _chunk = n_steps if chunk_size is None else int(chunk_size)

    # --- Dynamic vertex quantities ---
    for q in dynamic_quantities:
        if verbose:
            print('  Writing quantity: ', q)
        q_min =  num.inf
        q_max = -num.inf
        for t_start in range(0, n_steps, _chunk):
            t_end   = min(t_start + _chunk, n_steps)
            n_chunk = t_end - t_start
            q_chunk = num.zeros((n_chunk, number_of_global_nodes), num.float32)
            for filename in swwfiles:
                fl_nodes, f_node_l2g, _, _ = file_index_cache[filename]
                fid    = NetCDFFile(filename, netcdf_mode_r)
                q_data = num.array(fid.variables[q][t_start:t_end], dtype=num.float32)
                fid.close()
                q_chunk[:, f_node_l2g] = q_data[:, fl_nodes]
            fido.variables[q][t_start:t_end] = q_chunk
            q_min = min(q_min, float(num.min(q_chunk)))
            q_max = max(q_max, float(num.max(q_chunk)))

        # Update _range values
        q_range = fido.variables[q + Write_sww.RANGE][:]
        if q_min < q_range[0]:
            fido.variables[q + Write_sww.RANGE][0] = q_min
        if q_max > q_range[1]:
            fido.variables[q + Write_sww.RANGE][1] = q_max

    # --- Dynamic centroid quantities ---
    for q in dynamic_c_quantities:
        if verbose:
            print('  Writing quantity: ', q)
        for t_start in range(0, n_steps, _chunk):
            t_end   = min(t_start + _chunk, n_steps)
            n_chunk = t_end - t_start
            q_chunk = num.zeros((n_chunk, number_of_global_triangles), num.float32)
            for filename in swwfiles:
                _, _, ftri_ids, ftri_l2g = file_index_cache[filename]
                fid    = NetCDFFile(filename, netcdf_mode_r)
                q_data = num.array(fid.variables[q][t_start:t_end], dtype=num.float32)
                fid.close()
                q_chunk[:, ftri_l2g] = q_data[:, ftri_ids[0]]
            fido.variables[q][t_start:t_end] = q_chunk

    fido.close()

    if delete_old:
        import os
        for filename in swwfiles:
            if verbose:
                print('Deleting file ', filename, ':')
            os.remove(filename)


def _sww_merge_parallel_non_smooth(swwfiles, output, verbose=False, delete_old=False,
                                   chunk_size=None):
    """
    Merge a list of sww files into a single file.

    Used to merge files created by parallel runs stored in non-smooth format
    (3 × number_of_volumes == number_of_points; each triangle has its own
    vertex copies).

    The sww files to be merged must have exactly the same timesteps.

    swwfiles is a list of .sww files to merge.
    output is the output filename, including .sww extension.
    verbose True to log output information
    chunk_size maximum timesteps to hold in RAM at once; None means all.
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
            xllcorner  = fid.xllcorner
            yllcorner  = fid.yllcorner
            timezone   = getattr(fid, 'timezone', 'UTC')

            from anuga.coordinate_transforms.geo_reference import Geo_reference
            geo_reference = Geo_reference(NetCDFObject=fid)

            g_volumes = num.arange(number_of_global_triangles*3).reshape(-1,3)

            g_points = num.zeros((number_of_global_triangle_vertices,2),num.float32)

            # Cache per-file index arrays so they are not re-read for every quantity
            file_index_cache = {}

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
            # Deal with the centroid based variables (any _c variable in the file)
            #=======================================
            static_c_quantities  = []
            dynamic_c_quantities = []
            for quantity in sorted(v for v in fid.variables.keys() if v.endswith('_c')):
                shape = fid.variables[quantity].shape
                if len(shape) == 2 and shape[0] == n_steps:
                    dynamic_c_quantities.append(quantity)
                elif len(shape) == 1:
                    static_c_quantities.append(quantity)

            for quantity in static_c_quantities:
                out_s_c_quantities[quantity] = num.zeros((number_of_global_triangles,),num.float32)

            description = 'merged:' + fid.description
            first_file = False


        # Read in from files and add to global arrays

        tri_l2g       = fid.variables['tri_l2g'][:]
        tri_full_flag = fid.variables['tri_full_flag'][:]

        f_ids  = num.argwhere(tri_full_flag == 1).reshape(-1,)
        f_gids = tri_l2g[f_ids]

        g_vids = (3*f_gids.reshape(-1,1) + num.array([0,1,2])).reshape(-1,)
        l_vids = (3*f_ids.reshape(-1,1)  + num.array([0,1,2])).reshape(-1,)

        # Cache index arrays so the dynamic pass can reuse them
        file_index_cache[filename] = (f_ids, f_gids, g_vids, l_vids)

        l_x = num.array(fid.variables['x'][:], dtype=num.float32)
        l_y = num.array(fid.variables['y'][:], dtype=num.float32)

        g_points[g_vids, 0] = l_x[l_vids]
        g_points[g_vids, 1] = l_y[l_vids]

        # Read in static vertex quantities
        for quantity in static_quantities:
            out_s_quantities[quantity][g_vids] = \
                num.array(fid.variables[quantity][:]).astype(num.float32)[l_vids]

        # Read in static centroid quantities
        for quantity in static_c_quantities:
            out_s_c_quantities[quantity][f_gids] = \
                num.array(fid.variables[quantity][:]).astype(num.float32)[f_ids]

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
                             sww_precision=netcdf_float32,
                             timezone=timezone)


    sww.store_triangulation(fido, g_points, g_volumes,
                            points_georeference=Geo_reference())

    fido.order     = order
    fido.xllcorner = xllcorner
    fido.yllcorner = yllcorner
    geo_reference.write_NetCDF(fido)

    sww.store_static_quantities(fido, verbose=verbose, **out_s_quantities)
    sww.store_static_quantities_centroid(fido, verbose=verbose, **out_s_c_quantities)

    # Bulk write time
    fido.variables['time'][:] = times

    # ---------------------------------------------------------------
    # Chunked dynamic pass: process at most chunk_size timesteps at
    # a time so peak RAM is bounded regardless of n_steps.
    # ---------------------------------------------------------------

    _chunk = n_steps if chunk_size is None else int(chunk_size)

    for q in (dynamic_quantities + dynamic_c_quantities):
        if verbose:
            print('  Writing quantity: ', q)

        is_vertex_q = q in dynamic_quantities
        n_global_pts = (3 * number_of_global_triangles if is_vertex_q
                        else number_of_global_triangles)

        q_min =  num.inf
        q_max = -num.inf

        for t_start in range(0, n_steps, _chunk):
            t_end   = min(t_start + _chunk, n_steps)
            n_chunk = t_end - t_start
            q_chunk = num.zeros((n_chunk, n_global_pts), num.float32)

            # Read each file using cached index arrays; slice only the
            # current chunk of timesteps.
            for filename in swwfiles:
                f_ids, f_gids, g_vids, l_vids = file_index_cache[filename]
                fid    = NetCDFFile(filename, netcdf_mode_r)
                q_data = num.array(fid.variables[q][t_start:t_end], dtype=num.float32)
                fid.close()

                if is_vertex_q:
                    q_chunk[:, g_vids] = q_data[:, l_vids]
                else:
                    q_chunk[:, f_gids] = q_data[:, f_ids]

            fido.variables[q][t_start:t_end] = q_chunk

            if is_vertex_q:
                q_min = min(q_min, float(num.min(q_chunk)))
                q_max = max(q_max, float(num.max(q_chunk)))

        if is_vertex_q:
            # Update _range values (accumulated across all chunks)
            q_range = fido.variables[q + Write_sww.RANGE][:]
            if q_min < q_range[0]:
                fido.variables[q + Write_sww.RANGE][0] = q_min
            if q_max > q_range[1]:
                fido.variables[q + Write_sww.RANGE][1] = q_max

    fido.close()

    if delete_old:
        import os
        for filename in swwfiles:

            if verbose:
                print('Deleting file ', filename, ':')
            os.remove(filename)


