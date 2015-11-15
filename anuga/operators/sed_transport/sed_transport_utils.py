"""
Contains extra utilities created while making the sediment transport, vegetation, and spatially-varied precipitation operators.

Most of these should go in other classes within ANUGA but are here to not modify the main code

M. Perignon
06/2014
"""
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.geospatial_data.geospatial_data import Geospatial_data
import os
from scipy.interpolate import NearestNDInterpolator
import numpy as np
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Boundary, File_boundary, Dirichlet_boundary
from anuga.shallow_water.boundaries import Reflective_boundary 
from shallow_water_ext import rotate
from anuga.config import epsilon, g
import sed_transport_config as st
import numpy as np

# def get_distance_centroids(self):
#     """ Calculate the distance between each centroid and the
#     centroids of its neighboring triangles
#     """
#     
#     tri_coord = self.domain.get_centroid_coordinates()
#     tri_neigh = self.domain.neighbours
#     
#     tri_neigh[tri_neigh<0]=-1
#     
#     tri_x = tri_coord[:,0]
#     tri_x = tri_x[np.newaxis, :].transpose()
#     tri_y = tri_coord[:,1]
#     tri_y = tri_y[np.newaxis, :].transpose()
#     
#     neigh_x = tri_coord[tri_neigh,0]
#     neigh_y = tri_coord[tri_neigh,1]
#     
#     c_dist = np.sqrt(np.subtract(neigh_x,tri_x)**2 + \
#                    np.subtract(neigh_y,tri_y)**2) 
#     
#     c_dist[tri_neigh<0] = 0.0
#     
#     st.distances = c_dist # for unit tests
#     
#     return c_dist

"""
create_domain_from_regions_sed

Same as create_domain_from_regions in anuga/__init__.py but accepting
conserved_quantities, evolved_quantities and other_quantities as arguments
"""
def create_domain_from_regions_sed(bounding_polygon,
                               boundary_tags,
                               maximum_triangle_area=None,
                               mesh_filename=None,
                               interior_regions=None,
                               interior_holes=None,
                               hole_tags=None,
                               poly_geo_reference=None,
                               mesh_geo_reference=None,
                               minimum_triangle_angle=28.0,
                               conserved_quantities=None,
                               evolved_quantities=None,
                               other_quantities=None,
                               fail_if_polygons_outside=True,
                               use_cache=False,
                               verbose=True):


    # Build arguments and keyword arguments for use with caching or apply.
    args = (bounding_polygon,
            boundary_tags)
    
    kwargs = {'maximum_triangle_area': maximum_triangle_area,
              'mesh_filename': mesh_filename,
              'interior_regions': interior_regions,
              'interior_holes': interior_holes,
              'hole_tags': hole_tags,
              'poly_geo_reference': poly_geo_reference,
              'mesh_geo_reference': mesh_geo_reference,
              'minimum_triangle_angle': minimum_triangle_angle,
              'conserved_quantities': conserved_quantities,
              'evolved_quantities': evolved_quantities,
              'other_quantities': other_quantities,
              'fail_if_polygons_outside': fail_if_polygons_outside,
              'verbose': verbose} #FIXME (Ole): See ticket:14

    # Call underlying engine with or without caching
    if use_cache is True:
        try:
            from anuga.caching import cache
        except:
            msg = 'Caching was requested, but caching module'+\
                  'could not be imported'
            raise (msg)


        domain = cache(_create_domain_from_regions_sed,
                       args, kwargs,
                       verbose=verbose,
                       compression=False)
    else:
        domain = apply(_create_domain_from_regions_sed,
                       args, kwargs)

    return domain

        
def _create_domain_from_regions_sed(bounding_polygon,
                                boundary_tags,
                                maximum_triangle_area=None,
                                mesh_filename=None,                           
                                interior_regions=None,
                                interior_holes=None,
                                hole_tags=None,
                                poly_geo_reference=None,
                                mesh_geo_reference=None,
                                minimum_triangle_angle=28.0,
                                conserved_quantities=None,
                                evolved_quantities=None,
                                other_quantities=None,
                                fail_if_polygons_outside=True,
                                verbose=True):
    """_create_domain_from_regions_sed - internal function.

    See create_domain_from_regions in anuga/__init__.py for documentation.
    """

    from anuga.pmesh.mesh_interface import create_mesh_from_regions
    
    create_mesh_from_regions(bounding_polygon,
                             boundary_tags,
                             maximum_triangle_area=maximum_triangle_area,
                             interior_regions=interior_regions,
                             filename=mesh_filename,
                             interior_holes=interior_holes,
                             hole_tags=hole_tags,
                             poly_geo_reference=poly_geo_reference,
                             mesh_geo_reference=mesh_geo_reference,
                             minimum_triangle_angle=minimum_triangle_angle,
                             fail_if_polygons_outside=fail_if_polygons_outside,
                             use_cache=False,
                             verbose=verbose)

    domain = Domain(mesh_filename, 
                    conserved_quantities=conserved_quantities, 
                    evolved_quantities=evolved_quantities, 
                    other_quantities=other_quantities, 
                    use_cache=False, 
                    verbose=verbose)


    return domain



def set_quantity_NNeigh(self, name,
                           filename=None):
    """Set values for named quantity from a pts file
    using nearest neighbour interpolator. The quantity
    at each point in the mesh takes on the value of the
    raster cell the point falls within (same as the value
    of the nearest point in the quantity raster pts file).
    
    This is useful for importing maps of vegetation type
    where each value of vegetation type corresponds to
    a stem spacing and diameter in a lookup table
    
    Values are set for centroids only.
    Don't want to interpolate but just pull the value
    from the quantity raster.
    """
    
    L = [filename]
    msg = ('Filename must be present')
    assert L.count(None) == len(L)-1, msg
    
    msg = 'Extension should be .pts. Use generic_asc2dem and generic_dem2pts to convert it.'
    assert os.path.splitext(filename)[1] in ['.pts'], msg

    # Assign values
    set_values_NNeigh(self, name, filename)


def set_values_NNeigh(self, name, filename):

    """ Sets the values of the quantity 'name' at centroids,
    from raster 'filename' using
    a nearest neighbour interpolation. This extracts the exact
    value of the raster at those coordinates and does not
    interpolate the values.
    Do not set at vertices or edges - not used in veg calculations
    """

#     index = self.get_unique_vertices()
#     volume_id = [i / 3 for i in index]
#     vertex_id = [i % 3 for i in index]
#     
#     print volume_id
#     print vertex_id
#     
    coord = self.get_nodes(absolute=True)
    
    # extract the data from the pts file 
    G_data = Geospatial_data(filename)
    points = G_data.get_data_points(absolute=True)
    z = G_data.get_attributes(attribute_name=None)
    
    # create interpolator
    interp = NearestNDInterpolator( points, z )
    
    # set quantity at centroids
    z = interp( coord )
    z = z[np.newaxis, :].transpose()
    
#     z_c = np.concatenate((coord, z), axis=1 )

    
    
    self.quantities[name].set_values_from_array(z, location = 'unique vertices')


####################################################################
#
# Sediment transport - specific boundary conditions
#
####################################################################

class Reflective_boundary_Sed(Reflective_boundary):
    """Reflective boundary Sed returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    
    Specific to sed transport operator - reflects concentration
    without changing the sign
    """

    def __init__(self, domain=None):

        Reflective_boundary.__init__(self, domain)
    
        self.conc=domain.quantities['concentration']
    
        self.evolved_quantities = np.zeros(4, np.float)


    def __repr__(self):
        return 'Reflective_boundary Sed'

    def evaluate(self, vol_id, edge_id):
        """Calculate reflections (reverse outward momentum).

        vol_id   
        edge_id  
        """
    
        q = self.evolved_quantities
        q[0] = self.stage[vol_id, edge_id]
        q[1] = self.xmom[vol_id, edge_id]
        q[2] = self.ymom[vol_id, edge_id]
        q[3] = self.conc[vol_id,edge_id]

        normal = self.normals[vol_id, 2*edge_id:2*edge_id+2]

        r = rotate(q, normal, direction = 1)
    
        r[1] = -r[1]
        q = rotate(r, normal, direction = -1)

        return q



    def evaluate_segment(self, domain, segment_edges):
        """Apply reflective BC on the boundary edges defined by
        segment_edges
        """
    
        Reflective_boundary.evaluate_segment(self, domain, segment_edges)

        if segment_edges is None:
            return
        if domain is None:
            return


        ids = segment_edges
        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]

    
        Conc = domain.quantities['concentration']
    
        Conc.boundary_values[ids] = Conc.edge_values[vol_ids,edge_ids]
        
        
class Dirichlet_boundary_Sed(Dirichlet_boundary):
    """Dirichlet boundary returns constant values for the
    conserved quantities.
    
    Sed transport specific version to allow for concentration
    to be input as C instead of Ch
    
    Should NOT be used for outlets! Would not allow sediment to exit the domain
    
    """


    def __init__(self, dirichlet_values=None):
        Dirichlet_boundary.__init__(self, dirichlet_values=dirichlet_values)
        
        if dirichlet_values is None:
            msg = 'Must specify one value for each quantity'
            raise Exception(msg)

        self.dirichlet_values=np.array(dirichlet_values, np.float)



    def __repr__(self):
        return 'Dirichlet boundary Sed (%s)' %self.dirichlet_values


    def evaluate(self, vol_id=None, edge_id=None):
        return self.dirichlet_values

    def evaluate_segment(self, domain, segment_edges):

        if segment_edges is None:
            return
        if domain is None:
            return


        ids = segment_edges


        vol_ids  = domain.boundary_cells[ids]
        edge_ids = domain.boundary_edges[ids]

        q_bdry = self.dirichlet_values
        
        """
        Check if wrongly using Dirichlet_boundary_Sed for an outlet.
        If boundary elev <= min of vertex elevations, raise exception
        """
        elev_v = domain.quantities['elevation'].vertex_values
        
        msg = 'Do not use Dirichlet boundary Sed for outlet! '
        msg += 'Sediment can get trapped in the domain. Use a regular '
        msg += 'Dirichlet boundary instead.'
        assert q_bdry[0] > elev_v.min(), msg
        
        

        conserved_quantities = True
        if len(q_bdry) == len(domain.evolved_quantities):
            # enough dirichlet values to set evolved quantities
            conserved_quantities = False

        #--------------------------------------------------
        # First populate all the boundary values with
        # interior edge values
        #--------------------------------------------------
        if conserved_quantities:
            for j, name in enumerate(domain.evolved_quantities):
                Q = domain.quantities[name]
                Q.boundary_values[ids] = Q.edge_values[vol_ids,edge_ids]

        #--------------------------------------------------
        # Now over write with constant Dirichlet value
        #--------------------------------------------------
        if conserved_quantities:
            quantities = domain.conserved_quantities
        else:
            quantities = domain.evolved_quantities

        #-------------------------------------------------
        # Now update to Dirichlet values
        #-------------------------------------------------
        for j, name in enumerate(quantities):
            Q = domain.quantities[name]
            
            if name == 'concentration':
            
                if q_bdry[j] > 1:
                    msg = 'Concentration at the Dirichlet boundary must be '
                    msg += 'a fraction between 0 and 1.'
                    raise Exception(msg)
                    
                # convert from concentration to sediment volume
                elev = domain.quantities['elevation'].edge_values
                elev = elev[vol_ids,edge_ids]
                stage = domain.quantities['stage'].edge_values
                stage = stage[vol_ids,edge_ids]
                depth = stage - elev
                
                Q.boundary_values[ids] = q_bdry[j] * depth
                

                
            else:
                Q.boundary_values[ids] = q_bdry[j]   
                
                   
            