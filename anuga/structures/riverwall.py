
import os
from anuga import barrier, numprocs, myid
import numpy

class RiverWall(object):
    """Define the elevation of 'riverwalls'. 

    These are located along each cell edge, and can have an elevation different
    from the bed elevation.

    For the DE algorithms, they are used in computing the 'z_half' value [if they
    are greater than either edge elevation]. 

    In addition, the DE fluxes at riverwalls are adjusted to agree with a weir relation,
    so long as the riverwall is not too deeply submerged.

    As of 22/04/2014, they are only implemented for DE algorithms [the shallow water
    component would be more difficult to implement with other algorithms]
        
    How the flux over the riverwall is computed:
    # We have a riverwall, which is treated as a weir.
    #
    # Say the headwater-head is the upstream water depth above the weir elevation (min 0), and 
    #   the tailwater head is the downstream water depth above the weir elevation (min 0).
    #   By definition headwater-head > tailwater-head.
    #
    # Let s = (headwater head) / (tailwater head), and h = (tailwater head)/(weir height)
    #  where 'weir height' is the elevation of the weir crest above the minimum
    #  of the left/right bed elevations
    #
    # Denote ID as the 'ideal' weir flow, computed from the hydraulic
    # formula (including a submergence correction factor from Villemonte,
    # 1947), 
    #      Q1 = 2/3*headwater_head*sqrt(g*2/3*headwater_head)*Qfactor
    #      Q2 = 2/3*tailwater_head*sqrt(g*2/3*tailwater_head)*Qfactor
    #      ID = Q1*(1-Q2/Q1)**0.385
    #
    # Denote SW as the 'shallow-water' weir flux, computed from the approximate
    # reimann solver, where the mid-edge-elevation is the weir crest elevation.
    # This makes clear sense for DE0 and DE1. Cell centroid stage/height/bed_elevation
    # are used in the flux computation
    #
    # Then the flux over the weir is computed from:
    #   w1 = min( max(s-s1, 0.)/(s2-s1), 1.0) # Factor describing relative submergence
    #   w1' = min( max(h-h1,0.)/(h2-h1), 1.0) # Factor describing absolute submergence
    #  flux_over_weir = (w1*SW + (1-w1)*ID)*( 1-w1') + (w1')*SW 
    # where s1, s2, h1, h2 are user defined parameters
    #
    # The key idea is that if s<s1, h<h1, then the ideal weir solution is
    # used. Otherwise, we linearly blend with the SW solution,
    # and the SW solution is used completely if s>s2 or h>h2

    """
    
    def __init__(self, domain):
        """Riverwall data structure

        Allows reasonably efficient storage of riverwall elevations and hydraulic
        properties

        NOTE: In domain.edge_flux_type, riverwalls correspond to the value 1

        RiverWall variables are initialised to dummy values, because even
        if no riverwalls are used, some values have to be passed to the flux computation

        Hydraulic parameters are 
        Qfactor -- Multiplicative factor for ideal weir relation (calibration coef)
        s1 --  Submergence ratio at which we start blending with the shallow water solution (<s2)
        s2 -- Submergence ratio at which we entirely use the shallow water solution  (>s1)
        h1 -- Tailwater head / weir height at which we start blending with the shallow water solution (<h2)
        h2 -- Tailwater head / weir height at which we entirely use the shallow water solution (>h1)

        # Default riverwall hydraulic parameters 
        default_riverwallPar={'Qfactor':1.0,  
                              's1': 0.9,      
                              's2': 0.95,     
                              'h1': 1.0,      
                              'h2': 1.5       
                              }

        Other variables are:

            riverwall_elevation -- Variable to hold riverwall elevations. 
                                   len = number of riverwall edges in the domain. 
                                   First entry corresponds to first riverwall edge in domain.edge_coordinates,
                                   second entry to second riverwall edge in domain.edge_coordinates, etc

            hydraulic_properties_rowIndex --  Variable to hold the row-index of the hydraulic properties table
                                               len = number of riverwall edges
                                                     in the domain, ordered like riverwall_elevation

            riverwall_edges -- Holds indices of edges in domain which are riverwalls, ordered like riverwall_elevation

            names -- list with the names of the riverwalls
                     len = number of riverwalls which cover edges in the domain

            hydraulic_variable_names -- Variables to hold the names of variables in columns of the hydraulic
                                        properties table. THE ORDER IS IMPORTANT -- C code ASSUMES IT

            ncol_hydraulic_properties -- Number of columns in the hydraulic properties table 
                                        [ = len(hydraulic_variable_names) ]

            hydraulic_properties -- Array with the hydraulic parameters for each riverwall. 
                                      number of rows = number of riverwalls which cover edges in the domain
                                      number of cols = number of hydraulic variable names

                                    
            input_riverwall_geo, input_riverwall_par -- holds input information

        """
        self.domain=domain

        default_float=-9.0e+20
        default_int=-1_000_000_000
        self.riverwall_elevation=numpy.array([default_float])

        self.hydraulic_properties_rowIndex=numpy.array([default_int]).astype(int)

        self.names=[ ]

        # Default riverwall hydraulic parameters 
        self.default_riverwallPar={'Qfactor':1.0,  
                                   's1': 0.9,      
                                   's2': 0.95,     
                                   'h1': 1.0,      
                                   'h2': 1.5       
                                   }

        # DO NOT CHANGE THE ORDER OF hydraulic_variable_names
        # It needs to match hard-coded assumptions in C [compute_fluxes_central]
        # If you add a variable, append it to the end of hydraulic_variable_names
        self.hydraulic_variable_names=('Qfactor', 's1', 's2', 'h1', 'h2')

        self.ncol_hydraulic_properties=len(self.hydraulic_variable_names)
        # Variable to hold the riverwall hydraulic properties in a table
        #  number of rows = number of riverwalls which cover edges in the domain
        #  number of cols = number of hydraulic variable names
        self.hydraulic_properties=numpy.array([ [default_float] ])

        # Variable to hold the indices of riverwall edges
        #    len = number of riverwall edges in the domain
        self.riverwall_edges=numpy.array([default_int])

        # Input info
        self.input_riverwall_geo=None
        self.input_riverwallPar=None


    def create_riverwalls(self, riverwalls, riverwallPar={ },
                          default_riverwallPar={ },
                          tol=1.0e-4, verbose=True, 
                          output_dir=None):
        """Add riverwalls at chosen locations along the mesh

        As of 22/04/2014, these only work with DE algorithms [for which the concept is natural]

        The walls MUST EXACTLY COINCIDE with edges along the mesh 
        
        You can force the mesh to do this by passing riverwalls.values() 
        to the 'breaklines' argument in the function create_mesh_from_regions. You
        may also need to set the maximum_triangle_area for regions, if the breaklines
        split the region.  Do this with the regionPtArea argument in
        create_mesh_from_regions.

        As of 22/04/2014, the computational method used here is very 'brute-force'
        For each segment on each riverwall, every point in the mesh is checked to see
        if it is on the segment. A possibly faster/less memory method would be to
        'walk' through connected edges on the mesh.

        Inputs:
            riverwalls: Dictionary of '3D polygons', containing xyz locations of named riverwalls.

                exampleRiverWalls = { # Left bank n1 -- 
                                      'n1': [[1.0, 1000., 2.],
                                             [1.0, 50., 3.]],
                                      # Left bank n2 
                                       'n2': [[2.0, 30., 1.0],
                                              [3.0, 20., 2.], 
                                              [2.5, 15., 1.5]]
                                    }

            riverwallPar: Dictionary containing a dictionary of named hydraulic parameters for each named riverwall
                          If parameters are not provided, default values will be used. 
                          See the help for class 'RiverWall' for an explanation of these

                exampleRiverWallPar = {'n2': {'Qfactor':0.5} }
                    This would use a Qfactor of 0.5 for riverwall 'n2', while the other riverwall would have the default values

            default_riverwallPar:  Dictionary containing default values of the riverwall parameters, to be used
                                   if they are not explicitly set.
                                   If not provided, defaults from __init__ are used. See the help for class 'RiverWall' for more info

                example_default_riverwallPar = {'Qfactor':1.5,  
                                                's1': 0.9,      
                                                's2': 0.95,     
                                                'h1': 1.0,      
                                                'h2': 1.5       
                                               }

                example_default_riverwallPar = {'Qfactor':1.5,  
                                                's1': 10000.,      
                                                's2': 20000.     
                                               } # Here h1/h2 defaults will come from __init__
                 

            tol: Edges will be assigned a riverwall elevation if they are within 'tol' of
                 a segment in riverwalls. Round-off error means this should not be too small.

            verbose: TRUE/FALSE Print lots of information

            output_dir: Text representation of riverwalls is written to output_dir, unless output_dir=None

        Outputs:
            None, but it sets domain.edge_flux_type = 1 for every edge on a riverwall
            Also, it adds a RiverWall object domain.riverwallData to the domain
            The latter contains the riverwall heights, names, and hydraulic properties for each edge, in
              a way we can pass in/out of C code.

        """

        #if(verbose and myid==0):        
        #    print ' '
        #    print 'WARNING: Riverwall is an experimental feature'
        #    print '         At each riverwall edge, we place a thin "wall" between'
        #    print '         the 2 edges -- so each sees its neighbour edge as having'
        #    print '         bed elevation = max(levee height, neighbour bed elevation)'
        #    print '         We also force the discharge with a weir relation, blended with'
        #    print '         the shallow water flux solution as the ratio (min_head)/(weir_height) becomes '
        #    print '         large, or the ratio (downstream_head)/(upstream_head) becomes large'
        #    print ' '
        #    print '  It works in parallel, but you must use domain.riverwallData.create_riverwall AFTER distributing the mesh'
        #    print ' '

        # NOTE: domain.riverwallData is initialised in shallow_water_domain.py for DE algorithms
        domain=self.domain

        
        # Check flow algorithm
        if(not domain.get_using_discontinuous_elevation()):
            raise Exception('Riverwalls are currently only supported for discontinuous elevation flow algorithms')

        if(len(self.names)>0):
            # Delete any existing riverwall data
            # The function cannot presently be used to partially edit riverwall data
            if(verbose):
                print('Warning: There seems to be existing riverwall data')
                print('It will be deleted and overwritten with this function call')
            domain.riverwallData.__init__(domain)

        # Store input parameters
        # FIXME: Is this worth it?
        self.input_riverwall_geo=riverwalls
        self.input_riverwallPar=riverwallPar

        # Update self.default_riverwallPar (defined in __init__)
        for i in list(self.default_riverwallPar.keys()):
            if(i in default_riverwallPar):
                self.default_riverwallPar[i]=default_riverwallPar[i]

        # Check that all the keys in default_riverwallPar are allowed
        for i in list(default_riverwallPar.keys()):
            if(i not in self.default_riverwallPar):
                msg='Key ', i + ' in default_riverwallPar not recognized'
                raise Exception(msg)
        # Final default river-wall parameters
        default_riverwallPar=self.default_riverwallPar
        
        # Check that all named inputs in riverwallPar correspond to names in
        # riverwall
        for i in list(riverwallPar.keys()):
            if i not in riverwalls:
                msg= 'Key ', i, ' in riverwallPar has no corresponding key in riverwall'
                raise Exception(msg)
            #
            # Check that all hydraulic parameter names in riverwallPar correspond
            # to names in default_riverwallPar
            #
            for j in list(riverwallPar[i].keys()):
                if j not in default_riverwallPar:
                    msg = 'Hydraulic parameter named ', j ,\
                          ' not recognised in default_riverwallPar'
                    raise Exception(msg)
        
        if(verbose):
            print('Setting riverwall elevations (P'+str(myid)+')...')

        # Set up geometry
        exy=domain.edge_coordinates
        llx=domain.mesh.geo_reference.get_xllcorner()
        lly=domain.mesh.geo_reference.get_yllcorner()

        # Temporary variables
        from anuga.config import max_float
        riverwall_elevation=exy[:,0]*0. - max_float
        riverwall_Qfactor=exy[:,0]*0.
        riverwall_rowIndex=exy[:,0]*0 -1.
        riverwall_rowIndex.astype(int)

        # Loop over all segments in each riverwall, and set its elevation
        nw=list(range(len(riverwalls)))
        nw_names=list(riverwalls.keys()) # Not guarenteed to be in deterministic order

        if(verbose):
            # Use variable to record messages, allows cleaner parallel printing
            printInfo='' 

        for i in nw:
            # Name of riverwall
            riverwalli_name=nw_names[i]
            # Look at the ith riverwall
            riverwalli=riverwalls[riverwalli_name]

            ns=len(riverwalli)-1

            if(verbose): 
                printInfo=printInfo + '  Wall ' + str(i) +' ....\n'

            for j in range(ns):
                if(verbose):
                    printInfo=printInfo + '    Segment ' + str(j) +' ....\n'
                # Start and end xyz coordinates
                start=riverwalli[j]
                end=riverwalli[j+1]

                if(len(start)!=3 | len(end)!=3):
                    msg='Each riverwall coordinate must have at exactly 3 values [xyz]'
                    raise Exception(msg)

                # Find length
                segLen=( (start[0]-end[0])**2+(start[1]-end[1])**2)**0.5
                if(segLen<tol):
                    if(verbose):
                        printInfo=printInfo+'  Segment with length < tolerance ' + str(tol) +' ignored\n'
                    continue 
                
                # Find edge indices which are within 'tol' of the segment
                # We use a simple, inefficient method [but likely ok in practice
                # except for very complex riverwalls]
                
                # Unit vector along segment
                se_0=-(start[0]-end[0])/segLen
                se_1=-(start[1]-end[1])/segLen

                # Vector from 'start' to every point on mesh
                # NOTE: We account for georeferencing
                pv_0 = exy[:,0]-(start[0]-llx)
                pv_1 = exy[:,1]-(start[1]-lly)

                pvLen=( pv_0**2 + pv_1**2)**0.5

                # Dot product of pv and se == along-segment distance of projection
                # of each point onto segment 
                pv_dot_se = pv_0*se_0+pv_1*se_1
                # Normal distance^2 of each point to segment
                perp_len_sq = pvLen**2.-pv_dot_se**2.
               
                # Point is on a levee if the perpendicular distance is < tol,
                # AND it is between start and end [within a tolerance]
                onLevee=(perp_len_sq<tol**2)*(pv_dot_se > 0.-tol)*(pv_dot_se<segLen+tol)
                onLevee=onLevee.nonzero()
                onLevee=onLevee[0]
                if(len(onLevee)==0):
                    continue

                if(verbose):
                    printInfo=printInfo+'       Finding ' + str(len(onLevee)) + ' edges on this segment\n'
            
                # Levee has Edge_flux_type=1  
                domain.edge_flux_type[onLevee]=1
     
                # Get edge elevations as weighted averages of start/end elevations 
                w0=pv_dot_se[onLevee]/segLen
                w0=w0*(w0>=0.0) # Enforce min of 0
                w0=w0*(w0<=1.0) + 1.0*(w0>1.0) # Max of 1
                riverwall_elevation[onLevee]= start[2]*(1.0-w0)+w0*end[2]

                # Record row index
                riverwall_rowIndex[onLevee] = i

        # Now, condense riverwall_elevation to array with length = number of riverwall edges
        #
        #  We do this to avoid storing a riverwall_elevation for every edge in the mesh
        #  However, the data structure ends up being quite complex -- maybe there is a better way?
        #
        # The zeroth index in domain.edge_flux_type which = 1 will correspond to riverwall_elevation[0]
        # The first index will correspond to riverwall_elevation[1]
        # etc
        #
        riverwallInds=(domain.edge_flux_type==1).nonzero()[0]
        # elevation
        self.riverwall_elevation=\
            riverwall_elevation[riverwallInds]
        # corresponding row in the hydraulic properties table
        self.hydraulic_properties_rowIndex=\
            riverwall_rowIndex[riverwallInds].astype(int)
        # index of edges which are riverwalls 
        self.riverwall_edges=riverwallInds

        # Record the names of the riverwalls
        self.names=nw_names

       
        # Now create the hydraulic properties table

        # Temporary variable to hold hydraulic properties table
        # This will have as many rows are there are distinct riverwalls,
        # and as many columns as there are hydraulic variables
        hydraulicTmp=numpy.zeros((len(riverwalls), len(default_riverwallPar)))*numpy.nan
       
        if(verbose):
            print(' ') 
        # Loop over every riverwall / hydraulic parameter, and set its value
        for i in nw:
            # Get the riverwall's name and specified parameters
            name_riverwalli=nw_names[i]
            if(name_riverwalli in riverwallPar):
                riverwalli_Par=riverwallPar[name_riverwalli]
            else:
                riverwalli_Par=None

            # Set the ith riverwall's hydraulic properties 
            for j, hydraulicVar in enumerate(self.hydraulic_variable_names):    
                if((riverwalli_Par is not None) and (hydraulicVar in riverwalli_Par)):
                    if(verbose): 
                        printInfo=printInfo+ '  Using provided '+ str(hydraulicVar)+' '+\
                           str(riverwalli_Par[hydraulicVar])+ ' for riverwall '+ str(name_riverwalli)+'\n'
                    hydraulicTmp[i,j]=riverwalli_Par[hydraulicVar] 
                else:
                    if(verbose): 
                        printInfo=printInfo+ '  Using default '+ str(hydraulicVar)+' '+\
                            str(default_riverwallPar[hydraulicVar])+' for riverwall '+ str(name_riverwalli)+'\n'
                    hydraulicTmp[i,j]=default_riverwallPar[hydraulicVar]

        if(verbose):
            print(' ')

        # Check that s1 < s2 
        for i in nw:
            if(hydraulicTmp[i,1]>= hydraulicTmp[i,2]):
                msg = 's1 >= s2 on riverwall ' + nw_names[i] +'. This is not allowed' 
                raise Exception(msg)
            if( (hydraulicTmp[i,1]<0.) or (hydraulicTmp[i,2] < 0.)):
                raise Exception('s1 and s2 must be positive, with s1<s2')

        # Check that h1 < h2
        for i in nw:
            if(hydraulicTmp[i,3]>= hydraulicTmp[i,4]):
                msg = 'h1 >= h2 on riverwall ' + nw_names[i] +'. This is not allowed' 
                raise Exception(msg)
            if((hydraulicTmp[i,3]<0.) or (hydraulicTmp[i,4] < 0.)):
                raise Exception('h1 and h2 must be positive, with h1<h2')
       
        # Define the hydraulic properties 
        self.hydraulic_properties=hydraulicTmp
      
        # Check for riverwall 'connectedness' errors (e.g. theoretically possible
        # to miss an edge due to round-off)
        connectedness=self.check_riverwall_connectedness(verbose=verbose)

        self.export_riverwalls_to_text(output_dir=output_dir)
        
        # Pretty printing of riverwall information in parallel
        if(verbose): 
            if domain.parallel : barrier()
            for i in range(numprocs):
                if(myid==i):
                    print('Processor '+str(myid))
                    print(printInfo)
                    print(connectedness[0])
                    msg='Riverwall discontinuity -- possible round-off error in'+\
                         'finding edges on wall -- try increasing value of tol'
                    if(not connectedness[1]):
                        raise Exception(msg)
                if domain.parallel : barrier()
        return 
    
    #####################################################################################

    def get_centroids_corresponding_to_edgeInds(self, riverwalledgeInds):
        """ 
          Get indices of centroids containing edges with indices riverwalledgeInds
        """
        riverwallCentInds=numpy.floor(riverwalledgeInds/3.)
        riverwallCentInds=riverwallCentInds.astype(int)

        return riverwallCentInds

    #####################################################################################

    def get_vertices_corresponding_to_edgeInds(self, riverwalledgeInds, checkCoords=True):
        """
         Get indices of vertices corresponding to edges at index riverwalledgeInds
    
         Since each edge has 2 vertices, use V1 and V2
        
         There is indeed a simple relationship between the vertex and edge indices
        """


        #riverwallCentInds=self.get_centroids_corresponding_to_edgeInds(riverwalledgeInds)

        rwEI_mod3=riverwalledgeInds%3

        # Figure out the difference between each vertex index and the edge index. Is either
        # -2, -1, 1, 2
        rwV1_adjustment= -2*(rwEI_mod3==2) -1*(rwEI_mod3==1) +1*(rwEI_mod3==0)
        rwV2_adjustment= -1*(rwEI_mod3==2) +1*(rwEI_mod3==1) +2*(rwEI_mod3==0)
        riverwallV1Inds=riverwalledgeInds+rwV1_adjustment
        riverwallV2Inds=riverwalledgeInds+rwV2_adjustment

        if(checkCoords):
            ####################################################
            # Check that vertices and edges really do correspond
            domain=self.domain
            # X coordinates
            assert( numpy.allclose(
                    domain.edge_coordinates[riverwalledgeInds,0], 
                    0.5*(domain.vertex_coordinates[riverwallV1Inds,0]+domain.vertex_coordinates[riverwallV2Inds,0]))
                    )
            # Y coordinates
            assert( numpy.allclose(
                    domain.edge_coordinates[riverwalledgeInds,1], 
                    0.5*(domain.vertex_coordinates[riverwallV1Inds,1]+domain.vertex_coordinates[riverwallV2Inds,1]))
                    )
            ####################################################

        return riverwallV1Inds, riverwallV2Inds
    
    #####################################################################################
    def is_vertex_on_boundary(self, vertexIndices):
        """
            Determine whether a vertex is on the boundary of the domain
            (i.e. it's connected with an edge that is a boundary edge)

            INPUTS: self -- riverwallData
                    vertexIndices -- indices of vertices on the domain which are on the riverwall

            OUTPUT:
                    TRUE if the vertex is on a domain boundary, FALSE otherwise

        """
        domain=self.domain

        # Get edge/vertex indices for boundaries
        boundary_index_info=list(domain.boundary.keys())
        boundary_edges=[ boundary_index_info[i][0]*3+boundary_index_info[i][1] for i in range(len(boundary_index_info))]
        boundary_edges=numpy.array(boundary_edges)
        tmp=self.get_vertices_corresponding_to_edgeInds(boundary_edges, checkCoords=False)
        boundary_vertices=numpy.hstack([tmp[0], tmp[1]]).tolist()

        # Get 'unique' vertex coordinates on boundary 
        node_complex=domain.vertex_coordinates[boundary_vertices,0]+1j*domain.vertex_coordinates[boundary_vertices,1]

        # Get riverwall vertex coordinates as complex numbers (for equality testing)
        complex_vertex_coords=domain.vertex_coordinates[vertexIndices.tolist(),0]+\
                                1j*domain.vertex_coordinates[vertexIndices.tolist(),1]

        # Flag telling us if the vertex is on the boundary (1=True, 0=False)
        isOnBoundary=[ 1 if complex_vertex_coords[i] in node_complex else 0 for i in range(len(complex_vertex_coords))]
        isOnBoundary=numpy.array(isOnBoundary)

        return isOnBoundary

    #####################################################################################
    def check_riverwall_connectedness(self, verbose=True):
        """
            We expect riverwalls to be connected 
             (although they can pass through the bounding polygon several times, especially in parallel)
            Round-off can potentially cause riverwalls to be dis-connected
            Use this routine to check for that

            Basically, we search for edges which are connected to vertices which
                themselves are not connected to other edges. We ignore vertices on the domain's bounding-polygon
            
            For a continuous riverwall, there can be at most 2 endpoints inside the domain

            Otherwise, the riverwall is discontinuous, which is an error 

        """

        domain = self.domain

        # Preliminary definitions
        isConnected = True 
        printInfo = '' 

        if(len(self.names)==0):
            if(verbose):
                printInfo = printInfo+'  There are no riverwalls (P'+str(myid)+')\n'
            return [printInfo, isConnected]

        # Shorthand notation
        rwd = self

        for i, name in enumerate(rwd.names):
            # Get indices of edges on this riverwall
            riverwalledgeInds = rwd.riverwall_edges[(rwd.hydraulic_properties_rowIndex==i).nonzero()[0]]

            if(len(riverwalledgeInds)==0):
                printInfo = printInfo+'Riverwall '+name+' was not found on this mesh (if this is wrong, adjust tol in create_riverwalls)\n'
                continue
            # Get their corresponding vertices
            riverwallV1Inds, riverwallV2Inds = rwd.get_vertices_corresponding_to_edgeInds(riverwalledgeInds)

            # Flag telling us if vertex points are on the boundary of the model
            # Used to help detect disconnected riverwalls (due to round-off)
            v1_on_boundary = rwd.is_vertex_on_boundary(riverwallV1Inds)
            v2_on_boundary = rwd.is_vertex_on_boundary(riverwallV2Inds)

            # With discontinuous triangles, we expect edges to occur twice
            # Let's remove duplicates to simplify the analysis
            repeat = riverwalledgeInds*0
            lre = len(riverwalledgeInds)
            # edge coordinates as a complex number, for easy equality checking
            complex_edge_coordinates = domain.edge_coordinates[riverwalledgeInds,0]+\
                                       1j*domain.edge_coordinates[riverwalledgeInds,1]
            for j in range(lre-1):
                # Ignore if already checked
                if(repeat[j]==1):
                    continue 
                # Check for a dupulicate  
                dups = (complex_edge_coordinates[(j+1):lre]==complex_edge_coordinates[j]).nonzero()[0]
                if(len(dups)>0):
                    repeat[dups+j+1] = 1
            
            unique_riverwall_edge_indices = (repeat==0).nonzero()[0]

            # Finally, get 'unqiue' edges in the riverwall
            uEdges = riverwalledgeInds[unique_riverwall_edge_indices]
            uV1 = riverwallV1Inds[unique_riverwall_edge_indices]
            uV2 = riverwallV2Inds[unique_riverwall_edge_indices]
            uV1_boundary = v1_on_boundary[unique_riverwall_edge_indices]
            uV2_boundary = v2_on_boundary[unique_riverwall_edge_indices]

            # Next, count how many times each vertex value occurs. 
            # For a 'connected' riverwall, we only want 2 edges where a vertex occurs only once,
            #   unless the vertex is on the boundary of the domain
            lure = len(uEdges)
            complex_v1_coordinates = domain.vertex_coordinates[uV1,0]+\
                                     1j*domain.vertex_coordinates[uV1,1]
            complex_v2_coordinates = domain.vertex_coordinates[uV2,0]+\
                                     1j*domain.vertex_coordinates[uV2,1]
            v1Counter = uEdges*0
            v2Counter = uEdges*0
            for j in range(lure):
                v1Counter[j] = (complex_v1_coordinates==complex_v1_coordinates[j]).sum()+\
                               (complex_v2_coordinates==complex_v1_coordinates[j]).sum()
                v2Counter[j] = (complex_v1_coordinates==complex_v2_coordinates[j]).sum()+\
                               (complex_v2_coordinates==complex_v2_coordinates[j]).sum()

            num_disconnected_edges = ((v1Counter==1)*(1-uV1_boundary)).sum()+\
                                     ((v2Counter==1)*(1-uV2_boundary)).sum()
          
            if(verbose):    
                printInfo = printInfo+ '  On riverwall '+ str(name) +' there are '+ str(num_disconnected_edges)+\
                         ' endpoints inside the domain [ignoring points on the boundary polygon] (P'+str(myid)+')\n'

            if(num_disconnected_edges <= 2):
                if(verbose):
                    pass
                    #printInfo=printInfo+ "  This is consistent with a continuous wall \n"
            else:
                isConnected = False
                printInfo = printInfo + '  Riverwall ' + name +' appears to be discontinuous. (P'+str(myid)+')\n'+\
                    '  This suggests there is a gap in the wall, which should not occur\n'
        
        return [printInfo, isConnected]

    ###################################################################

    def export_riverwalls_to_text(self, output_dir=None):
        """
            Function for dumping riverwall outputs to text file (useful for QC)

            This will overwrite any existing files with the same location/name

            INPUT: output_dir = Directory where the outputs will be written

            OUTPUT:
                    None, but writes files as a side effect

        """
        if(output_dir is None):
            return
    
        if(myid == 0):
            # Make output directory
            try: 
                os.mkdir(output_dir)
            except:
                pass
            # Make output files with empty contents
            for i, riverWallFile in enumerate(self.names):
                newFile=open(output_dir + '/' + os.path.splitext(os.path.basename(riverWallFile))[0] + '.txt','w')
                # Write hydraulic variable information
                hydraulicVars=self.hydraulic_properties[i,:]
                newFile.write('## Hydraulic Variable values below ## \n')
                newFile.write(str(self.hydraulic_variable_names) + '\n')
                newFile.write(str(hydraulicVars) + '\n')
                newFile.write('\n')
                newFile.write('## xyElevation at edges below. Order may be erratic for parallel runs ##\n')
                newFile.close()
        else:
            pass
        



        domain = self.domain

        # The other processes might try to write into file
        # before process 0 has created file, so we need a 
        # barrier
        if domain.parallel: barrier()    

        # Now dump the required info to the files
        for i in range(numprocs):
            # Write 1 processor at a time
            if(myid == i):
                for j, riverWallname in enumerate(self.names):
                    # Get xyz data for riverwall j
                    riverWallInds = (self.hydraulic_properties_rowIndex==j).nonzero()[0].tolist()
                    riverWallDomainInds = self.riverwall_edges[riverWallInds].tolist()
                    myXCoords = domain.edge_coordinates[riverWallDomainInds,0] + domain.geo_reference.xllcorner
                    myYCoords = domain.edge_coordinates[riverWallDomainInds,1] + domain.geo_reference.yllcorner
                    myElev = self.riverwall_elevation[riverWallInds]
               
                    # Open file for appending data 
                    theFile = open(output_dir + '/' + os.path.splitext(os.path.basename(riverWallname))[0] + '.txt','a')
                    for k in range(len(myElev)):
                        theFile.write(str(myXCoords[k]) + ',' + str(myYCoords[k]) + ',' + str(myElev[k]) + '\n')
                    theFile.close()
                        
            else:
                pass
 

        return
