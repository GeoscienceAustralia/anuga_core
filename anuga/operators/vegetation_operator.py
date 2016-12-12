"""
Vegetation operators

"""

import numpy as num

from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

from math import sqrt, log
from anuga.config import epsilon, g

import anuga.utilities.log as log
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float


from anuga.geospatial_data.geospatial_data import Geospatial_data
import os
from scipy.interpolate import NearestNDInterpolator

#===============================================================================
# Specific Erosion operator trying to implement bed shear
#===============================================================================
class Vegetation_operator(Operator, object):
    """
    Vegetation operator that applies a drag on the flow due to the presence of veg
    """

    def __init__(self, domain,
                 indices=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                
        Operator.__init__(self, domain, description, label, logging, verbose)
        
        self.Cd = None
        self.veg = None
        
        try:
            # the value in quantity 'veg_diameter' should be the stem diameter in meters
            self.veg_diameter = self.domain.quantities['veg_diameter'].centroid_values
#             self.veg_diameter = self.veg_diameter * 0.00104725
            self.veg_diameter[self.veg_diameter < 0] = 0
            self.domain.quantities['veg_diameter'].\
            	set_values(self.veg_diameter, location = 'centroids')
        except:
            self.veg_diameter = None
            
        try:
            # the value in quantity 'veg_spacing' should be the stem spacing in meters
            self.veg_spacing = self.domain.quantities['veg_spacing'].centroid_values
        except:
            self.veg_spacing = None
            
            
        try:
            diff = self.domain.get_quantity('diffusivity')
        except:
            Quantity(domain, name='diffusivity', register=True)
            
        self.domain.set_use_kinematic_viscosity(True)
            
            
        self.xmom = self.domain.quantities['xmomentum'].centroid_values
        self.ymom = self.domain.quantities['ymomentum'].centroid_values
        self.depth = self.domain.quantities['height'].centroid_values


    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        
        
        """
        
        self.dt = self.get_timestep()
        
        if self.veg_diameter is None:
            self.veg_diameter = self.domain.quantities['veg_diameter'].centroid_values
            self.veg_diameter[self.veg_diameter < 0] = 0
            self.veg_diameter *= 1
            
            self.domain.quantities['veg_diameter'].\
            	set_values(self.veg_diameter, location = 'centroids')
            
        if self.veg_spacing is None:
            self.veg_spacing = self.domain.quantities['veg_spacing'].centroid_values
            self.veg_spacing[self.veg_spacing < 0] = 0
            
         
            self.domain.quantities['veg_spacing'].\
            	set_values(self.veg_spacing, location = 'centroids')
            
        if self.veg is None:
        
            self.veg = num.zeros(self.depth.shape)
        
            self.veg[self.veg_spacing > 0] = self.veg_diameter[self.veg_spacing > 0] / self.veg_spacing[self.veg_spacing > 0]**2
            self.ad = self.veg * self.veg_diameter
            
        if self.Cd is None:
            self.calculate_drag_coefficient()
            
        
        self.wet_cells = (self.depth > 0.2) & (self.ad > 0)
        
        if sum(self.wet_cells)>0:
        
            self.depth_w = self.depth[self.wet_cells]
        
            self.velocity, xvel, yvel = self.calculate_velocity()
#             self.calculate_diffusivity()
            
            
        
            Fd_x = self.Cd_veg[self.wet_cells] * xvel**2
            Fd_y = self.Cd_veg[self.wet_cells] * yvel**2
            
            
            
            dxv = num.sign(xvel) * num.abs(Fd_x) * self.dt
#             dxv[dxv > xvel] =  xvel[dxv > xvel]
            
            dyv = num.sign(yvel) * num.abs(Fd_y) * self.dt
#             dyv[dyv > yvel] =  yvel[dyv > yvel]

#             print self.ymom[self.wet_cells]
#             print yvel
#             print Fd_y
#             print dyv
#             print self.depth_w
#             print 10 * '-'

#             dxv[num.sign(dxv) != num.sign(xvel)] = -1. * xvel[num.sign(dxv) != num.sign(xvel)]
#             dyv[num.sign(dyv) != num.sign(yvel)] = -1. * yvel[num.sign(dyv) != num.sign(yvel)]
#         
            xvel_v = (xvel - dxv) * self.depth_w
            yvel_v = (yvel - dyv) * self.depth_w

            self.xmom[self.wet_cells] = xvel_v
            self.ymom[self.wet_cells] = yvel_v
            
#             print self.xmom.max(), self.xmom.min()
#             print '-' * 10

    #         xvel_v[dry_cells] = self.xmom[dry_cells]
    #         yvel_v[dry_cells] = self.ymom[dry_cells]
    
#             print self.ymom[~self.wet_cells].min(), self.ymom[self.wet_cells].min()
        
            self.domain.quantities['xmomentum'].\
                set_values(self.xmom, location = 'centroids')
            
            self.domain.quantities['ymomentum'].\
                set_values(self.ymom, location = 'centroids')
                

    def calculate_velocity(self):
        """
        Calcultes the magnitude of the flow velocity using flux limiting, as described
        in section 16.5 of the ANUGA docs
        """
        
        ho = 1e-6
        
        xvel = (self.xmom[self.wet_cells] * self.depth_w) / (self.depth_w**2 + ho)
        yvel = (self.ymom[self.wet_cells] * self.depth_w) / (self.depth_w**2 + ho)
        
        U = num.sqrt(xvel**2. + yvel**2)
        
        return U, xvel, yvel


    
    def calculate_drag_coefficient(self):
        '''
        Calculate the drag coefficient Cd as a function of ad using
        the curve fitted to Figure 6 in Nepf (1999)
        '''
        
        self.Cd = (56.11 * self.ad**2
                   - 15.28 * self.ad
                   + 1.3
                   - 0.0005465 / self.ad)
        self.Cd[self.ad < 0.006] = 1.2
        
        self.Cd_veg = 0.5 * self.Cd * self.veg
        

        
    def calculate_diffusivity(self):
        
#         self.momentum = num.sqrt(self.xmom**2 + self.ymom**2)
#         self.velocity = self.momentum / (self.depth + epsilon)
    
        # mixing length
        mix_length = num.zeros(self.veg_diameter.shape)
        diffusivity = num.zeros(self.veg_diameter.shape)
        
        # for all cells, the transition from mix_length = depth to a linear
        # form happens at veg_spacing = depth, which is ad = diameter**2 / depth**2
        
        ad_deltaS = self.veg_diameter[self.wet_cells]**2 / (self.depth_w)**2
        
        mix_length_slope = (self.veg_diameter[self.wet_cells] - self.depth_w) / (0.01 - ad_deltaS)
        
        mix_length = (mix_length_slope * self.ad[self.wet_cells] +
                    (self.depth_w - mix_length_slope * ad_deltaS))
        
        mix_length[self.veg_spacing[self.wet_cells] > self.depth_w] = \
                   self.depth_w[self.veg_spacing[self.wet_cells] > self.depth_w]
                   
        mix_length[self.ad[self.wet_cells] > 0.01] = self.veg_diameter[self.wet_cells][self.ad[self.wet_cells] > 0.01]
                        
        # turbulent kinetic energy
        Cb = 0.001
        k = ((1 - self.ad[self.wet_cells]) * Cb + (self.Cd[self.wet_cells] * self.ad[self.wet_cells])**0.66) * self.velocity**2
        
        # total diffusivity
        diffusivity[self.wet_cells] = (num.sqrt(k)**0.5 * mix_length +
                    self.ad[self.wet_cells] * self.velocity * self.veg_diameter[self.wet_cells])
                    
#         diffusivity[~self.wet_cells] = 0
        
        self.domain.quantities['diffusivity'].\
                set_values(diffusivity, location = 'centroids')



    def set_veg_quantity(self, name_in, quantity_name=None):
    
        print 'Converting asc files'
        
        Quantity(self.domain, name=quantity_name, register=True)
    
        self.generic_asc2dem(name_in + '.asc', quantity_name=quantity_name)
        self.generic_dem2pts(name_in + '.dem', quantity_name=quantity_name)
        
        
        self.set_quantity_NNeigh(quantity_name, name_in + '.pts')


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
        self.set_values_NNeigh(name, filename)


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
        coord = self.domain.get_nodes(absolute=True)
    
        # extract the data from the pts file 
        G_data = Geospatial_data(filename)
        points = G_data.get_data_points(absolute=True)
        z = G_data.get_attributes(attribute_name=None)
    
        # create interpolator
        interp = NearestNDInterpolator( points, z )
    
        # set quantity at centroids
        z = interp( coord )
        z = z[num.newaxis, :].transpose()
    
    #     z_c = np.concatenate((coord, z), axis=1 )

    
    
        self.domain.quantities[name].set_values_from_array(z, location = 'unique vertices')

                            

    def generic_dem2pts(self, name_in, name_out=None, quantity_name=None,
                easting_min=None, easting_max=None,
                northing_min=None, northing_max=None,
                use_cache=False, verbose=False,):
        """Read raster file from the following NetCDF format (.dem)
        Generic function, created from dem2pts

        Example:

        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........

        name_in may be a .asc or .dem file to be converted.

        Convert to NetCDF pts format which is

        points:  (Nx2) float array
        elevation: N float array
        """

        kwargs = {'name_out': name_out,
                  'quantity_name': quantity_name,
                  'easting_min': easting_min,
                  'easting_max': easting_max,
                  'northing_min': northing_min,
                  'northing_max': northing_max,
                  'verbose': verbose}

        if use_cache is True:
            from caching import cache
            result = cache(self._generic_dem2pts, name_in, kwargs,
                           dependencies = [name_in],
                           verbose = verbose)

        else:
            result = apply(self._generic_dem2pts, [name_in], kwargs)

        return result


    def _generic_dem2pts(self, name_in, name_out=None, quantity_name=None, verbose=False,
                easting_min=None, easting_max=None,
                northing_min=None, northing_max=None):
        """Read raster from the following NetCDF format (.dem)

        Internal function. See public function generic_dem2pts for details.
        """

        # FIXME: Can this be written feasibly using write_pts?

        import os
        from anuga.file.netcdf import NetCDFFile

        root = name_in[:-4]

        if name_in[-4:] == '.asc':
            intermediate = root + '.dem'
            if verbose:
                log.critical('Preconvert %s from asc to %s' % \
                                        (name_in, intermediate))
            asc2dem(name_in)
            name_in = intermediate
        elif name_in[-4:] != '.dem':
            raise IOError('Input file %s should be of type .asc or .dem.' % name_in)

        if name_out != None and basename_out[-4:] != '.pts':
            raise IOError('Input file %s should be of type .pts.' % name_out)

        # Get NetCDF
        infile = NetCDFFile(name_in, netcdf_mode_r) 

        if verbose: log.critical('Reading raster from %s' % (name_in))

        ncols = int(infile.ncols)
        nrows = int(infile.nrows)
        xllcorner = float(infile.xllcorner)  # Easting of lower left corner
        yllcorner = float(infile.yllcorner)  # Northing of lower left corner
        cellsize = float(infile.cellsize)
        NODATA_value = float(infile.NODATA_value)

        dem_elevation = infile.variables[quantity_name]

        zone = int(infile.zone)
        false_easting = float(infile.false_easting)
        false_northing = float(infile.false_northing)

        #print ncols, nrows, xllcorner,yllcorner, cellsize, NODATA_value, zone


        # Text strings
        projection = infile.projection
        datum = infile.datum
        units = infile.units

        #print projection, datum, units

        # Get output file
        if name_out == None:
            ptsname = root + '.pts'
        else:
            ptsname = name_out

        if verbose: log.critical('Store to NetCDF file %s' % ptsname)

        # NetCDF file definition
        outfile = NetCDFFile(ptsname, netcdf_mode_w)

        # Create new file
        outfile.institution = 'Geoscience Australia'
        outfile.description = 'NetCDF pts format for compact and portable ' \
                              'storage of spatial point data'

        # Assign default values
        if easting_min is None: easting_min = xllcorner
        if easting_max is None: easting_max = xllcorner + ncols*cellsize
        if northing_min is None: northing_min = yllcorner
        if northing_max is None: northing_max = yllcorner + nrows*cellsize


        #print easting_min, easting_max, northing_min, northing_max

        # Compute offsets to update georeferencing
        easting_offset = xllcorner - easting_min
        northing_offset = yllcorner - northing_min

        # Georeferencing
        outfile.zone = zone
        outfile.xllcorner = easting_min # Easting of lower left corner
        outfile.yllcorner = northing_min # Northing of lower left corner
        outfile.false_easting = false_easting
        outfile.false_northing = false_northing

        outfile.projection = projection
        outfile.datum = datum
        outfile.units = units

        # Grid info (FIXME: probably not going to be used, but heck)
        outfile.ncols = ncols
        outfile.nrows = nrows

        dem_elevation_r = num.reshape(dem_elevation, (nrows, ncols))
        totalnopoints = nrows*ncols



        #========================================
        # Do the preceeding with numpy
        #========================================
        y = num.arange(nrows,dtype=num.float)
        y = yllcorner + (nrows-1)*cellsize - y*cellsize

        x = num.arange(ncols,dtype=num.float)
        x = xllcorner + x*cellsize

        xx,yy = num.meshgrid(x,y)

        xx = xx.flatten()
        yy = yy.flatten()

    
        flag = num.logical_and(num.logical_and((xx <= easting_max),(xx >= easting_min)),
                               num.logical_and((yy <= northing_max),(yy >= northing_min)))

    
        dem = dem_elevation[:].flatten()


        id = num.where(flag)[0]

        xx = xx[id]
        yy = yy[id]
        dem = dem[id]


        clippednopoints = len(dem)
        #print clippedpoints
    
        #print xx
        #print yy
        #print dem

        data_flag = dem != NODATA_value

        data_id = num.where(data_flag)

        xx = xx[data_id]
        yy = yy[data_id]
        dem = dem[data_id]

        nn = clippednopoints - len(dem)

        nopoints = len(dem)


        if verbose:
            log.critical('There are %d values in the raster' % totalnopoints)
            log.critical('There are %d values in the clipped raster'
                         % clippednopoints)
            log.critical('There are %d NODATA_values in the clipped raster' % nn)

        outfile.createDimension('number_of_points', nopoints)
        outfile.createDimension('number_of_dimensions', 2) #This is 2d data

        # Variable definitions
        outfile.createVariable('points', netcdf_float, ('number_of_points',
                                                        'number_of_dimensions'))
        outfile.createVariable(quantity_name, netcdf_float, ('number_of_points',))

        # Get handles to the variables
        points = outfile.variables['points']
        elevation = outfile.variables[quantity_name]

        points[:,0] = xx - easting_min
        points[:,1] = yy - northing_min
        elevation[:] = dem


        infile.close()
        outfile.close()

    def generic_asc2dem(self, name_in, name_out=None, quantity_name=None,
                                      use_cache=False,
                                      verbose=False):
        """Read raster from the following ASCII format (.asc)
        Generic function, created from asc2dem

        Example:
        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........

        Convert name_in (.asc) to NetCDF format (.dem)
        mimicking the ASCII format closely.

        An accompanying file with same basename but extension .prj must exist
        and is used to fix the UTM zone, datum, false northings and eastings.

        The prj format is assumed to be as

        Projection    UTM
        Zone          56
        Datum         WGS84
        Zunits        NO
        Units         METERS
        Spheroid      WGS84
        Xshift        0.0000000000
        Yshift        10000000.0000000000
        Parameters
        """

        kwargs = {'name_out': name_out, 'quantity_name': quantity_name, 'verbose': verbose}

        if use_cache is True:
            from caching import cache
            result = cache(self._generic_convert_dem_from_ascii2netcdf, name_in, kwargs,
                           dependencies=[name_in,
                                         name_in[:-4] + '.prj'],
                           verbose=verbose)

        else:
            result = apply(self._generic_convert_dem_from_ascii2netcdf, [name_in], kwargs)

        return result


    def _generic_convert_dem_from_ascii2netcdf(self, name_in, name_out = None,
                                       quantity_name = None, verbose = False):
        """Read raster from the following ASCII format (.asc)

        Internal function. See public function convert_dem_from_ascii2netcdf
        for details.
        """

        import os
        from anuga.file.netcdf import NetCDFFile

        root = name_in[:-4]

        # Read Meta data
        if verbose: log.critical('Reading METADATA from %s' % (root + '.prj'))

        metadatafile = open(root + '.prj')
        metalines = metadatafile.readlines()
        metadatafile.close()

        L = metalines[0].strip().split()
        assert L[0].strip().lower() == 'projection'
        projection = L[1].strip()                   #TEXT

        L = metalines[1].strip().split()
        assert L[0].strip().lower() == 'zone'
        zone = int(L[1].strip())

        L = metalines[2].strip().split()
        assert L[0].strip().lower() == 'datum'
        datum = L[1].strip()                        #TEXT

        L = metalines[3].strip().split()
        assert L[0].strip().lower() == 'zunits'     #IGNORE
        zunits = L[1].strip()                       #TEXT

        L = metalines[4].strip().split()
        assert L[0].strip().lower() == 'units'
        units = L[1].strip()                        #TEXT

        L = metalines[5].strip().split()
        assert L[0].strip().lower() == 'spheroid'   #IGNORE
        spheroid = L[1].strip()                     #TEXT

        L = metalines[6].strip().split()
        assert L[0].strip().lower() == 'xshift'
        false_easting = float(L[1].strip())

        L = metalines[7].strip().split()
        assert L[0].strip().lower() == 'yshift'
        false_northing = float(L[1].strip())

        if name_in[-4:] != '.asc':
            raise IOError('Input file %s should be of type .asc.' % name_in)

        #Read DEM data
        datafile = open(name_in)

        if verbose: log.critical('Reading raster from %s' % (name_in))

        lines = datafile.readlines()
        datafile.close()

        if verbose: log.critical('Got %d lines' % len(lines))

        ncols = int(lines[0].split()[1].strip())
        nrows = int(lines[1].split()[1].strip())

        # Do cellsize (line 4) before line 2 and 3
        cellsize = float(lines[4].split()[1].strip())

        # Checks suggested by Joaquim Luis
        # Our internal representation of xllcorner
        # and yllcorner is non-standard.
        xref = lines[2].split()
        if xref[0].strip() == 'xllcorner':
            xllcorner = float(xref[1].strip()) # + 0.5*cellsize # Correct offset
        elif xref[0].strip() == 'xllcenter':
            xllcorner = float(xref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % xref[0].strip()
            raise Exception, msg

        yref = lines[3].split()
        if yref[0].strip() == 'yllcorner':
            yllcorner = float(yref[1].strip()) # + 0.5*cellsize # Correct offset
        elif yref[0].strip() == 'yllcenter':
            yllcorner = float(yref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % yref[0].strip()
            raise Exception, msg

        NODATA_value = int(float(lines[5].split()[1].strip()))

        assert len(lines) == nrows + 6

        if name_out == None:
            netcdfname = name_in[:-4]+'.dem'
        else:
            netcdfname = name_out + '.dem'

        if verbose: log.critical('Store to NetCDF file %s' % netcdfname)

        # NetCDF file definition
        fid = NetCDFFile(netcdfname, netcdf_mode_w)

        #Create new file
        fid.institution = 'Geoscience Australia'
        fid.description = 'NetCDF DEM format for compact and portable storage ' \
                          'of spatial point data'

        fid.ncols = ncols
        fid.nrows = nrows
        fid.xllcorner = xllcorner
        fid.yllcorner = yllcorner
        fid.cellsize = cellsize
        fid.NODATA_value = NODATA_value

        fid.zone = zone
        fid.false_easting = false_easting
        fid.false_northing = false_northing
        fid.projection = projection
        fid.datum = datum
        fid.units = units

        # dimension definitions
        fid.createDimension('number_of_rows', nrows)
        fid.createDimension('number_of_columns', ncols)

        # variable definitions
        fid.createVariable(quantity_name, netcdf_float, ('number_of_rows',
                                                       'number_of_columns'))

        # Get handles to the variables
        elevation = fid.variables[quantity_name]

        #Store data
        import numpy

        datafile = open(name_in)
        elevation[:,:] = numpy.loadtxt(datafile, skiprows=6)
        datafile.close()

    #    n = len(lines[6:])
    #    for i, line in enumerate(lines[6:]):
    #        fields = line.split()
    #        if verbose and i % ((n+10)/10) == 0:
    #            log.critical('Processing row %d of %d' % (i, nrows))
    #
    #        if len(fields) != ncols:
    #            msg = 'Wrong number of columns in file "%s" line %d\n' % (name_in, i)
    #            msg += 'I got %d elements, but there should have been %d\n' % (len(fields), ncols)
    #            raise Exception, msg
    #
    #        elevation[i, :] = num.array([float(x) for x in fields])

        fid.close()

        



    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return False
        

    def statistics(self):

        message = self.label + ': Veg_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message