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

import os
from scipy.interpolate import NearestNDInterpolator

#===============================================================================
# Vegetation operator applying drag to the flow
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
        
                     
        try:
            diff = self.domain.get_quantity('diffusivity')
        except:
            Quantity(domain, name='diffusivity', register=True)
            
            
        self.domain.set_use_kinematic_viscosity(True)
            
        self.xmom = self.domain.quantities['xmomentum'].centroid_values
        self.ymom = self.domain.quantities['ymomentum'].centroid_values
        self.depth = self.domain.quantities['height'].centroid_values
        
        self.num_cells = len(self.depth)
        self.mix_length = num.zeros((self.num_cells,))
        self.diffusivity = num.zeros((self.num_cells,))

        self.quantity_flag = False


    def __call__(self):
        """
        Apply vegetation drag according to veg_diameter and veg_spacing quantities
        """
        
        if not self.quantity_flag:
            self.check_quantities()    
        
        if self.quantity_flag:
        
            self.dt = self.get_timestep()
        
            self.ind = (self.depth > 0.2) & (self.ad > 0)
            self.update_quantities()
        
        
        
    def update_quantities(self):
        """
        Calculate the drag that vegetation imparts on the flow
        and update momentum quantities
        """
    
        if sum(self.ind)>0:
        
            self.mix_length[:] = 0
            self.diffusivity[:] = 0
        
            self.depth_w = self.depth[self.ind]
            self.veg_d_w = self.veg_diameter[self.ind]
            self.veg_s_w = self.veg_spacing[self.ind]
            self.ad_w = self.ad[self.ind]
        
            self.velocity, xvel, yvel = self.calculate_velocity()
            self.calculate_diffusivity()
            
            Fd_x = self.Cd_veg[self.ind] * xvel**2
            Fd_y = self.Cd_veg[self.ind] * yvel**2
            
            dxv = num.sign(xvel) * num.abs(Fd_x) * self.dt
            dyv = num.sign(yvel) * num.abs(Fd_y) * self.dt

            xvel_v = (xvel - dxv) * self.depth_w
            yvel_v = (yvel - dyv) * self.depth_w

            self.xmom[self.ind] = xvel_v
            self.ymom[self.ind] = yvel_v
            
            self.domain.set_quantity('xmomentum', self.xmom, location='centroids')
            self.domain.set_quantity('ymomentum', self.ymom, location='centroids')
            
                

    def calculate_velocity(self):
        """
        Calcultes the magnitude of the flow velocity using flux limiting, as described
        in section 16.5 of the ANUGA docs
        """
        
        ho = 1e-6
        
        xvel = (self.xmom[self.ind] * self.depth_w) / (self.depth_w**2 + ho)
        yvel = (self.ymom[self.ind] * self.depth_w) / (self.depth_w**2 + ho)
        
        U = num.sqrt(xvel**2. + yvel**2)
        
        return U, xvel, yvel


    
    def calculate_drag_coefficient(self):
        """
        Calculate the drag coefficient Cd as a function of ad using
        the curve fitted to Figure 6 in Nepf (1999)
        """
        
        self.Cd = (56.11 * self.ad**2
                   - 15.28 * self.ad
                   + 1.3
                   - 0.0005465 / self.ad)
        self.Cd[self.ad < 0.006] = 1.2
        
        self.Cd_veg = 0.5 * self.Cd * self.veg
        

        
    def calculate_diffusivity(self):
        """
        Calculates a value for the quantity diffusivity for use with
        kinematic viscosity operator.
        
        Based on the formulation of Nepf (1999)
        
        For all cells, the transition from mix_length = depth to a linear
        approximation of the mixing length happens at veg_spacing = depth,
        which is ad = diameter**2 / depth**2
        """
        
        ad_deltaS = self.veg_d_w**2 / self.depth_w**2
        
        mix_length_slope = (self.veg_d_w - self.depth_w) / (0.01 - ad_deltaS)
        
        self.mix_length = (mix_length_slope * self.ad_w +
                          (self.depth_w - mix_length_slope * ad_deltaS))
        
        ind = self.veg_s_w > self.depth_w
        self.mix_length[ind] = self.depth_w[ind]
        
        ind = self.ad_w > 0.01        
        self.mix_length[ind] = self.veg_d_w[ind]

        # turbulent kinetic energy
        Cb = 0.001
        k = ((1 - self.ad_w) * Cb + (self.Cd[self.ind] * self.ad_w)**0.66) * self.velocity**2
        
        
        # total diffusivity
        self.diffusivity[self.ind] = (num.sqrt(k)**0.5 * self.mix_length +
                                      self.ad_w * self.velocity * self.veg_d_w)
                    
                    
        self.domain.quantities['diffusivity'].\
                set_values(self.diffusivity, location = 'centroids')



    def check_quantities(self):
        """
        Check that the vegetation-related quantities have been created using
        sed_veg_quantity and calculate derived quantities.
            
        Don't calculate drag if flag is not set to True
        """
    
        try:
            self.veg_diameter = self.domain.quantities['veg_diameter'].centroid_values
            self.veg_spacing = self.domain.quantities['veg_spacing'].centroid_values
            
            self.veg = num.zeros(self.depth.shape)
            self.veg[self.veg_spacing > 0] = (self.veg_diameter[self.veg_spacing > 0] /
                                              self.veg_spacing[self.veg_spacing > 0]**2)
                                              
            self.ad = self.veg * self.veg_diameter
            self.calculate_drag_coefficient()
            
            self.quantity_flag = True
            
        except:
        
            print 'Vegetation quantities not yet defined. Continuing without veg'
            
            

    def set_veg_quantity(self, name_in, quantity_name=None, convert_file=False, load_interp=False):
        """
        Read raster file and sets a vegetation quantity.
    
        Veg quantities are: veg_diameter
                            veg_spacing
                            
        The values in the rasters should be in meters.                  
        Assumes all NODATA values in raster are zero vegetation
        
        Saves a (Nx3) .npy file of x,y,val of non-zero points
        """
                                
        Quantity(self.domain, name=quantity_name, register=True)
        points = None
    
        if load_interp:
            print 'Veg Load: loading', name_in + '_interp.npy' 
            z_ = num.load(name_in + '_interp.npy')
    
        else:
        
            if convert_file:
                print 'Veg Load: converting asc file'
                self.generic_asc2dem(name_in + '.asc', quantity_name=quantity_name)
                points = self.generic_dem2npy(name_in + '.dem', quantity_name=quantity_name)
        
            if points is None:
                print 'Veg Load: loading', name_in + '.npy' 
                points = num.load(name_in + '.npy')
            
            
            interp = NearestNDInterpolator(points[:,0:2], points[:,2])
        
            coord = self.domain.get_centroid_coordinates(absolute=True)
            z_ = interp( coord )
            
            print 'Veg Load: saving interpolated file: ', name_in + '_interp.npy'
            num.save(name_in + '_interp.npy', z_)
        
        print 'Veg Load: setting quantity', quantity_name
        self.domain.quantities[quantity_name].set_values(z_, location = 'centroids') 



    def generic_dem2npy(self, name_in, name_out=None, quantity_name=None,
                easting_min=None, easting_max=None,
                northing_min=None, northing_max=None,
                use_cache=False, verbose=False,):
        """
        Read raster file from the following NetCDF format (.dem)
        Generic function, created from dem2npy

        Example:

        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........

        name_in may be a .asc or .dem file to be converted.

        Convert to numpy array which is

        points:  (Nx3) float array
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
            result = cache(self._generic_dem2npy, name_in, kwargs,
                           dependencies = [name_in],
                           verbose = verbose)

        else:
            result = apply(self._generic_dem2npy, [name_in], kwargs)

        return result



    def _generic_dem2npy(self, name_in, name_out=None, quantity_name=None, verbose=False,
                easting_min=None, easting_max=None,
                northing_min=None, northing_max=None):
        """
        Read raster from the following NetCDF format (.dem)

        Internal function. See public function generic_dem2npy for details.
        """

        import os
        from anuga.file.netcdf import NetCDFFile

        root = name_in[:-4]

        if name_in[-4:] == '.asc':
            intermediate = root + '.dem'
            if verbose:
                log.critical('Preconvert %s from asc to %s' % \
                                        (name_in, intermediate))
            self.generic_asc2dem(name_in)
            name_in = intermediate
        elif name_in[-4:] != '.dem':
            raise IOError('Input file %s should be of type .asc or .dem.' % name_in)
            

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

        # Assign default values
        if easting_min is None: easting_min = xllcorner
        if easting_max is None: easting_max = xllcorner + ncols*cellsize
        if northing_min is None: northing_min = yllcorner
        if northing_max is None: northing_max = yllcorner + nrows*cellsize

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

        data_flag = dem != NODATA_value
        data_id = num.where(data_flag)
         
        points =  num.zeros((len(data_id[0]),3))
        points[:,0] =   xx[data_id] - easting_min
        points[:,1] = yy[data_id] - northing_min
        points[:,2] = dem[data_id]
        
        
        x_ = num.arange(0, num.ceil(points[:,0].max()))
        y_ = num.arange(0, num.ceil(points[:,1].max()))

        grid_x, grid_y = num.meshgrid(x_,y_)
        grid_z = num.zeros_like(grid_x)

        rloc = num.round(points[:,1]).astype(int)
        cloc = num.round(points[:,0]).astype(int)
        grid_z[rloc, cloc] = points[:,2]

        points = num.vstack((grid_x.flatten() + easting_min,
                             grid_y.flatten() + northing_min,
                             grid_z.flatten())).T
        
        print 'Veg Load: saving', name_in[:-4] + '.npy' 
        num.save(name_in[:-4], points)
        infile.close()
        
        return points
        
        
        
        
        
    def generic_asc2dem(self, name_in, name_out=None, quantity_name=None,
                                      use_cache=False,
                                      verbose=False):
        """
        Read raster from the following ASCII format (.asc)
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
        """
        Read raster from the following ASCII format (.asc)

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

        fid.close()



    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True
        

    def statistics(self):

        message = self.label + ': Veg_operator'
        message = message + ' on triangles '+ str(self.indices)
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message