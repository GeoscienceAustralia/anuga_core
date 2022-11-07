
import numpy as num
import anuga.utilities.log as log
from anuga.file.netcdf import NetCDFFile

import anuga

from anuga.config import max_float
from anuga.config import netcdf_float, netcdf_float32, netcdf_int
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     ensure_geo_reference
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a

class Write_sts(object):
    """ A class to write STS files.
    """
    sts_quantities = ['stage', 'xmomentum', 'ymomentum']
    RANGE = '_range'
    EXTREMA = ':extrema'

    def __init__(self):
        pass

    def store_header(self,
                     outfile,
                     times,
                     number_of_points,
                     description='Converted from URS mux2 format',
                     sts_precision=netcdf_float,
                     verbose=False):
        """Write a header to the underlying data file.

        outfile          handle to open file to write
        times            list of the time slice times *or* a start time
        number_of_points the number of URS gauge sites
        description      description string to write into the STS file
        sts_precision    format of data to write (netcdf constant ONLY)
        verbose          True if this function is to be verbose

        If 'times' is a list, the info will be made relative.
        """

        outfile.institution = 'Geoscience Australia'
        outfile.description = description

        try:
            revision_number = anuga.get_revision_number()
        except:
            revision_number = None

        # Allow None to be stored as a string
        outfile.revision_number = str(revision_number)

        # Start time in seconds since the epoch (midnight 1/1/1970)
        # This is being used to seperate one number from a list.
        # what it is actually doing is sorting lists from numeric arrays.
        if isinstance(times, (list, num.ndarray)):
            number_of_times = len(times)
            times = ensure_numeric(times)
            if number_of_times == 0:
                starttime = 0
            else:
                starttime = times[0]
                times = times - starttime  #Store relative times
        else:
            number_of_times = 0
            starttime = times

        outfile.starttime = starttime

        # Dimension definitions
        outfile.createDimension('number_of_points', number_of_points)
        outfile.createDimension('number_of_timesteps', number_of_times)
        outfile.createDimension('numbers_in_range', 2)

        # Variable definitions
        outfile.createVariable('permutation', netcdf_int, ('number_of_points',))
        outfile.createVariable('x', sts_precision, ('number_of_points',))
        outfile.createVariable('y', sts_precision, ('number_of_points',))
        outfile.createVariable('elevation', sts_precision, \
                                    ('number_of_points',))

        q = 'elevation'
        outfile.createVariable(q + Write_sts.RANGE, sts_precision,
                               ('numbers_in_range',))

        # Initialise ranges with small and large sentinels.
        # If this was in pure Python we could have used None sensibly
        outfile.variables[q + Write_sts.RANGE][0] = max_float  # Min
        outfile.variables[q + Write_sts.RANGE][1] = -max_float # Max

        self.write_dynamic_quantities(outfile, Write_sts.sts_quantities, times)

    def store_points(self,
                     outfile,
                     points_utm,
                     elevation, zone=None, new_origin=None,
                     points_georeference=None, verbose=False):

        """
        points_utm - currently a list or array of the points in UTM.
        points_georeference - the georeference of the points_utm

        How about passing new_origin and current_origin.
        If you get both, do a convertion from the old to the new.

        If you only get new_origin, the points are absolute,
        convert to relative

        if you only get the current_origin the points are relative, store
        as relative.

        if you get no georefs create a new georef based on the minimums of
        points_utm.  (Another option would be to default to absolute)

        Yes, and this is done in another part of the code.
        Probably geospatial.

        If you don't supply either geo_refs, then supply a zone. If not
        the default zone will be used.

        precondition:
             header has been called.
        """

        number_of_points = len(points_utm)
        points_utm = num.array(points_utm)

        # given the two geo_refs and the points, do the stuff
        # described in the method header
        points_georeference = ensure_geo_reference(points_georeference)
        new_origin = ensure_geo_reference(new_origin)

        if new_origin is None and points_georeference is not None:
            points = points_utm
            geo_ref = points_georeference
        else:
            if new_origin is None:
                new_origin = Geo_reference(zone, min(points_utm[:,0]),
                                                 min(points_utm[:,1]))
            points = new_origin.change_points_geo_ref(points_utm,
                                                      points_georeference)
            geo_ref = new_origin

        # At this stage I need a georef and points
        # the points are relative to the georef
        geo_ref.write_NetCDF(outfile)

        x = points[:,0]
        y = points[:,1]
        z = outfile.variables['elevation'][:]

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('More Statistics:')
            log.critical('  Extent (/lon):')
            log.critical('    x in [%f, %f], len(lat) == %d'
                         % (min(x), max(x), len(x)))
            log.critical('    y in [%f, %f], len(lon) == %d'
                         % (min(y), max(y), len(y)))
            log.critical('    z in [%f, %f], len(z) == %d'
                         % (min(elevation), max(elevation), len(elevation)))
            log.critical('geo_ref: %s' % str(geo_ref))
            log.critical('------------------------------------------------')

        #z = resize(bath_grid,outfile.variables['elevation'][:].shape)
        outfile.variables['x'][:] = points[:,0] #- geo_ref.get_xllcorner()
        outfile.variables['y'][:] = points[:,1] #- geo_ref.get_yllcorner()
        #outfile.variables['z'][:] = elevation
        outfile.variables['elevation'][:] = elevation  #FIXME HACK4

        # This updates the _range values
        q = 'elevation'
        outfile.variables[q + Write_sts.RANGE][0] = min(elevation)
        outfile.variables[q + Write_sts.RANGE][1] = max(elevation)

    def store_quantities(self, outfile, sts_precision=float,
                         slice_index=None, time=None,
                         verbose=False, **quant):
        """Write the quantity info.

        **quant is extra keyword arguments passed in. These must be
          the sts quantities, currently; stage.

        if the time array is already been built, use the slice_index
        to specify the index.

        Otherwise, use time to increase the time dimension

        Maybe make this general, but the viewer assumes these quantities,
        so maybe we don't want it general - unless the viewer is general

        precondition:
            triangulation and header have been called.
        """

        if time is not None:
            file_time = outfile.variables['time']
            slice_index = len(file_time)
            file_time[slice_index] = time

        # Write the conserved quantities from Domain.
        # Typically stage, xmomentum, ymomentum
        # other quantities will be ignored, silently.
        # Also write the ranges: stage_range
        for q in Write_sts.sts_quantities:
            if q not in quant:
                msg = 'STS file can not write quantity %s' % q
                raise Exception(msg)
            else:
                q_values = quant[q]
                outfile.variables[q][slice_index] = \
                                q_values.astype(sts_precision)

                # This updates the _range values
                q_range = outfile.variables[q + Write_sts.RANGE][:]
                q_values_min = num.min(q_values)
                if q_values_min < q_range[0]:
                    outfile.variables[q + Write_sts.RANGE][0] = q_values_min
                q_values_max = num.max(q_values)
                if q_values_max > q_range[1]:
                    outfile.variables[q + Write_sts.RANGE][1] = q_values_max



    def write_dynamic_quantities(self, outfile, quantities,
                    times, precis = netcdf_float, verbose = False):   
        """
            Write out given quantities to file.
        """
        for q in quantities:
            outfile.createVariable(q, precis, ('number_of_timesteps',
                                                      'number_of_points'))
            outfile.createVariable(q + Write_sts.RANGE, precis,
                                   ('numbers_in_range',))

            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q+Write_sts.RANGE][0] = max_float  # Min
            outfile.variables[q+Write_sts.RANGE][1] = -max_float # Max

        # Doing sts_precision instead of Float gives cast errors.
        outfile.createVariable('time', netcdf_float, ('number_of_timesteps',))

        if isinstance(times, (list, num.ndarray)):
            outfile.variables['time'][:] = times    # Store time relative

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('Statistics:')
            log.critical('    t in [%f, %f], len(t) == %d'
                         % (num.min(times), num.max(times), len(times.flat)))



def create_sts_boundary(sts_filename):
    """Create a list of points defining a boundary from an STS file.

    Create boundary segments from .sts file. Points can be stored in
    arbitrary order within the .sts file. The order in which the .sts points
    make up the boundary are given in order.txt file

    FIXME:
    Point coordinates are stored in relative eastings and northings.
    But boundary is produced in absolute coordinates
    """

    if sts_filename.endswith('.sts'):
        stsname_postfixed = sts_filename
    else:
        stsname_postfixed = sts_filename + '.sts'

    try:
        fid = NetCDFFile(stsname_postfixed, netcdf_mode_r)
    except IOError:
        msg = 'Cannot open %s' % stsname_postfixed
        raise IOError(msg)

    xllcorner = fid.xllcorner
    yllcorner = fid.yllcorner

    #Points stored in sts file are normalised to [xllcorner,yllcorner] but
    #we cannot assume that boundary polygon will be. At least the
    #additional points specified by the user after this function is called
    x = fid.variables['x'][:] + xllcorner
    y = fid.variables['y'][:] + yllcorner

    x = num.reshape(x, (len(x), 1))
    y = num.reshape(y, (len(y), 1))
    sts_points = num.concatenate((x,y), axis=1)

    return sts_points.tolist()



