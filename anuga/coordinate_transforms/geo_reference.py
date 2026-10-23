"""

"""


#FIXME: Ensure that all attributes of a georef are treated everywhere
#and unit test
import sys
import copy

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.anuga_exceptions import ANUGAError, TitleError, \
                                             ParsingError, ShapeError
from anuga.config import netcdf_float, netcdf_int, netcdf_float32
import anuga.utilities.log as log

import numpy as num


DEFAULT_ZONE = -1  # This signifies that simulation isn't located within a UTM framework.
                   # This is the case for hypothetical simulations or something relative
                   # to an arbitrary origin (e.g. the corner of a wavetank).

DEFAULT_PROJECTION = 'UTM'
DEFAULT_DATUM = 'wgs84'
DEFAULT_UNITS = 'm'

DEFAULT_SOUTHERN_FALSE_EASTING = 500000
DEFAULT_SOUTHERN_FALSE_NORTHING = 10000000

DEFAULT_NORTHERN_FALSE_EASTING = 500000
DEFAULT_NORTHERN_FALSE_NORTHING = 0

DEFAULT_HEMISPHERE = 'undefined'

# WGS84 UTM EPSG code ranges
# Northern hemisphere: 32601 (zone 1N) – 32660 (zone 60N)
# Southern hemisphere: 32701 (zone 1S) – 32760 (zone 60S)
_WGS84_UTM_NORTH_BASE = 32600
_WGS84_UTM_SOUTH_BASE = 32700

TITLE = '#geo reference' + "\n" # this title is referred to in the test format

class Geo_reference:
    """Coordinate reference system for an ANUGA domain.

    Attributes
    ----------
    zone : int
        UTM zone (1–60), or -1 (DEFAULT_ZONE) for non-UTM / arbitrary origin
        (e.g. a wave-tank simulation).
    hemisphere : str
        ``'southern'``, ``'northern'``, or ``'undefined'``.
    xllcorner : float
        X coordinate (easting) of the local origin relative to the UTM grid.
    yllcorner : float
        Y coordinate (northing) of the local origin relative to the UTM grid.
    datum : str
        Geodetic datum (default ``'wgs84'``).
    projection : str
        Map projection (default ``'UTM'``).
    units : str
        Units of measurement (default ``'m'``).
    false_easting : int
        False easting offset (500 000 m for WGS84 UTM).
    false_northing : int
        False northing offset (10 000 000 m southern, 0 m northern).
    epsg : int or None
        EPSG code for the coordinate reference system.  For WGS84 UTM this is
        auto-computed from *zone* and *hemisphere* (32600 + zone for northern,
        32700 + zone for southern).  Returns ``None`` when the zone is
        DEFAULT_ZONE or the CRS cannot be determined.
    """

    def __init__(self,
                 zone=None,
                 xllcorner=0.0,
                 yllcorner=0.0,
                 datum=DEFAULT_DATUM,
                 projection=DEFAULT_PROJECTION,
                 units=DEFAULT_UNITS,
                 false_easting=None,
                 false_northing=None,
                 hemisphere=DEFAULT_HEMISPHERE,
                 epsg=None,
                 NetCDFObject=None,
                 ASCIIFile=None,
                 read_title=None):
        """
        Parameters
        ----------
        zone : int, optional
            UTM zone (1–60) or -1 for no UTM framework.  Inferred from *epsg*
            when possible if not supplied.
        xllcorner : float, optional
            X (easting) of the local origin in metres.
        yllcorner : float, optional
            Y (northing) of the local origin in metres.
        datum : str, optional
            Geodetic datum.  Default ``'wgs84'``.
        projection : str, optional
            Map projection.  Default ``'UTM'``.
        units : str, optional
            Distance units.  Default ``'m'``.
        false_easting : int, optional
            Override the default false easting for the hemisphere.
        false_northing : int, optional
            Override the default false northing for the hemisphere.
        hemisphere : str, optional
            ``'southern'``, ``'northern'``, or ``'undefined'``.  Inferred from
            *epsg* when possible if not supplied.
        epsg : int, optional
            EPSG code.  For WGS84 UTM codes (32601–32660 northern,
            32701–32760 southern), *zone* and *hemisphere* are inferred
            automatically when they have not been set explicitly.
        NetCDFObject : file handle, optional
            Open NetCDF file to read geo-reference attributes from.
        ASCIIFile : file handle, optional
            Open text file to read geo-reference attributes from.
        read_title : str, optional
            Title line already read from *ASCIIFile* (pass if the caller has
            already consumed it).
        """

        if zone is None:
            zone = DEFAULT_ZONE

        self._epsg = None  # must exist before set_zone / set_hemisphere calls

        self.set_zone(zone)
        self.set_hemisphere(hemisphere)
        self.set_false_easting_northing(false_easting=false_easting, false_northing=false_northing)

        self.datum = datum
        self.projection = projection
        self.units = units
        self.xllcorner = float(xllcorner)
        self.yllcorner = float(yllcorner)

        if epsg is not None:
            self._set_epsg(int(epsg))

        if NetCDFObject is not None:
            self.read_NetCDF(NetCDFObject)

        if ASCIIFile is not None:
            self.read_ASCII(ASCIIFile, read_title=read_title)

        # Set flag for absolute points (used by get_absolute)
        # FIXME (Ole): It would be more robust to always use get_absolute()
        self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)

    def __eq__(self, other):

        # FIXME (Ole): Can this be automatically done for all attributes?
        # Anyway, it is arranged like this so one can step through and find out
        # why two objects might not be equal

        equal = True
        if self.false_easting != other.false_easting: equal = False
        if self.false_northing != other.false_northing: equal = False
        if self.datum != other.datum: equal = False
        if self.projection != other.projection: equal = False
        if self.zone != other.zone: equal = False
        if self.hemisphere != other.hemisphere: equal = False
        if self.units != other.units: equal = False
        if self.xllcorner != other.xllcorner: equal = False
        if self.yllcorner != other.yllcorner: equal = False
        if self.absolute != other.absolute: equal = False

        return(equal)

    def get_xllcorner(self):
        """Get the X coordinate of the origin of this georef."""
        return self.xllcorner

    def get_yllcorner(self):
        """Get the Y coordinate of the origin of this georef."""

        return self.yllcorner

    def set_zone(self, zone):
        """Set zone as an integer in [1,60] or -1 (DEFAULT_ZONE = unlocated).

        A negative zone in [-60, -2] is interpreted as a southern hemisphere
        shorthand: the zone number is taken as abs(zone) and hemisphere is set
        to 'southern' when it is currently 'undefined'.

        Note: zone=-1 is reserved for DEFAULT_ZONE (unlocated simulation) and
        is NOT interpreted as zone 1 southern hemisphere.  For zone 1 south,
        pass zone=1 with hemisphere='southern' explicitly.
        """
        zone = int(zone)

        if -60 <= zone <= -2:
            if self.hemisphere == 'undefined':
                self.hemisphere = 'southern'
            zone = abs(zone)

        assert (zone == -1 or (zone >= 1 and zone <= 60)), f'zone {zone} not valid.'

        self.zone = zone

    def get_zone(self):
        """Get the zone of this georef."""

        return self.zone

    def get_hemisphere(self):
        """Check if this georef has a defined hemisphere."""

        return self.hemisphere

    def set_hemisphere(self, hemisphere):

        msg = f"'{hemisphere}' not corresponding to allowed hemisphere values 'southern', 'northern' or 'undefined'"
        assert hemisphere in ['southern', 'northern', 'undefined'], msg

        self.hemisphere=str(hemisphere)

    def set_false_easting_northing(self, false_easting=None, false_northing=None):
        if false_easting is None:
            if self.hemisphere == 'southern':
                false_easting = DEFAULT_SOUTHERN_FALSE_EASTING
            elif self.hemisphere == 'northern':
                false_easting = DEFAULT_NORTHERN_FALSE_EASTING
            else:
                false_easting = 0
        if false_northing is None:
            if self.hemisphere == 'southern':
                false_northing = DEFAULT_SOUTHERN_FALSE_NORTHING
            elif self.hemisphere == 'northern':
                false_northing = DEFAULT_NORTHERN_FALSE_NORTHING
            else:
                false_northing = 0
        self.false_easting = int(false_easting)
        self.false_northing = int(false_northing)

    def _set_epsg(self, epsg):
        """Store EPSG code and infer zone/hemisphere for WGS84 UTM codes.

        For WGS84 UTM EPSG codes (32601–32660 northern, 32701–32760 southern),
        *zone* and *hemisphere* are inferred automatically when not already set.

        For all other EPSG codes (national grids, geographic CRS, etc.) the code
        is stored as-is.  If ``pyproj`` is available, *datum* and *projection*
        are populated from the CRS definition so the SWW metadata is accurate.

        Parameters
        ----------
        epsg : int
            EPSG code to store.
        """
        self._epsg = epsg
        used_pyproj = False

        if _WGS84_UTM_NORTH_BASE < epsg <= _WGS84_UTM_NORTH_BASE + 60:
            inferred_zone = epsg - _WGS84_UTM_NORTH_BASE
            inferred_hemisphere = 'northern'
        elif _WGS84_UTM_SOUTH_BASE < epsg <= _WGS84_UTM_SOUTH_BASE + 60:
            inferred_zone = epsg - _WGS84_UTM_SOUTH_BASE
            inferred_hemisphere = 'southern'
        else:
            # Not a WGS84 UTM code. Populate datum/projection/false easting-
            # northing from pyproj, and — for any UTM-based CRS on another datum
            # (e.g. GDA2020/GDA94 MGA EPSG 78xx/283xx, or NAD83 UTM) — infer the
            # UTM zone and hemisphere too. Genuinely non-UTM grids (British
            # National Grid 27700, RD New 28992, Albers 3577, geographic 4326)
            # yield (None, None) and keep zone/hemisphere unset.
            self._populate_from_pyproj(epsg)
            used_pyproj = True
            inferred_zone, inferred_hemisphere = self._utm_from_pyproj(epsg)

        if inferred_zone is None:
            return

        if self.zone == DEFAULT_ZONE:
            self.set_zone(inferred_zone)
        elif self.zone != inferred_zone:
            log.warning(f'EPSG {epsg} implies zone {inferred_zone} but zone '
                        f'{self.zone} is already set; zone unchanged.')

        if self.hemisphere == DEFAULT_HEMISPHERE:
            self.set_hemisphere(inferred_hemisphere)
            # The pyproj path already set false easting/northing precisely from
            # the CRS; only synthesise them for the WGS84 fast path.
            if not used_pyproj:
                self.set_false_easting_northing()
        elif self.hemisphere != inferred_hemisphere:
            log.warning(f'EPSG {epsg} implies {inferred_hemisphere} hemisphere but '
                        f'{self.hemisphere} is already set; hemisphere unchanged.')

    def _utm_from_pyproj(self, epsg):
        """Infer ``(zone, hemisphere)`` from a UTM-based EPSG code via pyproj.

        Recognises any UTM projection regardless of datum — WGS84 UTM, GDA2020/
        GDA94 MGA (EPSG 78xx / 283xx), NAD83 UTM, etc. — which the WGS84-only
        32601–32760 fast path in :meth:`_set_epsg` does not cover. PROJ reports
        these all as ``proj=utm`` with a ``zone``; the hemisphere is taken from
        the false northing (southern UTM uses 10 000 000 m), since the PROJ
        ``south`` flag is unreliable across the round-trip.

        Returns ``(None, None)`` for non-UTM CRS or when pyproj is unavailable.
        """
        try:
            from pyproj import CRS
            crs = CRS.from_epsg(epsg)
            d = crs.to_dict()
        except Exception:
            return None, None

        if d.get('proj') != 'utm' or d.get('zone') is None:
            return None, None
        try:
            zone = int(d['zone'])
        except (TypeError, ValueError):
            return None, None

        south = d.get('south')
        if south is None:
            try:
                fn = crs.to_cf().get('false_northing')
                south = fn is not None and float(fn) > 0
            except Exception:
                south = False
        return zone, ('southern' if south else 'northern')

    def _populate_from_pyproj(self, epsg):
        """Use pyproj to set CRS metadata from an EPSG code.

        Populates *projection*, *datum*, *false_easting*, and *false_northing*
        from the EPSG definition.  Does nothing if pyproj is not installed or
        the lookup fails.

        Parameters
        ----------
        epsg : int
            EPSG code to look up.
        """
        try:
            from pyproj import CRS
            crs = CRS.from_epsg(epsg)
        except Exception:
            return  # pyproj not available or invalid EPSG

        self.projection = crs.name

        if crs.datum is not None:
            self.datum = crs.datum.name
        elif crs.geodetic_crs is not None and crs.geodetic_crs.datum is not None:
            self.datum = crs.geodetic_crs.datum.name

        if not crs.is_projected:
            return

        # False easting/northing: try CF params first, fall back to PROJ dict.
        # Neither source covers all projections alone:
        #   to_cf()   works for BNG (27700) but not RD New (28992)
        #   to_dict() works for RD New (28992) but not some UTM variants
        fe = fn = None
        try:
            cf = crs.to_cf()
            fe = cf.get('false_easting')
            fn = cf.get('false_northing')
        except Exception:
            pass

        if fe is None or fn is None:
            try:
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    d = crs.to_dict()
                fe = d.get('x_0', fe)
                fn = d.get('y_0', fn)
            except Exception:
                pass

        if fe is not None:
            self.false_easting = int(round(float(fe)))
        if fn is not None:
            self.false_northing = int(round(float(fn)))

    @property
    def epsg(self):
        """EPSG code for this coordinate reference system.

        For WGS84 UTM projections the code is computed automatically from
        *zone* and *hemisphere* when it has not been set explicitly:

        * Northern hemisphere: ``32600 + zone`` (e.g. zone 55N → EPSG 32655)
        * Southern hemisphere: ``32700 + zone`` (e.g. zone 55S → EPSG 32755)

        Returns ``None`` when the zone is DEFAULT_ZONE (-1) or the CRS cannot
        be determined.

        Returns
        -------
        int or None
        """
        if self._epsg is not None:
            return self._epsg
        zone = self.zone
        if zone == DEFAULT_ZONE:
            return None
        # Negative zone implies southern hemisphere (legacy convention).
        if zone < 0:
            return _WGS84_UTM_SOUTH_BASE + abs(zone)
        if self.datum.lower() == 'wgs84' and self.projection.upper() == 'UTM':
            if self.hemisphere == 'southern':
                return _WGS84_UTM_SOUTH_BASE + zone
            if self.hemisphere == 'northern':
                return _WGS84_UTM_NORTH_BASE + zone
        return None

    @epsg.setter
    def epsg(self, value):
        if value is None:
            self._epsg = None
        else:
            self._set_epsg(int(value))

    def get_epsg(self):
        """Return the EPSG code for this coordinate reference system.

        Returns
        -------
        int or None
            EPSG code, or ``None`` if unknown.
        """
        return self.epsg

    def is_located(self):
        """Return True if this geo-reference describes a real, located CRS.

        A geo-reference is *located* when either:

        * *zone* is a valid UTM zone (1–60), or
        * an EPSG code has been set (covers national grids, geographic CRS, and
          any other projected CRS that does not use UTM zones, e.g.
          EPSG:28992 Netherlands RD New, EPSG:27700 British National Grid).

        A geo-reference with ``zone == DEFAULT_ZONE`` (-1) and no EPSG code is
        *unlocated* — this is the case for wavetank or other hypothetical
        simulations with an arbitrary local origin.

        Returns
        -------
        bool
        """
        return self.zone != DEFAULT_ZONE or self._epsg is not None

    def write_NetCDF(self, outfile):
        """Write georef attributes to an open NetCDF file.

        Parameters
        ----------
        outfile : file handle
            Handle to an open NetCDF file.
        """

        outfile.xllcorner = self.xllcorner
        outfile.yllcorner = self.yllcorner
        outfile.zone = self.zone
        outfile.hemisphere = self.hemisphere

        outfile.false_easting = self.false_easting
        outfile.false_northing = self.false_northing

        outfile.datum = self.datum
        outfile.projection = self.projection
        outfile.units = self.units

        epsg = self.epsg
        if epsg is not None:
            outfile.epsg = epsg

    def read_NetCDF(self, infile):
        """Set georef attributes from open NetCDF file.

        Parameters
        ----------
        infile : file handle
            Handle to an open NetCDF file.
        """

        self.xllcorner = float(infile.xllcorner)
        self.yllcorner = float(infile.yllcorner)
        self.zone = int(infile.zone)
        try:
            self.hemisphere = str(infile.hemisphere)
        except AttributeError:
            self.hemisphere = DEFAULT_HEMISPHERE

        self.false_easting = int(infile.false_easting)
        self.false_northing = int(infile.false_northing)

        self.datum = infile.datum
        self.projection = infile.projection
        self.units = infile.units

        # Read EPSG if present (old SWW files may not have it)
        try:
            self._epsg = int(infile.epsg)
        except AttributeError:
            self._epsg = None

        # Set flag for absolute points (used by get_absolute)
        # FIXME (Ole): It would be more robust to always use get_absolute()
        self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)

        if self.hemisphere == 'southern':
            if self.false_easting != DEFAULT_SOUTHERN_FALSE_EASTING:
                log.warning("WARNING: False easting of %f specified."
                             % self.false_easting)
                log.info("Default false easting is %f." % DEFAULT_SOUTHERN_FALSE_EASTING)
                log.info("ANUGA does not correct for differences in "
                             "False Eastings.")

            if self.false_northing != DEFAULT_SOUTHERN_FALSE_NORTHING:
                log.warning("WARNING: False northing of %f specified."
                             % self.false_northing)
                log.info("Default false northing is %f."
                             % DEFAULT_SOUTHERN_FALSE_NORTHING)
                log.info("ANUGA does not correct for differences in "
                             "False Northings.")

        if self.hemisphere == 'northern':
            if self.false_easting != DEFAULT_NORTHERN_FALSE_EASTING:
                log.warning("WARNING: False easting of %f specified."
                             % self.false_easting)
                log.info("Default false easting is %f." % DEFAULT_NORTHERN_FALSE_EASTING)
                log.info("ANUGA does not correct for differences in "
                             "False Eastings.")

            if self.false_northing != DEFAULT_NORTHERN_FALSE_NORTHING:
                log.warning("WARNING: False northing of %f specified."
                             % self.false_northing)
                log.info("Default false northing is %f."
                             % DEFAULT_NORTHERN_FALSE_NORTHING)
                log.info("ANUGA does not correct for differences in "
                             "False Northings.")



        # Suppress datum/projection warnings for non-UTM EPSG codes
        # (e.g. EPSG:28992 RD New has datum 'Amersfoort', not 'wgs84').
        # The EPSG code is the authoritative CRS identifier in that case.
        non_utm_epsg = (self._epsg is not None and
                        not (_WGS84_UTM_NORTH_BASE < self._epsg <= _WGS84_UTM_NORTH_BASE + 60 or
                             _WGS84_UTM_SOUTH_BASE < self._epsg <= _WGS84_UTM_SOUTH_BASE + 60))

        if not non_utm_epsg:
            if self.datum.upper() != DEFAULT_DATUM.upper():
                log.warning("WARNING: Datum of %s specified." % self.datum)
                log.info("Default Datum is %s." % DEFAULT_DATUM)
                log.info("ANUGA does not correct for differences in datums.")

            if self.projection.upper() != DEFAULT_PROJECTION.upper():
                log.warning("WARNING: Projection of %s specified."
                             % self.projection)
                log.info("Default Projection is %s." % DEFAULT_PROJECTION)
                log.info("ANUGA does not correct for differences in "
                             "Projection.")

        if self.units.upper() != DEFAULT_UNITS.upper():
            log.warning("WARNING: Units of %s specified." % self.units)
            log.info("Default units is %s." % DEFAULT_UNITS)
            log.info("ANUGA does not correct for differences in units.")



################################################################################
# ASCII files with geo-refs are currently not used
################################################################################

    def write_ASCII(self, fd):
        """Write georef attriutes to an open text file.

        fd  handle to open text file
        """

        fd.write(TITLE)
        fd.write(str(self.zone) + "\n")
        fd.write(str(self.xllcorner) + "\n")
        fd.write(str(self.yllcorner) + "\n")

    def read_ASCII(self, fd, read_title=None):
        """Set georef attribtes from open text file.

        fd  handle to open text file
        """

        try:
            if read_title is None:
                read_title = fd.readline()     # remove the title line
            if read_title[0:2].upper() != TITLE[0:2].upper():
                msg = ('File error.  Expecting line: %s.  Got this line: %s'
                       % (TITLE, read_title))
                raise TitleError(msg)
            self.zone = int(fd.readline())
            self.xllcorner = float(fd.readline())
            self.yllcorner = float(fd.readline())
        except SyntaxError:
            msg = 'File error.  Got syntax error while parsing geo reference'
            raise ParsingError(msg)

        # Fix some assertion failures
        if isinstance(self.zone, num.ndarray) and self.zone.shape == ():
            self.zone = self.zone[0]
        if (isinstance(self.xllcorner, num.ndarray) and
                self.xllcorner.shape == ()):
            self.xllcorner = self.xllcorner[0]
        if (isinstance(self.yllcorner, num.ndarray) and
                self.yllcorner.shape == ()):
            self.yllcorner = self.yllcorner[0]

        assert isinstance(self.xllcorner, float)
        assert isinstance(self.yllcorner, float)
        assert isinstance(self.zone, int)

################################################################################

    def change_points_geo_ref(self, points, points_geo_ref=None):
        """Change points to be absolute wrt new georef 'points_geo_ref'.

        points          the points to change
        points_geo_ref  the new georef to make points absolute wrt

        Returns the changed points data.
        If the points do not have a georef, assume 'absolute' values.
        """

        import copy

        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)

        # sanity checks
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            assert len(points) == 2, msg
            points = num.reshape(points, (1,2))

        msg = 'Points array must be two dimensional.\n'
        msg += 'I got %d dimensions' %len(points.shape)
        assert len(points.shape) == 2, msg

        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape
        assert points.shape[1] == 2, msg

        # FIXME (Ole): Could also check if zone, xllcorner, yllcorner
        # are identical in the two geo refs.
        if points_geo_ref is not self:
            # If georeferences are different
            points = copy.copy(points) # Don't destroy input
            if points_geo_ref is not None:
                # Convert points to absolute coordinates
                points[:,0] += points_geo_ref.xllcorner
                points[:,1] += points_geo_ref.yllcorner

            # Make points relative to primary geo reference
            points[:,0] -= self.xllcorner
            points[:,1] -= self.yllcorner

        if is_list:
            points = points.tolist()

        return points

    def is_absolute(self):
        """Test if points in georef are absolute.

        Return True if xllcorner==yllcorner==0 indicating that points in
        question are absolute.
        """

        # FIXME(Ole): It is unfortunate that decision about whether points
        # are absolute or not lies with the georeference object. Ross pointed this out.
        # Moreover, this little function is responsible for a large fraction of the time
        # using in data fitting (something in like 40 - 50%.
        # This was due to the repeated calls to allclose.
        # With the flag method fitting is much faster (18 Mar 2009).

        # FIXME(Ole): HACK to be able to reuse data already cached (18 Mar 2009).
        # Remove at some point
        if not hasattr(self, 'absolute'):
            self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)

        # Return absolute flag
        return self.absolute

    def get_absolute(self, points):
        """Given a set of points geo referenced to this instance, return the
        points as absolute values.
        """

        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)
        if len(points.shape) == 1:
            # One point has been passed
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError(msg)


        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape
        if not points.shape[1] == 2:
            raise ShapeError(msg)


        # Add geo ref to points
        if not self.is_absolute():
            points = copy.copy(points) # Don't destroy input
            points[:,0] += self.xllcorner
            points[:,1] += self.yllcorner


        if is_list:
            points = points.tolist()

        return points

    def get_relative(self, points):
        """Convert points to relative measurement.

        points Points to convert to relative measurements

        Returns a set of points relative to the geo_reference instance.

        This is the inverse of get_absolute().
        """

        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError(msg)

        if not points.shape[1] == 2:
            msg = ('Input must be an N x 2 array or list of (x,y) values. '
                   'I got an %d x %d array' % points.shape)
            raise ShapeError(msg)

        # Subtract geo ref from points
        if not self.is_absolute():
            points = copy.copy(points) # Don't destroy input
            points[:,0] -= self.xllcorner
            points[:,1] -= self.yllcorner

        if is_list:
            points = points.tolist()

        return points

    def reconcile_zones(self, other):

        if other is None:
            # FIXME(Ole): Why would we do this?
            other = Geo_reference()
            #raise Exception('Expected georeference object, got None')

        if (self.zone == other.zone):
            pass
        elif self.zone == DEFAULT_ZONE:
            self.zone = other.zone
        elif other.zone == DEFAULT_ZONE:
            other.zone = self.zone
        else:
            msg = ('Geospatial data must be in the same '
                   'ZONE to allow reconciliation. I got zone %d and %d'
                   % (self.zone, other.zone))
            raise ANUGAError(msg)

        # Should also reconcile hemisphere
        if (self.hemisphere == other.hemisphere):
            pass
        elif self.hemisphere == DEFAULT_HEMISPHERE:
            self.hemisphere = other.hemisphere
        elif other.hemisphere == DEFAULT_HEMISPHERE:
            other.hemisphere = self.hemisphere
        else:
            msg = ('Geospatial data must be in the same '
                   'HEMISPHERE to allow reconciliation. I got hemisphere %d and %d'
                   % (self.hemisphere, other.hemisphere))
            raise ANUGAError(msg)

    # FIXME (Ole): Do we need this back?
    #def easting_northing2geo_reffed_point(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    #def easting_northing2geo_reffed_points(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    def get_origin(self):
        """Get origin of this geo_reference."""

        return (self.zone, self.xllcorner, self.yllcorner)

    def __repr__(self):
        epsg = self.epsg
        is_utm = (epsg is not None and
                  (_WGS84_UTM_NORTH_BASE < epsg <= _WGS84_UTM_NORTH_BASE + 60 or
                   _WGS84_UTM_SOUTH_BASE < epsg <= _WGS84_UTM_SOUTH_BASE + 60))
        if epsg is not None and not is_utm:
            # Non-UTM CRS: zone/hemisphere are not meaningful; show CRS name instead
            return ('(crs=%s, easting=%f, northing=%f, epsg=%i)'
                    % (self.projection, self.xllcorner, self.yllcorner, epsg))
        if epsg is not None:
            return ('(zone=%i, easting=%f, northing=%f, hemisphere=%s, epsg=%i)'
                    % (self.zone, self.xllcorner, self.yllcorner, self.hemisphere, epsg))
        return ('(zone=%i, easting=%f, northing=%f, hemisphere=%s)'
                % (self.zone, self.xllcorner, self.yllcorner, self.hemisphere))

    #def __cmp__(self, other):
    #    """Compare two geo_reference instances.#
    #
    #    self   this geo_reference instance
    #    other  another geo_reference instance to compare against#
    #
    #    Returns 0 if instances have the same attributes, else returns 1.
    #
    #    Note: attributes are: zone, xllcorner, yllcorner.
    #    """

    #    # FIXME (DSG) add a tolerence
    #   if other is None:
    #        return 1
    #    cmp = 0
    #    if not (self.xllcorner == self.xllcorner):
    #        cmp = 1
    #    if not (self.yllcorner == self.yllcorner):
    #        cmp = 1
    #    if not (self.zone == self.zone):
    #        cmp = 1
    #    return cmp


def write_NetCDF_georeference(georef, outfile):
    """Write georeference info to a NetCDF file, usually a SWW file.

    georef   a georef instance or parameters to create a georef instance
    outfile  path to file to write

    Returns the normalised georef.
    """

    geo_ref = ensure_geo_reference(georef)
    geo_ref.write_NetCDF(outfile)
    return geo_ref


def ensure_geo_reference(origin):
    """Create a georef object from a tuple of attributes.

    origin  a georef instance or (zone, xllcorner, yllcorner)

    If origin is None, return None, so calling this function doesn't
    effect code logic.
    """

    if isinstance(origin, Geo_reference):
        geo_ref = origin
    elif origin is None:
        geo_ref = None
    else:
        if len(origin) == 1:
            geo_ref = Geo_reference(zone = origin)
        elif len(origin) == 2:
            geo_ref = Geo_reference(zone = -1, xllcorner=origin[0], yllcorner=origin[1])
        elif len(origin) == 3:
            geo_ref = Geo_reference(zone = origin[0], xllcorner=origin[1], yllcorner=origin[2])
        else:
            raise Exception(f'Invalid input {origin}, expected (zone), (xllcorner, yllcorner), or (zone, xllcorner, yllcorner).')


    return geo_ref


#-----------------------------------------------------------------------

if __name__ == "__main__":
    pass
