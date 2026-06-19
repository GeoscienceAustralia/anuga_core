#!/usr/bin/env python

# Lat Long - UTM, UTM - Lat Long conversions
#
# see http://www.pygps.org
#


from math import pi, sin, cos, tan, sqrt

#LatLong- UTM conversion..h
#definitions for lat/long to UTM and UTM to lat/lng conversions
#include <string.h>

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 2
_eccentricitySquared = 3

_ellipsoid = [
#  id, Ellipsoid name, Equatorial Radius, square of eccentricity
# first once is a placeholder only, To allow array indices to match id numbers
	[ -1, "Placeholder", 0, 0],
	[ 1, "Airy", 6377563, 0.00667054],
	[ 2, "Australian National", 6378160, 0.006694542],
	[ 3, "Bessel 1841", 6377397, 0.006674372],
	[ 4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
	[ 5, "Clarke 1866", 6378206, 0.006768658],
	[ 6, "Clarke 1880", 6378249, 0.006803511],
	[ 7, "Everest", 6377276, 0.006637847],
	[ 8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
	[ 9, "Fischer 1968", 6378150, 0.006693422],
	[ 10, "GRS 1967", 6378160, 0.006694605],
	[ 11, "GRS 1980", 6378137, 0.00669438],
	[ 12, "Helmert 1906", 6378200, 0.006693422],
	[ 13, "Hough", 6378270, 0.00672267],
	[ 14, "International", 6378388, 0.00672267],
	[ 15, "Krassovsky", 6378245, 0.006693422],
	[ 16, "Modified Airy", 6377340, 0.00667054],
	[ 17, "Modified Everest", 6377304, 0.006637847],
	[ 18, "Modified Fischer 1960", 6378155, 0.006693422],
	[ 19, "South American 1969", 6378160, 0.006694542],
	[ 20, "WGS 60", 6378165, 0.006693422],
	[ 21, "WGS 66", 6378145, 0.006694542],
	[ 22, "WGS-72", 6378135, 0.006694318],
	[ 23, "WGS-84", 6378137, 0.00669438]
]

#Reference ellipsoids derived from Peter H. Dana's website-
#http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
#Department of Geography, University of Texas at Austin
#Internet: pdana@mail.utexas.edu
#3/22/95

#Source
#Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
#1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency

def LLtoUTM(Lat, Long, ReferenceEllipsoid=23):
    """Convert latitude/longitude to UTM coordinates (WGS84).

    Parameters
    ----------
    Lat : float
        Latitude in decimal degrees (positive north).
    Long : float
        Longitude in decimal degrees (positive east).
    ReferenceEllipsoid : int, optional
        Ignored — retained for backward compatibility. pyproj always uses WGS84.

    Returns
    -------
    ZoneNumber : int
    UTMEasting : float
    UTMNorthing : float
    """
    # Compute UTM zone number, retaining special-zone handling for Norway/Svalbard.
    LongTemp = (Long + 180) - int((Long + 180) / 360) * 360 - 180
    ZoneNumber = int((LongTemp + 180) / 6) + 1
    if 56.0 <= Lat < 64.0 and 3.0 <= LongTemp < 12.0:
        ZoneNumber = 32
    if 72.0 <= Lat < 84.0:
        if   LongTemp <  9.0: ZoneNumber = 31
        elif LongTemp < 21.0: ZoneNumber = 33
        elif LongTemp < 33.0: ZoneNumber = 35
        elif LongTemp < 42.0: ZoneNumber = 37

    from anuga.coordinate_transforms.redfearn import ll_to_epsg
    epsg = 32700 + ZoneNumber if Lat < 0 else 32600 + ZoneNumber
    easting, northing = ll_to_epsg(Lat, Long, epsg)
    return (ZoneNumber, float(easting), float(northing))


def _UTMLetterDesignator(Lat):
#This routine determines the correct UTM letter designator for the given latitude
#returns 'Z' if latitude is outside the UTM limits of 84N to 80S
#Written by Chuck Gantz- chuck.gantz@globalstar.com

    if 84 >= Lat >= 72: return 'X'
    elif 72 > Lat >= 64: return 'W'
    elif 64 > Lat >= 56: return 'V'
    elif 56 > Lat >= 48: return 'U'
    elif 48 > Lat >= 40: return 'T'
    elif 40 > Lat >= 32: return 'S'
    elif 32 > Lat >= 24: return 'R'
    elif 24 > Lat >= 16: return 'Q'
    elif 16 > Lat >= 8: return 'P'
    elif  8 > Lat >= 0: return 'N'
    elif  0 > Lat >= -8: return 'M'
    elif -8> Lat >= -16: return 'L'
    elif -16 > Lat >= -24: return 'K'
    elif -24 > Lat >= -32: return 'J'
    elif -32 > Lat >= -40: return 'H'
    elif -40 > Lat >= -48: return 'G'
    elif -48 > Lat >= -56: return 'F'
    elif -56 > Lat >= -64: return 'E'
    elif -64 > Lat >= -72: return 'D'
    elif -72 > Lat >= -80: return 'C'
    else: return 'Z'	# if the Latitude is outside the UTM limits

def UTMtoLL(northing, easting, zone, isSouthernHemisphere=True,
            ReferenceEllipsoid=23):
    """Convert UTM coordinates to latitude/longitude (WGS84).

    Parameters
    ----------
    northing : float
        UTM northing in metres.
    easting : float
        UTM easting in metres.
    zone : int
        UTM zone number (1–60).
    isSouthernHemisphere : bool, optional
        True (default) for southern hemisphere, False for northern.
    ReferenceEllipsoid : int, optional
        Ignored — retained for backward compatibility. pyproj always uses WGS84.

    Returns
    -------
    lat : float
    lon : float
    """
    from anuga.coordinate_transforms.redfearn import epsg_to_ll
    epsg = 32700 + int(zone) if isSouthernHemisphere else 32600 + int(zone)
    lat, lon = epsg_to_ll(easting, northing, epsg)
    return float(lat), float(lon)

if __name__ == '__main__':
    (z, e, n) = LLtoUTM(-45.00, -75.00, 23)
    print(z, e, n)
    (lat, lon) = UTMtoLL(n, e, z, 23)
    print(lat, lon)
