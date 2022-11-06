"""point.py - Represents a generic point on a sphere as a Python object.

   See documentation of class Point for details.
   Ole Nielsen, ANU 2002
"""

from math import cos, sin, pi

def acos(c):
  """acos -  Safe inverse cosine
  
     Input argument c is shrunk to admissible interval
     to avoid case where a small rounding error causes
     a math domain error.
  """
  from math import acos

  if c > 1: c=1
  if c < -1: c=-1

  return acos(c)




class Point(object):
    """Definition of a generic point on the sphere.
  
    Defines a point in terms of latitude and longitude
    and computes distances to other points on the sphere.

    Initialise as
      Point(lat, lon), where lat and lon are in decimal degrees (dd.dddd)
      

    Public Methods: 
        distance_to(P)
        bearing_to(P)
        dist(P)

    Author: Ole Nielsen, ANU 2002        
    """
     
    # class constants
    R = 6372000  # Approximate radius of Earth (m)
    degrees2radians = pi/180.0
    

    def __init__(self, latitude, longitude):

        assert(latitude >= -90 and latitude <= 90.0)
        assert(longitude >= -180 and longitude <= 180.0)

        self.latitude =  float(latitude)
        self.longitude = float(longitude) 

        lat = latitude * self.degrees2radians    # Converted to radians
        lon = longitude * self.degrees2radians   # Converted to radians
        self.coslat = cos(lat)
        self.coslon = cos(lon)
        self.sinlat = sin(lat)
        self.sinlon = sin(lon)

    def BearingTo(self,P):
        """ Bearing (in degrees) to point P"""
        AZ = self.AZ(P)
        return int(round(AZ/self.degrees2radians))

    
    def DistanceTo(self,P):
        """ Distance to point P"""           
        GCA = self.GCA(P)
        return self.R*GCA


    def Dist(self,P):
        """ Very cheap and rough approximation to distance"""
        return max(abs(self.latitude-P.latitude),abs(self.longitude-P.longitude))


    #--------------------------------------------------------------------------

    def __repr__(self):
        d = 2
        lat = round(self.latitude,d)
        lon = round(self.longitude,d)
        return ' (' + str(lat)+ ', '+ str(lon)+')'


    def GCA(self,P):
        """ Compute the Creat Circle Angle (GCA) between current point and P.
        """
        
        alpha = P.coslon*self.coslon + P.sinlon*self.sinlon
        # The original formula is alpha = cos(self.lon - P.lon)
        # but rewriting lets us make us of precomputed trigonometric values.
        
        x = alpha*self.coslat*P.coslat + self.sinlat*P.sinlat
        return acos(x)

    
    def AZ(self,P):
        """ Compute Azimuth bearing (AZ) from current point to P"""

        # Compute cosine(AZ), where AZ is the azimuth angle
        GCA = self.GCA(P)
        c = P.sinlat - self.sinlat*cos(GCA)
        c = c/self.coslat/sin(GCA)

        AZ = acos(c)

        # Reverse direction if bearing is westward, i.e. sin(self.lon - P.lon) > 0
        # Without this correction the bearing due west, say, will be 90 degrees
        # because the formulas work in the positive direction which is east.
        #
        # Precomputed trigonometric values are used to rewrite the formula:

        if self.sinlon*P.coslon - self.coslon*P.sinlon > 0: AZ = 2*pi - AZ
    
        return AZ



#-----------------------------------------------------------------------

if __name__ == "__main__":
    pass
