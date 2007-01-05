#!/usr/bin/env python

import sys

from types import ListType, TupleType

import anuga.mesh_engine.mesh_engine_c_layer as triang

from Numeric import array, Float, Int

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.anuga_exceptions import ANUGAError
    
def generate_mesh(points=None,
                  segments=None,holes=None,regions=None,
                  pointatts=None,segatts=None,
                  mode=None, dummy_test=None):
    """
    
    """
    #FIXME (DSG-DSG): Catch parameters that are lists,
    #instead of lists of lists
    # check shape[1] is 2 etc

    if points is None:
        points = []

    if segments is None:
        segments = []

    if holes is None:
        holes = []
        
    if regions is None:
        regions = []

    if dummy_test is None:
        dummy_test  = []
        
    try:
        points =  ensure_numeric(points, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent points array.'
        raise ANUGAError, msg
    if points.shape[1] <>2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError, msg

    # This is after points is numeric
    if pointatts is None or pointatts == []:
        pointatts = [[] for x in range(points.shape[0])]
        
    try:
        segments = ensure_numeric(segments, Int)
    except ValueError:
        msg = 'ERROR: Inconsistent segments array.'
        raise ANUGAError, msg
    
    # This is after segments is numeric
    if segatts is None or segatts == []:
        segatts = [0 for x in range(segments.shape[0])]
        
    try:
        holes = ensure_numeric(holes, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent holess array.'
        raise ANUGAError, msg

   
    regions = add_area_tag(regions)
    try:
        regions = ensure_numeric(regions, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent regions array.'
        raise ANUGAError, msg
        
    if not regions.shape[0] == 0 and regions.shape[1] <= 2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError, msg
    
    try:
        pointatts = ensure_numeric(pointatts, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError, msg

    if pointatts.shape[0] <> points.shape[0]:
        msg = """ERROR: Point attributes array not the same shape as
        point array."""
        raise ANUGAError, msg
    
    try:
        segatts = ensure_numeric(segatts, Int)
    except ValueError:
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError, msg
    if segatts.shape[0] <> segments.shape[0]:
        msg = """ERROR: Segment attributes array not the same shape as
        segment array."""
        raise ANUGAError, msg
    return triang.genMesh(points,segments,holes,regions,
                          pointatts,segatts, mode, points)

def add_area_tag(regions):
    """
    So, what is the format?
    A list with
    [x,y,region_tag,area] OR [x,y,region_tag]
    if it's [x,y,region_tag], add a 4th element, value of 0.0.
    """
    if isinstance(regions, ListType):
        for i, region in enumerate(regions):
            if len(region) == 3:
                if isinstance(region, TupleType):
                    #FIXME: How do you convert a tuple to a list?
                    # I can do it a stupid way..
                    tuple = region[:]
                    regions[i] = []
                    for j in tuple:
                        regions[i].append(j)
                    regions[i].append(0.0)
                else:
                    regions[i].append(0.0)
                    
                # let ensure numeric catch this..    
                #len(region) <= 2:
                #msg = 'ERROR: Inconsistent regions array.'
                #raise msg
            #elif
    return regions

if __name__ == "__main__":
    pass 
