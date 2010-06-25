

def load_ungenerate(ofile):
    """
    import a file, ofile, with the format
    [poly]
    poly format:
    First line:  <# of vertices> <x centroid> <y centroid>
    Following lines: <x> <y> 
    last line:  "END"

    Note: These are clockwise.
    
    Returns a dict containing "points", "segments", and "polygons".
    """
    fd = open(ofile,'r')
    Dict = readUngenerateFile(fd)
    fd.close()
    return Dict


def readUngenerateFile(fd):
    """
    import a file, ofile, with the format
    [poly]
    poly format:
    First line:  <# of polynomial> <x centroid> <y centroid>
    Following lines: <x> <y> 
    last line:  "END"
    """
    
    END_DELIMITER = 'END'
    
    points = []
    segments = []
    polygons = []
    
    isEnd = False
    line = fd.readline() #not used <# of polynomial> <x> <y>
    while not isEnd:
        poly = []
        line = fd.readline()
        fragments = line.split()
        x = float(fragments.pop(0))
        y = float(fragments.pop(0))
        points.append([x,y])
        poly.append([x,y])
        PreviousVertIndex = len(points)-1
        firstVertIndex = PreviousVertIndex
        
        line = fd.readline() #Read the next line
        while not line.startswith(END_DELIMITER): 
            #print "line >" + line + "<"
            fragments = line.split()
            x = float(fragments.pop(0))
            y = float(fragments.pop(0))
            points.append([x,y])
            poly.append([x,y])
            thisVertIndex = len(points)-1
            segment = [PreviousVertIndex,thisVertIndex]
            segments.append(segment)
            PreviousVertIndex = thisVertIndex
            line = fd.readline() #Read the next line
        # If the last and first segments are the same,
        # Remove the last segment and the last vertex
        # then add a segment from the second last vert to the 1st vert
        thisVertIndex = len(points)-1
        firstVert = points[firstVertIndex]
        thisVert = points[thisVertIndex]
        #print "firstVert",firstVert
        #print "thisVert",thisVert
        if (firstVert[0] == thisVert[0] and firstVert[1] == thisVert[1]):
            points.pop()
            segments.pop()
            poly.pop()
            thisVertIndex = len(points)-1
        segments.append([thisVertIndex, firstVertIndex])
        
        line = fd.readline() # read <# of polynomial> <x> <y> OR END
        #print "line >>" + line + "<<"
        # do poly stuff here
        polygons.append(poly)
        if line.startswith(END_DELIMITER):
            isEnd = True
    
    #print "points", points       
    #print "segments", segments
    ungenerated_dict = {}
    ungenerated_dict['points'] = points
    ungenerated_dict['segments'] = segments
    ungenerated_dict['polygons'] = polygons
    return ungenerated_dict
