"""
 FILE :  MODEL_BUILD_TOOLS_04.py
 DATE:  28/05/2013
 
 Last Change:  Can now define directories by using underscore _ 
 so mannings = 0_05



This script is meant to work with a set of standardised directories so that an ANUGA model script
can be kept extremely minimalist, this module provides functions that will semi-automate the population
of polygons and attributes to:

1. define MESH Refinement
2. define Variation in ROughness in the model
3. define the additional of buildings of various heights
4. define rainfall spatial variation via polygons

To ensure these functions operate it is important to use the
STANDARD_MODEL_CREATE_SCRIPT that will create a set of standardised directories to store model data

The concept is to populate this standard directory structure with the required 'csv' files 
describing the required polygons.

MESH REFINE
EG: c:\Model\Data\Mesh_Refine\10000 and c:\Model\Data\Mesh-Refine\1000
Contain the polygons that describe the mesh refine for triangles no larger 
than 10,000m2 and 1,000m2

BUILDINGS
Similarly buildings within numbered building directories will signify the building height
EG: c:\Model\Data\Buildings\8  and c:\Model\Data\Buildings\12
will contain polygons to describe buildings that are 8 and 12 metres high

ROUGHNESS
The directory structure will be based on identifying for example the roughness directory
under which numbered directories will signify the rougness value.
EG: c:\Model\Data\Roughness\100 and c:\Model\Data\Roughness\10
will contain polygons top describe surface roughness in the model of 0.100 and 0.010

RAINFALL
EG: c:\Model\Data\Rainfall_Polys\Gauge_name_01 and c:\Model\Data\Rainfall_Polys\Gauge_name_02
Will list the polygons that apply rainfall described by the raingauge files in the directory name
options: Either the directory Name is the Raingauge file name, or the Raingauge file (with extension tms) 
is located in the directory and associated with the polygons in the same directory.... or
the raingauge is kept in another directory under FORCEFUNC and there is a pointer to it ???

TO BE RESOLVED 

METHOD OF CALLING:
Identify the directory in which sub directories will appear that have names as numeric values

USAGE FOR REFINE (example):
mesh_refine_directory =join('POLYS','Cycleway')
mesh_refine_defined = get_REFINE_polygon_value_list(mesh_refine_directory) # Returns List of Polygons and Value from Directory Name
print 'Mesh Refine Definitions in Directory....'
print mesh_refine_defined
# Need to Reverse the Lists so that mesh Refine is applied from Coarse to Fine. As default directory read is from smallest to largest number....
mesh_refine_defined.reverse()
print 'Mesh Refine REVERSED HERE....'
#mesh_refine_reversed = mesh_refine_defined.reverse()
print mesh_refine_defined


buildings_directory =join('02POLYS','04_BLDGS')
buildings_defined = get_BUILDING_polygon_value_list(buildings_directory)

such that all polygons under a directory named 45, under directory:- 04_BLDGS  will be assigned a height of 4.5m

Need to add the following lines to Scripts:


from anuga.utilities.model_tools import get_polygon_list_from_files
from anuga.utilities.model_tools import get_polygon_dictionary
from anuga.utilities.model_tools import get_REFINE_polygon_value_list
from anuga.utilities.model_tools import get_ROUGHNESS_polygon_value_list
from anuga.utilities.model_tools import get_BUILDING_polygon_value_list
"""

import os
import glob
import anuga.geometry.polygon.read_polygon as read_polygon


def get_polygon_list_from_files(dir):
    """Read all polygons found in specified dir and return them in a list
       Called by:
       get_polygon_dictionary
       Purpose:
       To fill a list with all of the polygons read under a specified directory
       CALLS:
       anuga.utilities.polygon.read_polygon
    """
    
    from anuga.geometry.polygon import read_polygon
    print 'Reading polygons from ' + dir
    polylist = []
    for filename in os.listdir(dir):
        print filename
        polylist.append(read_polygon(os.path.join(dir, filename)))
    return polylist

def get_polygon_dictionary(dir):
    """Create dictionary of polygons with directory names 
       indicating associated attribute values 
       Called by:
       get_polygon_value_list
       Purpose:
       To Fill a Dictionary with sets of poygons and attribute, from a list of polygons 
       and using the directory name as the attribute
       For Example used to read Mesh Size Directory 1500, using all polygons in the directory
       to create mesh refinement to 1500m2
       CALLS:
       get_polygon_list_from_files
    """
    
    try:
        attribute_values = os.listdir(dir)  # Create the Attribute from the Directory Name
    except:
        msg = 'Directory %s was not found' % dir
        raise Exception(msg)
    D = {}   # Create Empty Dictionary
    for a in attribute_values:
        # If a ROUGNESS Directory Need to Devide value by 1000   How can this be done ???   Or Do I create a New function ??
        # Similarly if Rainfall File Name how is it handled ??
        D[a] = get_polygon_list_from_files(os.path.join(dir, a)) # Fill Item [a] in the Dictionary with FIle name and attribute
    return D

# Create seperate Functions for get_REFINE_polygon_value_list   and get_ROUGHNESS_polygon_value_list

def get_REFINE_polygon_value_list(dir):
    """Create MESH Refinement Polygon-value list with values obtained from directory names 
       that is a List of Polygons attributed with the Value read from the directory name...
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    
    print 'Mesh Refinement Triangular Area Values are assigned from Directory name with underscore being decimal point'    
    print 'Note So using 0_10 limits Triangles to 0.1m2 etc..'
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            print 'Target Mesh Cell Area = ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    return polygon_value_list
    
def get_ROUGHNESS_polygon_value_list(dir):
    """Create Mannings ROUGHNESS Polygon-value list with values obtained from directory names (x1000)
       that is a List of Polygons attributed with the Value read from the directory name...
        For this case of ROUGHNESS the directory name should be called Mannings x 1000.0
        so that 0.015 has a directory name of 15, etc....
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """

    print 'Roughness Values are assigned from Directory name :'
    print 'so that 0_085 = mannings rougness of 0.085'
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            print 'Mannings Roughness = ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    return polygon_value_list    

def get_STAGE_polygon_value_list(dir):
    """Create STAGE Polygon-value list with values obtained from directory names (/100.0)
       that is a List of Polygons attributed with the Value read from the directory name...
        For this case of STAGE the directory name should be called Stage Level x 100.0
        so that a directory name of 45 results in a STAGE Level of 0.45m, etc....       
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    print 'STAGE LEVEL Values are assigned from Directory name :'
    print 'So with 0_36 directory name = 0.36m in height'
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            print 'Stage Set to = ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    return polygon_value_list

    
def get_BUILDING_polygon_value_list(dir):
    """Create BUILDING Polygon-value list with values obtained from directory names (/10.0)
       that is a List of Polygons attributed with the Value read from the directory name...
        For this case of BUILDINGS the directory name should be called Building height x 10.0
        so that a directory name of 45 results in a building height of 4.5m, etc....       
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """

    print 'BUILDING HEIGHT Values are assigned from Directory name:'
    print 'So with directory name 36_5 height = 36.5m '
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            print 'Building Height set to ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    return polygon_value_list
    
    
def get_POLYS_from_Multi_polyfile(dir):
    """Create List of Polygons from a Directory with a File containing Multiple-Polygons
       
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    print 'Getting File with Multiple POLYS from Directory:'
    print 'The one file will be used to create multiple Polys for ANUGA'
    filepattern='*.xyz'   # get a list of matching filenames in the directory, MAKE SURE only 1 csv file is in the DIR
    pattern = os.path.join(dir, filepattern)
    holes_file_list = glob.glob(pattern)  # List of Files has only 1 in it
    print 'Holes file_List:-',holes_file_list
    for fd in holes_file_list:
        holes_file=fd
        print 'holes file:-',holes_file
        infid = open(holes_file, 'r')
    # Got the Multi Poly file now process
    polylist = []
    check_pts_list=[]
    Poly_count=0
    lines = infid.readlines() # Reads ALL Lines in file infid
    total_lines_in_file= len(lines)
    print "total number of lines in the Holes FILE is: ",total_lines_in_file
    for i, line in enumerate(lines): # ==================================================== FOR LOOP ===========================
        if line.strip().startswith('####'):
            check_pts_list=[]
            if Poly_count==0:
                pass
            else:
                polylist.append(polygon)
                outfid.close()
            Poly_count+=1
            # Create A poly File for each polygon defined in the Multi-Poly file
            print 'Polygon #',Poly_count
            path_poly=os.path.dirname(os.path.dirname(holes_file))
            print path_poly
            poly_write_file="Poly_"+str(Poly_count)+".csv"
            outfid = open(os.path.join(path_poly,poly_write_file), 'w')
            polygon=[]
            # Instead write to a List
        else:
            print line
            fields = line.split(',')
            polygon.append([float(fields[0]),float(fields[1])])
            #line.rstrip('\n'))  # Add Points [x,y]
            if line in check_pts_list:
                pass
            else:
                outfid.write(line)
                check_pts_list.append(line)
            
    outfid.close()  
    polylist.append(polygon)    # Add polygon to the list of Polys
    #print polylist
    return polylist    

def get_POLYS_from_Mid_Mif(dir):
    """Create List of Polygons from a Directory with a File containing Multiple-Polygons
       
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    print 'Getting File with Multiple POLYS from Directory:'
    print 'The one file will be used to create multiple Polys for ANUGA'
    filepattern='*.mif'   # get a list of matching filenames in the directory, MAKE SURE only 1 csv file is in the DIR
    pattern = os.path.join(dir, filepattern)
    holes_file_list = glob.glob(pattern)  # List of Files has only 1 in it
    print 'From DIR: ',dir
    print 'Holes file_List:-',holes_file_list
    for fd in holes_file_list:
        holes_file=fd
        print 'holes file:-',holes_file
        infid = open(holes_file, 'r')
    # Got the Multi Poly file now process
    polylist = []
    check_pts_list=[]
    Poly_count=0
    Poly_line_count=0
    Points_in_Poly = 100
    lines = infid.readlines() # Reads ALL Lines in file infid
    total_lines_in_file= len(lines)
    print "total number of lines in the Holes FILE is: ",total_lines_in_file
    for i, line in enumerate(lines): # ==================================================== FOR LOOP ===========================


        if line.strip().startswith('Region'):

            Poly_line_count=0
            check_pts_list=[]
            if Poly_count==0:
                pass
            else:
                polylist.append(polygon)
                outfid.close()
            Poly_count+=1
            # Create A poly File for each polygon defined in the Multi-Poly file
            print 'Polygon #',Poly_count
            path_poly=os.path.dirname(os.path.dirname(holes_file))
            print path_poly
            poly_write_file="Poly_"+str(Poly_count)+".csv"
            outfid = open(os.path.join(path_poly,poly_write_file), 'w')
            polygon=[]
            # Instead write to a List
        elif line.strip().startswith('    Pen'):
            pass
        elif line.strip().startswith('    Brush'):
            pass
        else:
            Poly_line_count+=1
            if Poly_line_count > 1 and Poly_line_count <= (Points_in_Poly+1) and Poly_count<>0:
                print line, Points_in_Poly,Poly_line_count
                fields = line.split(' ')
                polygon.append([float(fields[0]),float(fields[1])])
                #line.rstrip('\n'))  # Add Points [x,y]
                if line in check_pts_list:   # Get rid of any doubled up points
                    pass
                else:
                    outfid.write("%.3f,%.3f\n" % (float(fields[0]),float(fields[1])))
 
                    check_pts_list.append(line)
            elif Poly_line_count==1 and Poly_count<>0:
                # read number of points in poly
                print 'line=',line
                Points_in_Poly=int(line)
                print 'Number Points in Poly =',Points_in_Poly
                
        #raw_input('Check Output...')            
    outfid.close()  
    
    polylist.append(polygon)    # Add polygon to the list of Polys
    print polylist
    return polylist    




"""
Here are the tools that Peter Milevski  has contributed.
"""







def read_polygon_dir(weight_dict, directory, filepattern='*.csv'):
    """
    In a directory directory looks at all files matching filepattern
    and returns a list of tuples consisting of polygon and a weight 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # check that the dictionary contains *all* the files
    
    errors = []
    for f in files:
        try:
            _ = weight_dict[f]
        except KeyError:
            errors.append(f)
            
    if errors:
        msg = ''
        for f in errors:
            msg = msg + ', ' + f
        raise KeyError, 'Files not defined in dictionary: %s' % msg[2:]

    # now get the result list
    result = []
    for f in files:
        result.append((read_polygon(f), weight_dict[f]))
    return result



#define a function with without an attribute
def read_hole_dir(directory, filepattern='*.csv'):
    """
    In a directory directory looks at all files matching filepattern
    and returns a list of polygon 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # now get the result list
    result = []
    for f in files:
        result.append(read_polygon(f))
    return result
    
    
# Define a function to read Single File with Multi-polygons
def read_multi_poly_file(multi_P_file):
    """
    Reads a file with multiple polygons, formatted as 
    x,y
    x,y
    x,y
    
    x,y
    x,y
    x,y ...
    
    I.e each poly is defined by x,y position of vertices. New polygon starts after
    a space.
    
    Returns a list of polygons 
    """
    delimiter = ','
    fid = open(multi_P_file)
    lines = fid.readlines()
    fid.close()
    polygon = []
    polygons = []
    for line in lines:
        fields = line.split(delimiter)
        try:
            polygon.append([float(fields[0]), float(fields[1])])
        except:
            # Found a line without correct data, assume this signifies the start of a new polygon
            polygons.append(polygon)
            polygon = []
        
    # Pickup the last polygon
    polygons.append(polygon)

    #print len(polygons)

    #print polygons    
    return polygons


# Define a function to read Single File with Multi-polygons and attribute a value
def read_multi_poly_file_value(multi_P_file,attribute):
    """
    Reads a file with multiple polygons, formatted as 
    x,y
    x,y
    x,y
    
    x,y
    x,y
    x,y ...
    
    I.e each poly is defined by x,y position of vertices. New polygon starts after
    a space.
    
    Returns a list of tuples (polygon, attribute) 
    """
    delimiter = ','
    fid = open(multi_P_file)
    lines = fid.readlines()
    fid.close()
    polygon = []
    polygon_value_list = []
    for line in lines:
        fields = line.split(delimiter)
        try:
            polygon.append([float(fields[0]), float(fields[1])])
        except:
            # Found a line without correct data, assume this signifies the start of a new polygon
            pair = [polygon, attribute] # create polygon , value pair....
            polygon_value_list.append(pair) # add it to the list....
            polygon = []
        
    # Pickup the last polygon
    pair = [polygon, attribute]

    #print '================================='
    polygon_value_list.append(pair)

    #print len(polygon_value_list)

    #print polygon_value_list
    return polygon_value_list 


