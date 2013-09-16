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
import numpy
from anuga.geometry.polygon import read_polygon
from anuga.structures.boyd_box_operator import Boyd_box_operator
from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator


def get_polygon_list_from_files(dir):
    """Read all polygons found in specified dir and return them in a list
       Called by:
       get_polygon_dictionary
       Purpose:
       To fill a list with all of the polygons read under a specified directory
       CALLS:
       anuga.utilities.polygon.read_polygon
    """
    
    #print 'Reading polygon files from ' + dir
    #print 'This will check the file for Multiple Polygons or read mutiple files with a single polygon per file...' # Need to read files with multiple polys also....
    polylist = []
    for filename in os.listdir(dir):
        Rfile = dir +'/'+filename
        print Rfile
        #print filename
        
        # Check if file contains blank lines, if so multi poly's
        # Read file, check number of lines, if blank before last = multi 
        fid = open(Rfile)
        lines = fid.readlines()
        fid.close()
        polygon = []
        polycount = 0
        for line in lines:
            fields = line.split(',')
            #print line
            if line in ('\n', '\r\n'): # Must be a blank Line....
                # Found a line without INcorrect data, assume this signifies the start of a new polygon
                polycount+=1
                #print 'Polygon '+str(polycount)
                polylist.append(polygon)
                #print polygon
                polygon =[]
            else:
                polygon.append([float(fields[0]), float(fields[1])])
            
            """
            try:
                polygon.append([float(fields[0]), float(fields[1])])
            except:
                # Found a line without INcorrect data, assume this signifies the start of a new polygon
                polycount+=1
                print 'Polygon '+str(polycount)
                polylist.append(polygon)
                polygon =[]
            """
        polylist.append(polygon)
    #print polylist
    #raw_input('hold at polylist..')
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
        # How to read a file with multiple polygons ??
        D[a] = get_polygon_list_from_files(os.path.join(dir, a)) # Fill Item [a] in the Dictionary with FIle name and attribute
    return D

# ---- GENERIC POLYGON VALUE LIST Generator
def get_polygon_value_list(dir):
    """Create list of multiple Polygons attributed with a value
       Where the values are obtained from directory names 
       that is a List of Polygons attributed with the Value read from the directory name...
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    
    #print 'Read directories of polygons and attributing DIR NAME to polygon'    
    #print 'Naming convention uses the underscore as decimal point eg:0_015, 1000_0'
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            #print 'Polygon Attribute = ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    #print polygon_value_list
    return polygon_value_list


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
    #print 'From DIR: ',dir
    #print 'Holes file_List:-',holes_file_list
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
def read_hole_dir_multi_files_with_single_poly(directory, filepattern='*.csv'):
    """
    Looks in a directory, and reads all .csv files as polygons
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



#define a function with without an attribute
def read_hole_dir_single_file_with_multi_poly(directory, filepattern='*.csv'):
    """
    Looks in a directory, and reads 1 .csv file
    containing muliple polygons
    and returns a list of polygon 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # now get the result list
    result = []
    for f in files: # For the 1 file
        result.append(read_multi_poly_file(multi_P_file)) # Get the multiple Polygons
    return result


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



# Define a function to read Culvert and Bridge data from Files in Directory
def Create_culvert_bridge_Operator(domain,culvert_bridge_file):
    """This script reads in culvert and bridge data files
    and populate Operator parameters.    
    
    """
    #print culvert_bridge_file
    globals={}
    locals={}    
    
    execfile(culvert_bridge_file, globals, locals)
    #print locals
    if 'diameter' in locals:
        culvert = Boyd_pipe_operator(domain, **locals)
    elif 'height' in locals:
        culvert = Boyd_box_operator(domain, **locals)
    else:
        raise Exception, 'Cant create culvert'
    #print culvert.description

"""
# Define a function to read Culvert and Brdige data from Files in Directory
def TESTCreate_culvert_bridge_Operator(domain,culvert_bridge_file):
    # Read file and populate Operator parameters
    print culvert_bridge_file
    delimiter = ','
    fid = open(culvert_bridge_file)
    lines = fid.readlines()
    fid.close()
    line_counter=0
    for line in lines:
        line_counter+=1
        print line
        if line_counter ==1:
            pass
        elif line_counter ==2:
            raw_input('I got here.. 2')
            fields = line.split(delimiter)
            if len(fields) > 4:
                # Two lines
                el0=numpy.array([[fields[0],fields[1]],[fields[2],fields[3]]])
                el1=numpy.array([[fields[4],fields[5]],[fields[6],fields[7]]])
                #el0 = numpy.array([[305945.955,6193836.293] , [305945.125,6193835.387]])
                exchange_lines=[el0,el1]
                Exchange_Type = 'LINES'
                raw_input('I got here.. 3')
            else: # two Points
                #ep0 = numpy.array([296653.0,6180014.9])
                #ep1 = numpy.array([296642.5,6180036.3]) 
                ep0=numpy.array([float(fields[0]),float(fields[1])])
                ep1=numpy.array([float(fields[2]),float(fields[3])])
                exchange_points=[ep0,ep1]  
                Exchange_Type = 'POINTS'
                raw_input('I got here.. 4')
        # Continue if....        
        elif line.strip().startswith('width'):
            delimiter = '='        
            fields = line.split(delimiter)
            width = float(fields[1])
        elif line.strip().startswith('height'):
            delimiter = '='        
            fields = line.split(delimiter)
            height = float(fields[1])
            Culvert_Type = 'BOX'
            raw_input('@ height...')
        elif line.strip().startswith('diameter'):
            delimiter = '='        
            fields = line.split(delimiter)
            diameter = float(fields[1])            
            Culvert_Type = 'PIPE'
            raw_input('@ diam...')
        elif line.strip().startswith('apron'):
            delimiter = '='        
            fields = line.split(delimiter)
            apron = float(fields[1])            

        elif line.strip().startswith('enquiry_gap'):
            delimiter = '='        
            fields = line.split(delimiter)
            enquiry_gap= float(fields[1])            

        elif line.strip().startswith('Mannings'):
            delimiter = '='        
            fields = line.split(delimiter)
            manning = float(fields[1])            
        elif line.strip().startswith('use_momentum_jet'):
            delimiter = '='        
            fields = line.split(delimiter)
            use_momentum_jet = fields[1]            

        elif line.strip().startswith('use_velocity_head'):
            delimiter = '='        
            fields = line.split(delimiter)
            use_velocity_head = fields[1]            

        elif line.strip().startswith('verbose'):
            delimiter = '='        
            fields = line.split(delimiter)
            verbose = fields[1]            
        elif line.strip().startswith('logging'):
            delimiter = '='        
            fields = line.split(delimiter)
            logging = fields[1]            
            
        elif line.strip().startswith('losses'):
            delimiter = '='        
            fields = line.split(delimiter)
            losses = fields[1][1:-1]            
            print losses
        else:
            pass
            #Local_Defaults:
            #losses={'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
       
    print Exchange_Type
    raw_input('hold...')
    # ----- Now Create Operator
    if Culvert_Type =='BOX' and Exchange_Type =='LINES':
        culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=width,
                                    exchange_lines=[el0, el1],
                                    height=height,
                                    apron=apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=use_momentum_jet,
                                    use_velocity_head=use_velocity_head,
                                    manning=manning,
                                    logging=logging,
                                    label=culvert_bridge_file[0:-4],
                                    verbose=verbose)  
    elif Culvert_Type =='BOX' and Exchange_Type == 'POINTS':
            culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=width,
                                    end_points=[ep0, ep1],
                                    height=height,
                                    apron=apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=use_momentum_jet,
                                    use_velocity_head=use_velocity_head,
                                    manning=manning,
                                    logging=logging,
                                    label=culvert_bridge_file[0:-4],
                                    verbose=verbose)  
    elif Culvert_Type =='PIPE' and Exchange_Type == 'LINES':
        culvert = Boyd_pipe_operator(domain,
                                    losses=losses,
                                    exchange_lines=[el0, el1],
                                    diameter=diameter,
                                    apron=apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=use_momentum_jet,
                                    use_velocity_head=use_velocity_head,
                                    manning=manning,
                                    logging=logging,
                                    label=culvert_bridge_file[0:-4],
                                    verbose=verbose)  
    elif Culvert_Type =='PIPE' and Exchange_Type == 'POINTS':
            culvert = Boyd_pipe_operator(domain,
                                    losses=losses,
                                    end_points=[ep0, ep1],
                                    diameter=diameter,
                                    apron=apron,
                                    enquiry_gap=enquiry_gap,
                                    use_momentum_jet=use_momentum_jet,
                                    use_velocity_head=use_velocity_head,
                                    manning=manning,
                                    logging=logging,
                                    label=culvert_bridge_file[0:-4],
                                    verbose=verbose)  
    else:
        pass
"""
