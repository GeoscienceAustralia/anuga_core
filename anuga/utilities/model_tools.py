r"""
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
from anuga import Boyd_box_operator
from anuga import Boyd_pipe_operator
from anuga import Weir_orifice_trapezoid_operator
from anuga import Inlet_operator


def get_polygon_from_single_file(Rfile):
    fid = open(Rfile)
    lines = fid.readlines()
    fid.close()
    polylist = []
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
    polylist.append(polygon)
    return polylist


def get_polygons_from_Mid_Mif(Rfile):
    """Create List of Polygons from a Directory with a File containing Multiple-Polygons
       
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    #print 'in get_polygon_from_Mid_Mif -line 126'
    fid = open(Rfile)
    lines = fid.readlines()
    fid.close()
    # Got the Multi Poly file now process
    polylist = []
    polygon = []
    check_pts_list=[]
    Poly_count=0
    Poly_line_count=0
    Points_in_Poly = 100
    total_lines_in_file= len(lines)
    #print "total number of lines in the Polygons FILE is: ",total_lines_in_file

    for i, line in enumerate(lines): 
        if line.strip().startswith('Region'):
            Poly_line_count=0
            check_pts_list=[]
            if Poly_count==0:
                pass
            else:
                polylist.append(polygon)
                #outfid.close()
                #print polygon
                polygon=[]
            Poly_count+=1
            # Create A poly File for each polygon defined in the Multi-Poly file
            #print 'Polygon #',Poly_count
            #poly_write_file="Poly_"+str(Poly_count)+".csv"
            #outfid = open(poly_write_file, 'w')
            #raw_input('Check Output... -line 155')  
            # Instead write to a List
        elif line.strip().startswith('    Pen'):
            pass
        elif line.strip().startswith('    Brush'):
            pass
        else:
            Poly_line_count+=1
            if Poly_line_count > 1 and Poly_line_count <= (Points_in_Poly+1) and Poly_count!=0:
                #print line, #Points_in_Poly,#Poly_line_count
                fields = line.split(' ')
                if line in check_pts_list and Poly_line_count != Points_in_Poly+1:   # Get rid of any doubled up points NOTE this gets rid of last line !!!
                    #print Poly_line_count, Points_in_Poly+1
                    pass
                else:
                    #outfid.write("%.3f,%.3f\n" % (float(fields[0]),float(fields[1])))
                    polygon.append([float(fields[0]),float(fields[1])])
                    check_pts_list.append(line)
            elif Poly_line_count==1 and Poly_count!=0:
                # read number of points in poly
                #print 'line=',line
                Points_in_Poly=int(line)
                #print 'Number Points in Poly =',Points_in_Poly
                #print polygon
                #raw_input('Check Poly...')
    #print 'End For Loop....-line178'
    polylist.append(polygon)
    #print polylist
    #outfid.close()          
    return polylist          


def get_polygon_list_from_files(dir, verbose = False):
    """Read all polygons found in specified dir and return them in a list
    
       Called by:
       get_polygon_dictionary
    
       Purpose:
       To fill a list with all of the polygons read under a specified directory
    
       Calls:
       anuga.utilities.polygon.read_polygon
    """
    
    #print 'Reading polygon files from ' + dir
    #print 'This will check the file for Multiple Polygons or read mutiple files with a single polygon per file...' # Need to read files with multiple polys also....
    polylist = []
    for filename in os.listdir(dir):
        Rfile = dir +'/'+filename
        if verbose: print(Rfile)
        
        if Rfile[-4:] == '.csv':
            polygon = get_polygon_from_single_file(Rfile)
        if Rfile[-4:] == '.mif':
            polygon = get_polygons_from_Mid_Mif(Rfile)
        
        polylist = polylist + polygon

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
       
       Calls:
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


def get_polygon_value_list(dir):
    """Create list of multiple Polygons attributed with a value
       Where the values are obtained from sub directory names based on number and decimal at underscore
       So: Passing Directory ROUGHNESS containing, subs, 0_015, and 0_06 for example
       
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
            print('Non numerical attributes not yet implemented. I got %s' % key)
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    #print polygon_value_list
    return polygon_value_list


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
        raise KeyError('Files not defined in dictionary: %s' % msg[2:])

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
    and populates Operator parameters.    
    
    """
    #print culvert_bridge_file
    globals={}
    locals={}    
    
    exec(open(culvert_bridge_file).read(), globals, locals) 
    #print locals
    #if 'height' and 'z1' and 'z2' in locals:
    if 'z1' and 'z2' in locals:
        culvert = Weir_orifice_trapezoid_operator(domain, **locals)
    elif 'diameter' in locals:
        culvert = Boyd_pipe_operator(domain, **locals)
    elif 'height' in locals:
        culvert = Boyd_box_operator(domain, **locals)
    else:
        raise Exception('Cant create culvert')



#-----------------------------------------------------------------------------------------
#          FUNCTION FOR BLOCKAGE
#-----------------------------------------------------------------------------------------
def get_WCC_2016_Blockage_factor(Structure,Event,Scenario, long_result=False, verbose=True):
    """
    If the Structure has a single dimension it is assumed a Diameter (hence Pipe)
    Else it is a Box with w & h
    The Event is grouped 1,2,5= small, 10,20 = med, 50,100,pmp = large
    The Scenario can be a Design or Risk Management Outlook
    Based on these there are two arrays containng the Blockage Factor to be applied
    
    
    --------------------------------------------------------------------
    2017 - 02 - 03
    
    Wollongong Blockage Factor Calculator
    
    Author: Rudy Van Drie
    
    --------------------------------------------------------------------
    Class P1. 
    Pipes 1.2 m internal diameter or smaller.  
    
    Class P2. 
    Pipes  greater  than  1.2  m  internal  diameter.
    
    Class B1
    Box culverts or bridges with a diagonal opening less than 1.5 m,  
    and a width or height less than 0.9 m. 
    
    Class B2. 
    Box  culverts  or  bridges  with a diagonal opening of more than or equal to 1.5 m, less than 3 m 
    and minimum dimension of 0.9 m for both width and height. 
    >= 0.9m w and h
    Class 3. 
    Box culverts or bridges with a diagonal opening of more than or equal 
    to  3  m,  less  than  6  m,  
    and  a  minimum  dimension  of  1.2  m  for  both  width and height.   
    
    Class 4. 
    Box culverts or bridges with a diagonal opening greater than or equal 
    to 6 m, and a minimum dimension of 2.5 m for both width and height.   
    
        CLASSP1   Diam =< 1.2m
        CLASSP2   Diam > 1.2m
    
        CLASSB1:    diag < 1.5m and W or H < 0.9m 
                    
        CLASSB2:    diag >= 1.5m AND diag < 3.0m AND both W and H >= 0.9m
                                    
        CLASSB3:    diag >= 3.0m AND diag < 6.0m AND both W and H >= 1.2m 
    
        CLASSB4:    diag >= 6.0m AND W and H >= 2.5m 
    
    
    DESIGN BLOCKAGE FACTORS
                      CLP1,CLP2
    event,            CLB1,CLB2,CLB3,CLB4
    1,2,5,small,      0.35,0.25,0.15,0.00 
    10,20,medium,     0.50,0.40,0.30,0.05 
    50,100,pmp,large, 0.70,0.50,0.40,0.10 
    
    RISK MANAGEMENT BLOCKAGE FACTORS
                      CLP1,CLP2
    event,            CLB1,CLB2,CLB3,CLB4
    1,2,5,small,      0.60,0.50,0.35,0.05 
    10,20,medium,     0.75,0.65,0.50,0.10 
    50,100,pmp,large, 0.95,0.75,0.60,0.15    
    
    """
    
    # REQUIRED DATA FOR small,medium,large and class 1,2,3,4 for two Scenarios
    BF_DES = [[0.35,0.25,0.15,0.00],[0.50,0.40,0.30,0.05],[0.70,0.50,0.40,0.10]]
    BF_RMN = [[0.60,0.50,0.35,0.05],[0.75,0.65,0.50,0.10],[0.95,0.75,0.60,0.15]]
    
    
    if len(Structure) > 1:# ====== FOR BOX =================
        h = float(Structure[0])
        w = float(Structure[1])
        diag = (w**2+h**2)**0.5
                
        if diag >= 6.00 and w >= 2.5 and h >= 2.5:
            BF_clss = 'CLASS B4'                                    
            cclass = 3
        elif diag >= 3.0 and w >= 1.2 and h >= 1.2:
            BF_clss = 'CLASS B3'            
            cclass = 2
        elif diag >= 1.5 and w >= 0.9 and h >= 0.9:
            BF_clss = 'CLASS B2'            
            cclass = 1
        elif diag < 1.5 or w < 0.9 or h < 0.9:
            BF_clss = 'CLASS B1'
            cclass = 0
    else:   # ====== FOR PIPE ================
        d = float(Structure[0])
        if d < 1.2:
            diag = d
            BF_clss =  'CLASS P1'
            cclass = 0
        else:
            diag = d
            BF_clss =  'CLASS P2'
            cclass = 1
            
    if Event in [1,2,5]: 
        Ev_Row = 0
        Ev_mag = 'Small'
    elif Event in [10,20]: 
        Ev_Row = 1 
        Ev_mag = 'Medium'
    elif Event in [50,100,9999]: 
        Ev_Row = 2
        Ev_mag = 'Large'
    
    if Scenario == 'D':
        Scenario = 'DESIGN'
        BF = BF_DES[Ev_Row][cclass]
    elif Scenario == 'R':
        Scenario = 'RISKMAN'
        BF = BF_RMN[Ev_Row][cclass]

    if verbose:
        print('       Importing Culverts')
        print('   Culvert Size ([H,W] or [d]): ', Structure)
        print('                    Event Size: ', Ev_mag)
        print('             Blockage Scenario: ', Scenario)
        print('               Blockage Factor: ', BF)
        print('')

    if long_result:
        return(Scenario, Ev_mag,BF_clss,diag,BF)
    else:
        return(BF)



def get_WCC_2002_Blockage_factor(Structure, verbose=True):
    """
    If the Structure has a single dimension it is assumed a Diameter (hence Pipe)
    Else it is a Box with w & h

    --------------------------------------------------------------------
    2017 - 06 - 22
    
    Wollongong Blockage Factor Calculator for 2002 Blockage Policy
    
    Author: Petar Milevski
    
    --------------------------------------------------------------------
    For all design storm events    
     
    if diag >= 6.0m, blockage factor = 0.25, otherwise blockage factor is 100%

    """
    
    if len(Structure) > 1:# ====== FOR BOX =================
        h = float(Structure[0])
        w = float(Structure[1])
        diag = (w**2+h**2)**0.5
                
    else:   # ====== FOR PIPE ================
        d = float(Structure[0])
        diag = d
     
    if diag >= 6.0:
        BF = 0.25
    else:
        BF = 1.0

    if verbose:
        print('       Importing Culverts')
        print('   Culvert Size ([H,W] or [d]): ', Structure)
        print('                      Diagonal: ', diag)
        print('               Blockage Factor: ', BF)
        print('')

    return(BF)
