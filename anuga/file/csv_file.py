"""
    A set of functions which extend the capabilities of the Python csv
    module.
    
    CSV files have the extension .csv, which stands for Comma Separated Value
    file. There is no standardised form for this format, so the user is provided
    with a variety of options for parsing different styles of csv files.
    
    These have been left as functions to aviod confusion with the standard
    csv module.
"""


import csv
import numpy as num
import anuga.utilities.log as log


def load_csv_as_dict(file_name, title_check_list=None, delimiter=',',
                        d_type = str):
    """
    Load in the csv as a dictionary, title as key and column info as value.
    Also, create a dictionary, title as key and column index as value,
    to keep track of the column order.

    file_name The path to the file to read.
    
    title_check_list List of titles that *must* be columns in the file.

    delimiter is the delimiter used to separate the fields
    
    format is one of float, str, int

    return 2 dictionaries: ({key:column}, {title:index}).

    WARNING: Values are returned as strings.
             Do this to change a list of strings to a list of floats
                 time = [float(x) for x in time]
    """

    # FIXME(Ole): Consider dealing with files without headers
    
    attribute_dic = {}
    title_index_dic = {}
    titles_stripped = [] # List of titles

    fid = open(file_name)
    reader = csv.reader(fid, delimiter=delimiter)

    # Read in and manipulate the title info
    titles = next(reader)
    for i, title in enumerate(titles):
        header = title.strip()
        titles_stripped.append(header)
        title_index_dic[header] = i
    title_count = len(titles_stripped)

    # Check required columns
    if title_check_list is not None:
        for title_check in title_check_list:
            if title_check not in title_index_dic:
                msg = 'Reading error. This row is not present %s' % title_check
                raise IOError(msg)


    # Create a dictionary of column values, indexed by column title
    for line in reader:
        n = len(line) # Number of entries
        if n < title_count:
            msg = 'Entry in file %s had %d columns ' % (file_name, n)
            msg += 'although there were %d headers' % title_count
            raise IOError(msg)
        for i, val in enumerate(line[:title_count]):  # skip trailing data
            attribute_dic.setdefault(titles_stripped[i], []).append(d_type(val))

    return attribute_dic, title_index_dic


          
def load_csv_as_array(file_name, delimiter = ','):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a dictionary of numeric arrays.

    file_name The path to the file to read.
    delimiter is the delimiter used to separate the fields    

    See underlying function load_csv_as_dict for more details.
    """

    X, _ = load_csv_as_dict(file_name, delimiter=delimiter)


    # Return result as a dict of arrays
    ret = {}
    for key in list(X.keys()):
        ret[key] = num.array([float(x) for x in X[key]])
            
    return ret


def load_csv_as_matrix(file_name, delimiter = ','):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a numeric matrix.

    file_name The path to the file to read.
    delimiter is the delimiter used to separate the fields    

    See underlying function load_csv_as_dict for more details.
    """

    X, title_indices = load_csv_as_dict(file_name, delimiter=delimiter)

    col_titles = list(title_indices.keys())

    # Return result as a 2D array
    ret = num.zeros((len(X[col_titles[0]]), len(title_indices)), float)

    header = []
    for col_title in col_titles:
        index = title_indices[col_title]
        header.append(col_title)
        for i, x in enumerate(X[col_title]):
            ret[i, index] = float(x)

    return header, ret



def store_parameters(verbose=False, **kwargs):
    """
    Store "kwargs" into a temp csv file, if "completed" is in kwargs,
    csv file is kwargs[file_name] else it is kwargs[output_dir]+details_temp.csv

    Must have a file_name keyword arg, this is what is writing to.
    might be a better way to do this using CSV module Writer and writeDict.

    writes file to "output_dir" unless "completed" is in kwargs, then
    it writes to "file_name" kwargs
    """

    # Check that kwargs is a dictionary
    if not isinstance(kwargs, dict):
        raise TypeError

    # is 'completed' in kwargs?
    completed = 'completed' in kwargs

    # get file name and removes from dict and assert that a file_name exists
    if completed:
        try:
            file_name = str(kwargs['file_name'])
        except:
            raise Exception('kwargs must have file_name')
    else:
        # write temp file in output directory
        try:
            file_name = str(kwargs['output_dir']) + 'detail_temp.csv'
        except:
            raise Exception('kwargs must have output_dir')

    # extracts the header info and the new line info
    line = ''
    header = ''
    count = 0
    keys = list(kwargs.keys())
    keys.sort()

    # used the sorted keys to create the header and line data
    for k in keys:
        header += str(k)
        line += str(kwargs[k])
        count += 1
        if count < len(kwargs):
            header += ','
            line += ','
    header += '\n'
    line += '\n'

    # checks the header info, if the same, then write, if not create a new file
    # try to open!
    try:
        fid = open(file_name, 'r')
        file_header = fid.readline()
        fid.close()
        if verbose: log.critical('read file header %s' % file_header)
    except Exception:
        msg = 'try to create new file: %s' % file_name
        if verbose:
            log.critical(msg)
        #tries to open file, maybe directory is bad
        try:
            fid = open(file_name, 'w')
            fid.write(header)
            fid.close()
            file_header=header
        except:
            msg = 'cannot create new file: %s' % file
            raise Exception(msg)

    # if header is same or this is a new file
    if file_header == str(header):
        fid = open(file_name, 'a')
        fid.write(line)
        fid.close()
    else:
        # backup plan,
        # if header is different and has completed will append info to
        # end of details_temp.cvs file in output directory
        file_name = str(kwargs['output_dir']) + 'detail_temp.csv'
        fid = open(file_name, 'a')
        fid.write(header)
        fid.write(line)
        fid.close()

        if verbose:
            log.critical('file %s', file_header.strip('\n'))
            log.critical('head %s', header.strip('\n'))
        if file_header.strip('\n') == str(header):
            log.critical('they equal')

        msg = 'WARNING: File header does not match input info, ' \
              'the input variables have changed, suggest you change file name'
        log.critical(msg)



def load_csv_as_building_polygons(file_name,
                          floor_height=3):
    """
    Convert CSV files of the form:

    easting,northing,id,floors
    422664.22,870785.46,2,0
    422672.48,870780.14,2,0
    422668.17,870772.62,2,0
    422660.35,870777.17,2,0
    422664.22,870785.46,2,0
    422661.30,871215.06,3,1
    422667.50,871215.70,3,1
    422668.30,871204.86,3,1
    422662.21,871204.33,3,1
    422661.30,871215.06,3,1

    to a dictionary of polygons with id as key.
    The associated number of floors are converted to m above MSL and 
    returned as a separate dictionary also keyed by id.
    
    Optional parameter floor_height is the height of each building story.
    Optional parameter clipping_olygons is a list of polygons selecting
    buildings. Any building not in these polygons will be omitted.
    
    See csv2polygons for more details
    """

    polygons, values = load_csv_as_polygons(file_name,
                                    value_name='floors',
                                    clipping_polygons=None)    

    
    heights = {}
    for key in list(values.keys()):
        v = float(values[key])
        heights[key] = v*floor_height
        
    return polygons, heights                
            

def load_csv_as_polygons(file_name,
                 value_name='value',
                 clipping_polygons=None):
    """
    Convert CSV files of the form:

    easting,northing,id,value
    422664.22,870785.46,2,0
    422672.48,870780.14,2,0
    422668.17,870772.62,2,0
    422660.35,870777.17,2,0
    422664.22,870785.46,2,0
    422661.30,871215.06,3,1
    422667.50,871215.70,3,1
    422668.30,871204.86,3,1
    422662.21,871204.33,3,1
    422661.30,871215.06,3,1

    to a dictionary of polygons with id as key.
    The associated values are returned as a separate dictionary also keyed by id.


    easting: x coordinate relative to zone implied by the model
    northing: y coordinate relative to zone implied by the model    
    id: tag for polygon comprising points with this tag
    value: numeral associated with each polygon. These must be the same for all points in each polygon.
   
    The last header, value, can take on other names such as roughness, floors, etc - or it can be omitted 
    in which case the returned values will be None
    
    Eastings and Northings will be returned as floating point values while
    id and values will be returned as strings.

    Optional argument: clipping_polygons will select only those polygons that are
    fully within one or more of the clipping_polygons. In other words any polygon from
    the csv file which has at least one point not inside one of the clipping polygons
    will be excluded 
    
    See underlying function load_csv_as_dict for more details.
    """

    X, _ = load_csv_as_dict(file_name)

    msg = 'Polygon csv file must have 3 or 4 columns'
    assert len(list(X.keys())) in [3, 4], msg
    
    msg = 'Did not find expected column header: easting'
    assert 'easting' in list(X.keys()), msg
    
    msg = 'Did not find expected column header: northing'    
    assert 'northing' in list(X.keys()), msg
    
    msg = 'Did not find expected column header: northing'        
    assert 'id' in list(X.keys()), msg
    
    if value_name is not None:
        msg = 'Did not find expected column header: %s' % value_name        
        assert value_name in list(X.keys()), msg    
    
    polygons = {}
    if len(list(X.keys())) == 4:
        values = {}
    else:
        values = None

    # Loop through entries and compose polygons
    excluded_polygons={}
    past_ids = {}
    last_id = None
    for i, poly_id in enumerate(X['id']):

        # Check for duplicate polygons
        if poly_id in past_ids:
            msg = 'Polygon %s was duplicated in line %d' % (id, i)
            raise Exception(msg)
        
        if poly_id not in polygons:
            # Start new polygon
            polygons[poly_id] = []
            if values is not None:
                values[poly_id] = X[value_name][i]

            # Keep track of previous polygon ids
            if last_id is not None:
                past_ids[last_id] = i
            
        # Append this point to current polygon
        point = [float(X['easting'][i]), float(X['northing'][i])]

        if clipping_polygons is not None:
            exclude=True
            for clipping_polygon in clipping_polygons:
                if inside_polygon(point, clipping_polygon):
                    exclude=False
                    break
                
            if exclude is True:
                excluded_polygons[poly_id]=True

        polygons[poly_id].append(point)    
            
        # Check that value is the same across each polygon
        msg = 'Values must be the same across each polygon.'
        msg += 'I got %s in line %d but it should have been %s' % \
                            (X[value_name][i], i, values[poly_id])
        assert values[poly_id] == X[value_name][i], msg

        last_id = poly_id

    # Weed out polygons that were not wholly inside clipping polygons
    for poly_id in excluded_polygons:
        del polygons[poly_id]
        
    return polygons, values


            


