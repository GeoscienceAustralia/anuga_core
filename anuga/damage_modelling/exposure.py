
import csv

from anuga.anuga_exceptions import TitleValueError, \
                                    DataMissingValuesError

from anuga.file.csv_file import load_csv_as_dict

from anuga.geospatial_data.geospatial_data import Geospatial_data,\
     ensure_absolute

LAT_TITLE = 'LATITUDE'
LONG_TITLE = 'LONGITUDE'
X_TITLE = 'x'
Y_TITLE = 'y'

def cmp_0(a, b):
    return (a > b) - (a < b)

class Exposure(object):
    """Class for National Exposure Database storage (NEXIS).
    Returns a csv file handle
    """

    def __init__(self,file_name, latitude_title=LAT_TITLE,
                 longitude_title=LONG_TITLE, is_x_y_locations=None,
                 x_title=X_TITLE, y_title=Y_TITLE,
                 refine_polygon=None, title_check_list=None):
        """
        This class is for handling the exposure csv file.
        It reads the file in and converts the lats and longs to a geospatial
        data object.
        Use the methods to read and write columns.

        The format of the csv files it reads is;
           The first row is a title row.
           comma's are the delimiters
           each column is a 'set' of data

        Feel free to use/expand it to read other csv files.

        It is not for adding and deleting rows

        Can geospatial handle string attributes? It's not made for them.
        Currently it can't load and save string att's.

        So just use geospatial to hold the x, y and georef? Bad, since
        different att's are in diferent structures.  Not so bad, the info
        to write if the .csv file is saved is in attribute_dic

        The location info is in the geospatial attribute.
        """

        self._file_name = file_name
        self._geospatial = None #

        # self._attribute_dic is a dictionary.
        #The keys are the column titles.
        #The values are lists of column data

        # self._title_index_dic is a dictionary.
        #The keys are the column titles.
        #The values are the index positions of file columns.
        self._attribute_dic, self._title_index_dic = \
            load_csv_as_dict(self._file_name, \
            title_check_list=title_check_list)
        try:
            #Have code here that handles caps or lower
            lats = self._attribute_dic[latitude_title]
            longs = self._attribute_dic[longitude_title]
        except KeyError:
            # maybe a warning..
            #Let's see if this works..
            if False != is_x_y_locations:
                is_x_y_locations = True
            pass
        else:
            self._geospatial = Geospatial_data(latitudes=lats,
                                               longitudes=longs)

        if is_x_y_locations is True:
            if self._geospatial is not None:
                pass #fixme throw an error
            try:
                xs = self._attribute_dic[x_title]
                ys = self._attribute_dic[y_title]
                points = [[float(i),float(j)] for i,j in zip(xs,ys)]
            except KeyError:
                # maybe a warning..
                msg = "Could not find location information."
                raise TitleValueError(msg)
            else:
                self._geospatial = Geospatial_data(data_points=points)

        # create a list of points that are in the refining_polygon
        # described by a list of indexes representing the points

    def __cmp__(self, other):
        """Compare this and another object.

        self   this object
        other  the other object

        Returns True if objects are the 'same'.
        """

        # FIXME: Deprecate this method

        #check that 'other' is an instance of this class
        if isinstance(self, type(other)):
            result = cmp_0(self._attribute_dic, other._attribute_dic)
            if result != 0:
                return result

            # The order of the columns is important. Therefore..
            result = cmp_0(self._title_index_dic, other._title_index_dic)
            if result != 0:
                return result
            for self_ls, other_ls in zip(self._attribute_dic,
                                         other._attribute_dic):
                result = cmp_0(self._attribute_dic[self_ls],
                             other._attribute_dic[other_ls])
                if result != 0:
                    return result
            return 0
        else:
            return 1

    def __eq__(self, other):
        """Compare this and another object.

        self   this object
        other  the other object

        Returns True if objects are the 'same'.
        """

        if isinstance(self, type(other)):
            #print(self._attribute_dic)
            #print(other._attribute_dic)

            # Note (Ole) Dictionaries are now sorted in Python3
            #result = self._attribute_dic == other._attribute_dic
            #if result:
            #    return(result)

            # However, to work also with Python use this
            if self._title_index_dic != other._title_index_dic:
                return False

            for title in self._title_index_dic:
                if self._attribute_dic[title] != other._attribute_dic[title]:
                    return False

            # All matched up
            return True

            #if result != 0:
            #    return result

            # The order of the columns is important. Therefore..
            #result = cmp_0(self._title_index_dic, other._title_index_dic)
            #if result != 0:
            #    return result
            #for self_ls, other_ls in zip(self._attribute_dic,
            #                             other._attribute_dic):
            #    result = cmp_0(self._attribute_dic[self_ls],
            #                 other._attribute_dic[other_ls])
            #    if result != 0:
            #        return result
            #return False
        else:
            return False



    def get_column(self, column_name, use_refind_polygon=False):
        """Get a list of column values given a column name.

        column_name         The name of the column to get values from
        use_refind_polygon  if True, only return values in the refined polygon
                            (not yet implemented)

        Note, the type of the values will be String!
        do this to change a list of strings to a list of floats
        time = [float(x) for x in time]
        """

        if column_name not in self._attribute_dic:
            msg = 'There is no column called %s!' % column_name
            raise TitleValueError(msg)

        return self._attribute_dic[column_name]

    def get_value(self, value_column_name, known_column_name,
                  known_values, use_refind_polygon=False):
        """
        Do linear interpolation on the known_colum, using the known_value,
        to return a value of the column_value_name.
        """

        pass

    def get_location(self, use_refind_polygon=False):
        """
        Return a geospatial object which describes the
        locations of the location file.

        Note, if there is not location info, this returns None.

        Not implemented:
        if use_refind_polygon is True, only return values in the
        refined polygon
        """

        return self._geospatial

    def set_column(self, column_name, column_values, overwrite=False):
        """
        Add a column to the 'end' (with the right most column being the end)
        of the csv file.

        Set overwrite to True if you want to overwrite a column.

        Note, in column_name white space is removed and case is not checked.
        Precondition
        The column_name and column_values cannot have comma's in it.
        """

        # sanity checks
        value_row_count = \
                len(self._attribute_dic[list(self._title_index_dic.keys())[0]])
        if len(column_values) != value_row_count:
            msg = 'The number of column values must equal the number of rows.'
            raise DataMissingValuesError(msg)

        # check new column name isn't already used, and we aren't overwriting
        if column_name in self._attribute_dic:
            if not overwrite:
                msg = 'Column name %s already in use!' % column_name
                raise TitleValueError(msg)
        else:
            # New title.  Add it to the title index.
            self._title_index_dic[column_name] = len(self._title_index_dic)

        self._attribute_dic[column_name] = column_values

    def save(self, file_name=None):
        """Save the exposure csv file.

        file_name  if supplied write this to filename, not original
        """

        if file_name is None:
            file_name = self._file_name

        fd = open(file_name, 'w', newline="")
        writer = csv.writer(fd)

        #Write the title to a cvs file
        line = [None] * len(self._title_index_dic)
        for title in self._title_index_dic.keys():
            line[self._title_index_dic[title]] = title
        writer.writerow(line)

        # Write the values to a cvs file
        value_row_count = \
                len(self._attribute_dic[list(self._title_index_dic.keys())[0]])
        for row_i in range(value_row_count):
            line = [None] * len(self._title_index_dic)
            for title in self._title_index_dic.keys():
                line[self._title_index_dic[title]] = \
                     self._attribute_dic[title][row_i]
            writer.writerow(line)
