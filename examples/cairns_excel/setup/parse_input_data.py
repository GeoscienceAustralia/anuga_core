#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Read xls(x) data to define ANUGA scenario

Gareth Davies, Geoscience Australia 2014+

"""


import os
import xlrd
import glob
from anuga.utilities import spatialInputUtil as su


class ProjectData(object):

    """Class to hold input data that previously occurred in project.py

    We read the data from an excel file via the AnugaXls class

    """

    def __init__(self, filename):
        """Call the right method to read the filetype

        """

        file_extension = os.path.splitext(filename)[1]

        if file_extension in ['.xls', '.xlsx']:
            self.get_data_from_excel(filename)
        else:
            msg = 'File type not supported'
            raise ValueError(msg)

    def get_data_from_excel(self, filename):
        """Get data from xls or xlsx file

            print statements are directed to self.print_info,
            so we print to a logfile in the output directory.
            (as the name of the output directory depends on parameters read
            from the config file)


        """

        # Send print statements to print_info so we can control when they print
        self.print_info = ['---------------------',
                           'PARSING CONFIG FILE ',
                           '--------------------',
                           '']

        # Parse the xls data

        self.config_filename = filename
        self.data_source = AnugaXls(filename)
        data_source = self.data_source

        # Worksheet names -- define here so they can easily be changed

        project_ws = 'project_settings'
        mesh_ws = 'mesh'
        init_ws = 'initial_conditions'
        init_add_ws = 'initial_condition_additions'
        boundary_ws = 'boundary_conditions'
        inlet_ws = 'inlets'
        rain_ws = 'rainfall'
        bridge_ws = 'bridges'
        pumping_station_ws = 'pumping_stations'

        # #####################################################################
        #
        # Generic 'project' type variables
        #
        # #####################################################################

        self.scenario = data_source.get_var(project_ws, 'scenario', [0, 1],
                                            post_process=str)

        self.output_basedir = data_source.get_var(
            project_ws, 'output_base_directory', [
                0, 1], post_process=str)

        self.yieldstep = data_source.get_var(project_ws, 'yieldstep', [0, 1],
                                             post_process=float)

        self.finaltime = data_source.get_var(project_ws, 'finaltime', [0, 1],
                                             post_process=float)

        self.projection_information = data_source.get_var(
            project_ws, 'projection_information', [
                0, 1])

        # Could be unicode string or float -- if float, we need int

        if isinstance(self.projection_information, float):
            self.projection_information = int(self.projection_information)
        if isinstance(self.projection_information, unicode):
            self.projection_information = str(self.projection_information)

        self.flow_algorithm = data_source.get_var(
            project_ws, 'flow_algorithm', [0, 1], post_process=str)

        # Coerce this to a logical variable

        self.use_local_extrapolation_and_flux_updating = \
            data_source.get_var(project_ws,
                                'use_local_extrapolation_and_flux_updating',
                                offset=[0, 1], post_process=bool)

        self.output_tif_cellsize = data_source.get_var(
            project_ws,
            'output_tif_cellsize',
            offset=[0, 1],
            post_process=float)

        self.output_tif_bounding_polygon = data_source.get_var(
            project_ws,
            'output_tif_bounding_polygon',
            offset=[0, 1],
            post_process=empty_string_to_none)

        self.max_quantity_update_frequency = data_source.get_var(
            project_ws,
            'max_quantity_update_frequency',
            offset=[0, 1],
            post_process=int)

        self.max_quantity_collection_start_time = data_source.get_var(
            project_ws,
            'max_quantity_collection_starttime',
            offset=[0, 1],
            post_process=float)

        if self.max_quantity_collection_start_time >= self.finaltime:
            msg = 'max_quantity_collection_starttime should be < finaltime'
            raise Exception(msg)

        self.store_vertices_uniquely = data_source.get_var(
            project_ws,
            'store_vertices_uniquely',
            offset=[0, 1],
            post_process=bool)

        self.store_elevation_every_timestep = data_source.get_var(
            project_ws,
            'store_elevation_every_timestep',
            offset=[0, 1],
            post_process=bool)

        self.spatial_text_output_dir = data_source.get_var(
            project_ws,
            'spatial_text_output_dir',
            offset=[0, 1],
            post_process=str)

        self.report_mass_conservation_statistics = data_source.get_var(
            project_ws,
            'report_mass_conservation_statistics',
            offset=[0, 1],
            post_process=bool)

        self.report_peak_velocity_statistics = data_source.get_var(
            project_ws,
            'report_peak_velocity_statistics',
            offset=[0, 1],
            post_process=bool)

        self.report_smallest_edge_timestep_statistics = data_source.get_var(
            project_ws,
            'report_smallest_edge_timestep_statistics',
            offset=[0, 1],
            post_process=bool)

        # #####################################################################
        #
        # Mesh type variables
        #
        # #####################################################################

        self.use_existing_mesh_pickle = data_source.get_var(
            mesh_ws, 'use_existing_mesh_pickle', [0, 1], post_process=bool)

        self.bounding_polygon_and_tags_file = data_source.get_var(
            mesh_ws, 'bounding_polygon', [1, 1], post_process=str)

        self.default_res = data_source.get_var(mesh_ws, 'bounding_polygon',
                                               [2, 1], post_process=str)

        self.interior_regions_data = data_source.get_paired_list(
            mesh_ws, 'interior_regions', [1, 1], post_process=string_or_float)

        # breaklines, riverwalls, point movement threshold
        # We allow wildcard type names as input files

        breakline_files = data_source.get_list(
            mesh_ws, 'breakline_files', [0, 1], post_process=str)
        # Allow wildcard names as input for breaklines
        self.breakline_files = glob_files(breakline_files)

        riverwall_csv_files = data_source.get_list(
            mesh_ws, 'riverwall_csv_files', [0, 1], post_process=str)
        # Allow wildcards for riverwalls
        self.riverwall_csv_files = glob_files(riverwall_csv_files)

        self.break_line_intersect_point_movement_threshold = \
            data_source.get_var(mesh_ws, 'breakline_intersection_threshold',
                                [0, 1], post_process=string_or_float)

        self.pt_areas = data_source.get_var(
            mesh_ws, 'region_areas_file', [0, 1], empty_string_to_none)

        # Determine whether the attribute of the pt_areas file defines
        # the maximum mesh triangle area, or an approximate maximum length

        length_or_area = data_source.get_var(
            mesh_ws, 'region_areas_file', [1, 1], str)
        if self.pt_areas is not None:
            if length_or_area == 'length':
                self.region_resolutions_from_length = True
            elif length_or_area == 'area':
                self.region_resolutions_from_length = False
            else:
                msg = 'Invalid value defining attribute from ' + \
                      'region_areas_file. Should be *length* or *area*' + \
                      '. I got *' + length_or_area + '*'
                raise ValueError(msg)

        # Check we are only using interior_region OR breakline methods

        use_ir = self.interior_regions_data != []
        use_bl = self.breakline_files != [] or self.riverwall_csv_files != [] \
            or self.pt_areas is not None

        if use_ir and use_bl:
            msg = ' In worksheet ' + mesh_ws \
                + 'Cannot have both *interior_regions_data* and *breaklines ' \
                + '/ riverwalls / region point areas* being non-empty.'
            raise ValueError(msg)

        # #####################################################################
        #
        # Initial conditions type variables
        #
        # #####################################################################

        self.elevation_data, self.elevation_clip_range, self.elevation_mean = \
            get_initial_condition_data(data_source, init_ws, 'Elevation',
                                       self.print_info)

        self.friction_data, self.friction_clip_range, self.friction_mean = \
            get_initial_condition_data(data_source, init_ws, 'Friction',
                                       self.print_info)

        self.stage_data, self.stage_clip_range, self.stage_mean = \
            get_initial_condition_data(data_source, init_ws, 'Stage',
                                       self.print_info)

        self.xmomentum_data, self.xmomentum_clip_range, self.xmomentum_mean = \
            get_initial_condition_data(
                data_source, init_ws, 'X-velocity_times_depth', self.print_info)

        self.ymomentum_data, self.ymomentum_clip_range, self.ymomentum_mean = \
            get_initial_condition_data(
                data_source, init_ws, 'Y-velocity_times_depth', self.print_info)

        # Get the 'additions' to the initial conditions
        self.elevation_additions, _, _ = \
            get_initial_condition_data(data_source, init_add_ws, 'Elevation',
                                       self.print_info)

        self.friction_additions, _, _ = \
            get_initial_condition_data(data_source, init_add_ws, 'Friction',
                                       self.print_info)

        self.stage_additions, _, _ = \
            get_initial_condition_data(data_source, init_add_ws, 'Stage',
                                       self.print_info)

        self.xmomentum_additions, _, _ = get_initial_condition_data(
            data_source, init_add_ws, 'X-velocity_times_depth',
            self.print_info)

        self.ymomentum_additions, _, _ = get_initial_condition_data(
            data_source, init_add_ws, 'Y-velocity_times_depth',
            self.print_info)

        # #####################################################################
        #
        # Boundary condition type variables
        #
        # #####################################################################

        # This can fail if we don't coerce to str

        self.boundary_tags_attribute_name = data_source.get_var(
            boundary_ws, 'boundary_tags_attribute_name', [0, 1], str)

        self.boundary_data = data_source.get_subtable(
            boundary_ws, 'Boundary_Condition_Tag', [1, 0],
            post_process=string_or_float)

        # #####################################################################
        #
        # Inlet type variables
        #
        # #####################################################################

        self.inlet_data = data_source.get_subtable(
            inlet_ws, 'Inlet_name', [1, 0], post_process=string_or_float)

        # #####################################################################
        #
        # Rainfall type variables
        #
        # #####################################################################

        self.rain_data = data_source.get_subtable(
            rain_ws, 'Time-series_file', [1, 0], post_process=string_or_float)

        # #####################################################################

        # #####################################################################
        #
        # Bridges
        #
        # #####################################################################

        self.bridge_data = []
        bridges_supported = True
        if bridges_supported:
            # Append bridge data to self.bridge_data
            bridge_data = data_source.get_subtable(
                bridge_ws, 'bridges', [1, 1], post_process=string_or_float)
            # Remove rows which are switched off
            for i in range(len(bridge_data)):
                # Only include the bridge data if first column != 1
                if bridge_data[i][0] != 1:
                    new_bridge = bridge_data[i][1:len(bridge_data[i])]
                    self.bridge_data.append(new_bridge)

            for i in range(len(self.bridge_data)):
                bd = self.bridge_data[i]
                # Check that bridge deck + exchange line files exist
                for j in [1, 3, 4, 6]:
                    msg = 'Cannot find file ' + bd[j]
                    assert os.path.exists(bd[j]), msg
                # Add bridge deck to breaklines
                self.breakline_files.extend([bd[1]])
                # Add bridge deck to elevation data
                self.elevation_data = [[bd[1], bd[2]]] + self.elevation_data
                self.elevation_clip_range = [ [-1e+100, 1e+100] ] + \
                    self.elevation_clip_range
        else:
            raise Exception('bridges not supported')

        # ####################################################################
        #
        # Pumping Stations
        #
        # ####################################################################
        pumping_station_data = data_source.get_subtable(
            pumping_station_ws, 'pumping stations', [1, 1], 
            post_process=string_or_float)

        # Remove data with 1 in the 'switch off' column
        self.pumping_station_data = []
        for i in range(len(pumping_station_data)):
            ps = pumping_station_data[i]
            if (ps[0] != 1):
                self.pumping_station_data.append(ps[1:len(ps)])
        

        # Put pumping station basins in the mesh/elevation data
        # and perform checks on files
        for i in range(len(self.pumping_station_data)):
            ps = self.pumping_station_data[i]
            self.breakline_files.extend([ps[6]])
            self.elevation_data = [[ ps[6], ps[7] ]] + self.elevation_data
            self.elevation_clip_range = [ [-1e+100, 1e+100] ] + \
                    self.elevation_clip_range

            for j in [6, 8, 9]:
                msg = 'Cannot find file ' + ps[j]
                assert os.path.exists(ps[j]), msg
               
##############################################################################
#
# END OF CLASS
#
##############################################################################


class AnugaXls(object):

    """Read an xls or xlsx file with the ANUGA input data

    This holds a list of worksheet_names, and
    a dict containing the worksheet_data as lists

    """

    def __init__(self, filename):

        workbook = xlrd.open_workbook(filename)

        self.worksheet_names = workbook.sheet_names()

        self.worksheet_data = {}

        for worksheet_name in self.worksheet_names:
            worksheet = workbook.sheet_by_name(worksheet_name)
            self.worksheet_data[worksheet_name] = \
                get_worksheet_as_list(worksheet)

        return

    def get_var(
        self,
        worksheet_name,
        flag,
        offset=[0, 0],
        post_process=None,
    ):
        """Find the cell in the first column of worksheet_name that (partially)
        matches 'flag', and return the nearby cell value which is 'offset'
        from this by [nrow, ncol]

        Optionally, post_process is a function applied to the cell value
        before returning. This can be useful to e.g. coerce integers to
        bool

        """

        sheet_data = self.worksheet_data[worksheet_name]

        matching_row = find_matching_row(sheet_data, flag, worksheet_name)

        # Check that offset[0] would not exceed the number of rows

        if matching_row + offset[0] >= len(sheet_data):
            return ''

        # If the col offset exceeds the data length, return a blank

        if offset[1] < len(sheet_data[matching_row]):
            output = sheet_data[matching_row + offset[0]][offset[1]]
            if post_process is not None:
                output = post_process(output)
            return output
        else:
            return ''

        return

    def get_paired_list(
        self,
        worksheet_name,
        flag,
        offset=[0, 0],
        post_process=None,
    ):
        """Find the cell in the first column of worksheet_name that (partially)
        matches 'flag', and read data from 2 rows of the spreadsheet,
        starting at a location offset from the matching cell by
        offset = [nrow, ncol].

        The data is returned as a paired list:
        [ [r00, r10], [r01, r11], [r02, r12]]

        """

        paired_list = []

        counter = -1
        while True:
            counter = counter + 1
            cell_x = offset[0]
            cell_y = counter + offset[1]
            next_row1 = self.get_var(worksheet_name, flag, [cell_x, cell_y],
                                     post_process)
            if next_row1 == '':

                # We are at the end of the data

                break
            else:
                cell_x = offset[0] + 1
                next_row2 = self.get_var(
                    worksheet_name, flag, [
                        cell_x, cell_y], post_process)

            paired_list.append([next_row1, next_row2])

        return paired_list

    def get_list(
        self,
        worksheet_name,
        flag,
        offset=[0, 0],
        post_process=None,
    ):
        """Find the cell in the first column of worksheet_name that (partially)
        matches 'flag', and read data from 1 column of the spreadsheet,
        starting at a location offset from the matching cell by
        offset = [nrow, ncol].

        The data is returned as a list:
        [ r00, r01, r02, r03, ...]

        """

        out_list = []

        counter = -1
        while True:
            counter = counter + 1
            cell_x = offset[0]
            cell_y = counter + offset[1]
            next_row1 = self.get_var(worksheet_name, flag, [cell_x, cell_y],
                                     post_process)
            if next_row1 == '':

                # We are at the end of the data

                break
            out_list.append(next_row1)

        return out_list

    def get_subtable(
        self,
        worksheet_name,
        flag,
        offset=[0, 0],
        post_process=None,
    ):
        """Find the cell in the first column of worksheet_name that (partially)
        matches 'flag', and read data from a nearby 'sub-table' in the
        spreadsheet as a list of lists.

        The 'sub-table' starts at a location
        offset fom the matching cell by offset = [nrow, ncol]

        Rows can have different numbers of non-empty columns, but each row
        ends when its first empty cell is found
        """

        out_list = []

        counter = -1
        while True:
            counter = counter + 1
            cell_x = offset[0] + counter
            cell_y = offset[1]
            next_row = self.get_list(worksheet_name, flag, [cell_x, cell_y],
                                     post_process)
            if next_row == []:

                break
            out_list.append(next_row)

        return out_list

    def get_fixed_size_subtable_by_columns(
        self,
        worksheet_name,
        flag,
        dimensions,
        offset=[0, 0],
        post_process=None,
    ):
        """Find the cell in the first column of worksheet_name that (partially)
        matches 'flag', and read data from a nearby 'sub-table' in the
        spreadsheet, which has dimensions=[nrow,ncol] (given by 'dimensions')

        """
        sheet_data = self.worksheet_data[worksheet_name]

        matching_row = find_matching_row(sheet_data, flag, worksheet_name)

        output_data = []

        # Loop column by column and collect the data
        for col in range(dimensions[1]):
            output_data.append([])
            for row in range(dimensions[0]):
                cell_x = offset[0] + matching_row + row
                cell_y = offset[1] + col
                new_val = sheet_data[cell_x][cell_y]
                if post_process is not None:
                    new_val = post_process(new_val)
                output_data[col].append(new_val)

        return output_data
##############################################################################
#
# END OF CLASS
#
##############################################################################

# Functions below here


def find_matching_row(sheet_data, flag, worksheet_name = ""):
    """Find row index of cell in the first column of sheet_data which
       matches 'flag'

    """
    matching_row = None

    # Loop through the rows and look for a match in the first column

    for i in range(len(sheet_data)):
        if flag in sheet_data[i][0]:
            matching_row = i
            break

    if matching_row is None:
        msg = 'Could not find a cell in the first column of worksheet *' \
            + worksheet_name + '* matching with *' + flag + '*'
        raise Exception(msg)

    return matching_row


def get_worksheet_as_list(worksheet):
    """Read a particular xls (xlsx) file worksheet as a list of lists

    [
    [cell 0 0, cell 0 1, cell 0 2, ...],
    [cell 1 0, cell 1 1, cell 1 2, ...],
    ...
    ]

    """

    worksheet_data = []

    for i in range(worksheet.nrows):
        worksheet_data.append([])
        for j in range(worksheet.ncols):
            worksheet_data[i].extend([worksheet.cell_value(i, j)])

    return worksheet_data


def empty_string_to_none(string):
    """Given a string, return None if string == "", and the string value
    (with type str) otherwise
    """

    if string == '':
        return None
    else:
        return str(string)


def string_or_float(value):
    """If value is a number, return it with type float
    If value is anything else, return value of type str
    """

    if isinstance(value, (int, long, float, complex)):
        return float(value)
    else:
        return str(value)


def reformat_clip_range(clip_range):
    """Given the clip_range defined in excel, reformat it to
        go into ANUGA. This includes treating strings
        appropriately using very large hard coded numbers
    """
    assert len(clip_range[0]) == 2

    l = len(clip_range)

    output = list()

    for i in range(l):
        min_clip = clip_range[i][0]
        if isinstance(min_clip, str):
            min_clip = -1.0e+100
        max_clip = clip_range[i][1]
        if isinstance(max_clip, str):
            max_clip = 1.0e+100
        output.append([min_clip, max_clip])

    return output


def glob_files(globname):
    """Given a list of strings to glob certain files,
        return a list of all the files

        We check that some files are alway matched to reduce the risk of input
        errors
    """
    file_list = []
    if len(globname) > 0:
        for i in range(len(globname)):
            matching_files = glob.glob(globname[i])
            if len(matching_files) == 0:
                msg =  'Could not match this file: ', globname[i]
                raise Exception(msg)
            file_list = file_list + \
                glob.glob(globname[i])

    return file_list


def get_initial_condition_data(data_source, worksheet, flag, print_info):
    """Convenience function to extract the initial condition data
       (and initial_condition_additions) from the xls worksheet

       The needs have become more elaborate over time, e.g.
        to support combining 2 line files into a polygon

       Given a character string referring to a quantity which has initial
       conditions in the xls worksheet (e.g. 'Elevation'),
       extract the associated data from the
       'data_source' (an AnugaXls object) on worksheet 'worksheet'

       This assumes a particular format in the excel sheet
    """

    # Read the polygon / value pairs
    quantity_data = data_source.get_paired_list(
        worksheet, flag, [1, 1], post_process=string_or_float)

    # If the polygon is a wildcard, assume it matches 2 lines, read them in,
    # and join them to make a polygon. This is a convenient shorthand
    # for when we have lkl
    for i in range(len(quantity_data)):
        polygon_files = glob.glob(quantity_data[i][0])

        # Check it only matches 0 or 1 or 2 files
        msg = 'Polygon:' + str(i) + ' : ' + quantity_data[i][0] + \
              '  for ' + flag + ' on  worksheet' + \
              worksheet + ' matches > 2 files. We can join at most 2 lines' + \
              'to make a polygon'
        assert len(polygon_files) <= 2, msg

        if len(polygon_files) == 0:
            # Check it is valid
            msg = 'Polygon:' + str(i) + ' : ' + quantity_data[i][0] + \
                  '  for ' + flag + ' on  worksheet' + \
                  worksheet + ' matches no files, and is not All or None ' + \
                  'or Extent (for a raster)'
            assert ((quantity_data[i][0] == 'All') |
                    (quantity_data[i][0] is None) |
                    (quantity_data[i][0] == 'Extent')), msg
        elif len(polygon_files) == 2:
            # If it matches 2, try to combine to 1.
            # This is often required to use pairs of breaklines as polygons
            # Do this by:
            # 1) Setting up the 2 lines as though they were in a
            #    breakline object
            # 2) Using su.polygon_from_matching_breakLines
            print_info.append('Initial ' + flag)
            print_info.append('Combining these files to a polygon: ')
            print_info.append(str(polygon_files))
            print_info.append('')

            l0 = su.read_polygon(polygon_files[0])
            l1 = su.read_polygon(polygon_files[1])
            fake_breakline = {polygon_files[0]: l0, polygon_files[1]: l1}
            fake_match = quantity_data[i][0].split('*')[0]
            out_poly = su.polygon_from_matching_breaklines(
                fake_match, fake_breakline)
            quantity_data[i][0] = out_poly

    # Get the clip_range for each polygon / function pair
    quantity_clip_range = data_source.get_fixed_size_subtable_by_columns(
        worksheet, flag,
        dimensions=[2, len(quantity_data)],
        offset=[3, 1], post_process=string_or_float)

    new_quantity_clip_range = reformat_clip_range(quantity_clip_range)

    # Get sub-grid size for spatial averaging, if applicable
    spatial_average = data_source.get_var(
        worksheet, flag, [5, 1], post_process=string_or_float)

    if type(spatial_average) == str:
        spatial_average = None

    return quantity_data, new_quantity_clip_range, spatial_average
