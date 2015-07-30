#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Setup base data for ANUGA run

Gareth Davies, Geoscience Australia 2014+

"""
import sys
import glob
import os
import time
import shutil
import hashlib
import json
import numpy
import anuga
from anuga.utilities import spatialInputUtil as su
from anuga.parallel import myid, barrier, send, receive, numprocs

# Local modules
from read_boundary_tags_line_shapefile import \
    read_boundary_tags_line_shapefile
from setup.parse_input_data import ProjectData


class Logger(object):

    """Makes it simple to get stdout to flush, and seems to work in
        parallel
    """

    def __init__(self, logfile):
        self.log = open(logfile, "a")

    def write(self, message):
        self.log.write(message)
        self.log.flush()


class PrepareData(ProjectData):

    """Converts the information in the configuration file to inputs
        suitable for ANUGA

    """

    def __init__(self, filename, make_directories=True, output_log=None):
        """Parse the input data then process it for ANUGA

            @param filename = configuration file (xls)
            @param make_directories Create output directories for simulation           
            @param output_log filename to redirect stdout (inside output
                directories)

        """
        # Get the 'raw' data
        ProjectData.__init__(self, filename)

        # Create a unique output directory, and redirect stdout there
        self.define_output_directory_and_redirect_stdout(
            make_directories=make_directories,
            output_log=output_log)

        # Read files / pre-process
        self.process_project_data()

        if make_directories:
            # Make other required directories, copy code files, etc
            self.setup_directory_structure()

            # Export spatial data to some text
            self.export_spatial_data_as_text()

        barrier()

    def define_output_directory_and_redirect_stdout(self,
                                                    make_directories=True, 
                                                    output_log=None):
        """Make the main output directory, and redirect stdout to a file there

            @param output_log Name of file (stored inside the output directory)
                    which will hold stdout

        """
        if myid == 0:
            runtime = time.strftime('%Y%m%d_%H%M%S', time.localtime())
            for i in range(1, numprocs):
                send(runtime, i)
        else:
            runtime = receive(0)
        barrier()

        self.output_dir = self.output_basedir + 'RUN_' + str(runtime) +\
            '_' + self.scenario

        if ((make_directories) and (myid == 0)):
            try:
                os.mkdir(self.output_basedir)
            except:
                pass

            # Make the output directory

            try:
                os.mkdir(self.output_dir)
            except:
                pass

            print 'OUTPUT_DIRECTORY: ' + str(self.output_dir)

        # Send stdout to a file inside the output directory
        if output_log is not None:
            if make_directories:
                stdout_file = self.output_dir + '/' + output_log
            else:
                stdout_file = output_log

            if myid == 0:
                print 'Redirecting output now to ' + stdout_file
                sys.stdout = Logger(stdout_file)
            barrier()

            if myid != 0:
                sys.stdout = Logger(stdout_file)

        return

    def process_project_data(self):
        """Process the input data ready for ANUGA

        """

        # Print messages from the ProjectData.__init__ call
        # This allows us to log those messages without refactoring
        # (Consider refactoring though)
        if myid == 0:
            for p in self.print_info:
                print p
            print ''
            print '---------------------'
            print 'PROCESS_PROJECT_DATA'
            print '---------------------'
            print ''
            # Record the time and broadcast to other processers
            time_number = time.time()
            if numprocs > 1:
                for i in range(1, numprocs):
                    send(time_number, i)
        else:
            time_number = receive(0)

        # We can either use interior regions, or breaklines

        if not self.interior_regions_data == []:
            assert self.pt_areas is None, \
                'Cannot define both ptAreas and non-empty interior regions'

        bounding_polygon_and_tags = \
            read_boundary_tags_line_shapefile(
                self.bounding_polygon_and_tags_file,
                self.boundary_tags_attribute_name)
        self.bounding_polygon = bounding_polygon_and_tags[0]
        self.boundary_tags = bounding_polygon_and_tags[1]

        self.breaklines = su.readListOfBreakLines(self.breakline_files)
        (self.riverwalls, self.riverwall_par) = \
            su.readListOfRiverWalls(self.riverwall_csv_files)

        if self.pt_areas is not None:
            self.region_point_areas = su.readRegionPtAreas(
                self.pt_areas,
                convert_length_to_area=self.region_resolutions_from_length)
        else:
            self.region_point_areas = None

        # Hack to override resolution
        # region_point_areas=\
        # [ region_point_areas[i][0:2]+[150*150*0.5] for i in \
        #                      range(len(region_point_areas))]

        # Redefine interior_regions to contain the polygon data + resolutions

        self.interior_regions = [[su.read_polygon(ir[0]), ir[1]] for ir in
                                 self.interior_regions_data]

        # Deal with intersections in the bounding polygon / breaklines /
        # riverwalls. At the moment we cannot add points to the bounding
        # polygon because the boundary tags are not adjusted -- so check that
        # the length of the bounding polygon doesn't change
        lbp = len(self.bounding_polygon)
        if type(self.break_line_intersect_point_movement_threshold) is not str:
            (self.bounding_polygon, self.breaklines, self.riverwalls) = \
                su.add_intersections_to_domain_features(
                    self.bounding_polygon,
                    self.breaklines,
                    self.riverwalls,
                    point_movement_threshold=self.break_line_intersect_point_movement_threshold,
                    verbose=True)

        msg = 'Bounding polygon had points added or dropped because of ' + \
              'intersections --' + \
              'This is not yet properly supported.  Please add ' + \
              ' the intersection points to the bounding polygon'
        assert lbp == len(self.bounding_polygon), msg

        # Here we make a unique ID based on the all the mesh geometry inputs
        # This tells us if we need to regenerate partitions, or use old ones
        mesh_dependency_information = [
                self.bounding_polygon,
                self.interior_regions,
                self.riverwalls,
                self.breaklines,
                self.region_point_areas,
                self.default_res,
                self.boundary_tags
            ]

        if not self.use_existing_mesh_pickle:
            # Append the time to the mesh dependency so we don't reuse old
            # meshes
            mesh_dependency_information.append([time_number])

        self.mesh_id_hash = hashlib.md5(json.dumps(mesh_dependency_information)).hexdigest()

        # Fix the output tif bounding polygon
        if self.output_tif_bounding_polygon is None:
            self.output_tif_bounding_polygon = self.bounding_polygon
        else:
            self.output_tif_bounding_polygon = \
                su.read_polygon(self.output_tif_bounding_polygon)

        # Make proj4string from projection information
        #

        if isinstance(self.projection_information, int):

            # projection_information describes a UTM zone
            # e.g. '+units=m +ellps=WGS84 +zone=47 +south=False +proj=utm '

            if self.projection_information < 0:
                self.proj4string = '+proj=utm +zone=' \
                    + str(abs(self.projection_information)) \
                    + ' +south +datum=WGS84 +units=m +no_defs'
            else:
                self.proj4string = '+proj=utm +zone=' \
                    + str(self.projection_information) \
                    + ' +datum=WGS84 +units=m +no_defs'
        elif isinstance(self.projection_information, str):
            self.proj4string = self.projection_information
        else:
            msg = 'Invalid projection information ' + \
                ' --  must be a proj4string, or an integer' + \
                ' defining a UTM zone [positive for northern hemisphere,' + \
                ' negative for southern hemisphere]'
            raise Exception(msg)

        # Set up directories etc

        self.partition_basedir = 'PARTITIONS/'
        self.partition_dir = self.partition_basedir + 'Mesh_' +\
            str(self.mesh_id_hash)
        self.meshname = self.output_dir + '/mesh.tsh'

    def setup_directory_structure(self):
        """Make 'detailed' output directories, send important data to files

        This has no arguments, and returns nothing, but has lots of important
        side effects
        """
        if myid == 0:

            # Make the partitions directory (2 steps as it is recursive)

            try:
                os.mkdir(self.partition_basedir)
            except:
                pass

            try:
                os.mkdir(self.partition_dir)
            except:
                pass

            # Make the spatialTxt directory

            try:
                os.mkdir(self.spatial_text_output_dir)
            except:
                pass

            # Copy code files to output dir. This should be 'nearly the first'
            # thing we do
            code_output_dir = self.output_dir + '/code'
            setup_code_output_dir = code_output_dir + '/setup'
            try:
                os.mkdir(code_output_dir)
                os.mkdir(setup_code_output_dir)
            except:
                pass

            model_files = glob.glob('*.py') + glob.glob('*.sh') + \
                [self.config_filename]
            for model_file in model_files:
                anuga.copy_code_files(code_output_dir, model_file)

            model_files = glob.glob('setup/*.py')
            for model_file in model_files:
                anuga.copy_code_files(setup_code_output_dir, model_file)

        return

    def export_spatial_data_as_text(self):
        """Export the imported spatial data as text
           This is useful to check that the data was imported correctly

        """
        # #################################################################
        # Export breaklines to txt for QC

        def sav_lines(lines, filename='check.txt'):
            """quick convenience function
            """
            new_a = numpy.array(lines)
            numpy.savetxt(filename, new_a)
            return

        if myid == 0:
            if self.interior_regions is not []:
                for i in range(len(self.interior_regions)):
                    sav_lines(self.interior_regions[i][0],
                              filename=self.spatial_text_output_dir
                              + '/interior_region_' + str(i) + '.txt')
            if self.breaklines is not {}:
                for i in range(len(self.breaklines)):
                    sav_lines(self.breaklines.values()[i],
                              filename=self.spatial_text_output_dir +
                              '/breakline_' + str(i) + '.txt')
            if self.riverwalls is not {}:
                for i in range(len(self.riverwalls)):
                    sav_lines(self.riverwalls.values()[i],
                              filename=self.spatial_text_output_dir +
                              '/riverwall_' + str(i) + '.txt')

            sav_lines(self.bounding_polygon,
                      filename=self.spatial_text_output_dir
                      + '/boundingPoly.txt')

            if self.region_point_areas is not None:
                sav_lines(self.region_point_areas,
                          filename=self.spatial_text_output_dir
                          + '/regionPtAreas.txt')

            # Copy all the txt representations of these files to the output
            # directory

            shutil.copytree(self.spatial_text_output_dir,
                            self.output_dir +
                            '/' + self.spatial_text_output_dir)
        return
