"""
Parse the boundary tags line shapefile

Gareth Davies, Geoscience Australia 2014+
"""

import os
import ogr
import numpy
import subprocess

############################################################
def check_output(list_of_commands):
    """
        Adding check_output for python < 2.7
    """
    process = subprocess.Popen(list_of_commands, stdout=subprocess.PIPE)

    output = process.communicate()[0]
    print output

    #output, unused_err = process.communicate()
    #retcode = process.poll()
    #if retcode:
    #    cmd = kwargs.get("args")
    #    if cmd is None:
    #        cmd = popenargs[0]
    #    error = subprocess.CalledProcessError(retcode, cmd)
    #    error.output = output
    #    raise error
    return output

##############################################################
def parse_ogr_info_text(ogr_info, tag_attribute):
    """

    """

    boundary_line_and_tags = []

    for i in range(len(ogr_info)):
        s = ogr_info[i]

        # Find feature definitions, and extract the linestring and tag from
        # each
        if s.startswith('OGRFeature'):
            # Get boundary tag. This will start with two spaces, followed by
            # the tag attribute, followed by ' ('
            feature_tag = None
            counter = 0
            tag_match = '  ' + tag_attribute + ' ('

            while feature_tag is None:
                counter = counter + 1

                if (i + counter) == len(ogr_info):
                    print ogr_info
                    msg = 'Failed to parse the above output from ogr_info' + \
                          '\n Check that your boundary tag attribute name' + \
                          ' is correctly specified in the input file'
                    raise Exception(msg)

                if ogr_info[i + counter].startswith(tag_match):
                    feature_tag = ogr_info[i + counter].split(' = ')[1]
       
            # We now have feature_tag

            # Get the coordinates. They will start with '  LINESTRING ('
            line_coordinates = None
            counter = 0
            line_coordinates_match = '  LINESTRING ('

            while line_coordinates is None:
                counter = counter + 1

                if (i + counter) == len(ogr_info):
                    print ogr_info
                    msg = 'Failed to parse the above output from ogr_info' + \
                          '\n Could not find the linestring of every ' + \
                          'Feature. \n Check that all geometries are part ' + \
                          'of the boundary and have a non-empty LINESTRING'
                    raise Exception(msg)

                if ogr_info[i + counter].startswith(line_coordinates_match):
                    # Hack the coordinates out of the string, as a 2 column
                    # list
                    line_coordinates = ogr_info[i + counter].split('(')[1]
                    line_coordinates = line_coordinates.split(')')[0]
                    line_coordinates = line_coordinates.split(',')
                    line_coordinates = \
                        [ [float(x.split(' ')[0]), \
                           float(x.split(' ')[1])] for x in line_coordinates ]

            # We now have the line_coordinates
            boundary_line_and_tags.append([line_coordinates, feature_tag])

    return boundary_line_and_tags



def get_boundary_tags_from_ogrinfo(shapefile_name,
                                   tag_attribute='bndryTag'):
    """

    Get the boundary tags geometry information out of ogrinfo

    This is a 'fall-back' routine to be used if the user has
    an incompatible python gdal interface

    """

    # Hack work around for python < 2.7
    try:
        ogr_info = subprocess.check_output('ogrinfo -al ' + shapefile_name)
    except:
        ogr_info = check_output(['ogrinfo', '-al', shapefile_name])

    ogr_info = ogr_info.split(os.linesep)

    boundary_line_and_tags = parse_ogr_info_text(ogr_info, tag_attribute)


    return boundary_line_and_tags


def read_boundary_tags_line_shapefile(shapefile_name,
                                      tag_attribute='bndryTag'):
    """
    Read in the boundary lines + tags from an appropriately structured line
    shapefile. (Update: It can now alternatively be a csv file)

    Return the bounding polygon + dict with tags + indices
    """

    if not os.path.exists(shapefile_name):
        msg = 'Cannot find file ' + shapefile_name
        raise ValueError(msg)

    # Step 1: Read the data
    try:
        # Read from a vector GIS file using gdal python interface
        data_source = ogr.Open(shapefile_name)

        assert data_source.GetLayerCount() == 1, \
            ' Bounding polygon + boundary tags file can only have 1 layer'

        layer = data_source.GetLayer()

        boundary_line_and_tags = []

        for feature in layer:
            line = feature.GetGeometryRef().GetPoints()
            line = [list(l) for l in line]
            tag = feature.GetField(tag_attribute)
            boundary_line_and_tags.append([line, tag])
    except:
        # Try to use ogrinfo to get the answer

        boundary_line_and_tags = get_boundary_tags_from_ogrinfo(shapefile_name,
            tag_attribute)


    # Step 2: Convert to bounding polygon + boundary tags in ANUGA format

    boundary_segnum = len(boundary_line_and_tags)

    # Initial values

    bounding_polygon = boundary_line_and_tags[0][0]
    boundary_tags = \
        {boundary_line_and_tags[0][1]: range(len(bounding_polygon) - 1)}

    # Treat the case of a 'closed' polygon with first point == last point
    # by dropping the last point
    if ((numpy.allclose(bounding_polygon[0][0], bounding_polygon[-1][0])) and\
        (numpy.allclose(bounding_polygon[0][1], bounding_polygon[-1][1]))):

        lbp = len(bounding_polygon) - 1
        bounding_polygon = [bounding_polygon[k] for k in range(lbp)]

    for i in range(boundary_segnum - 1):
        # Loop over all the 'other' boundary segments, and find the line which
        # starts at the end point of bounding_polygon

        found_match = False  # Flag for error if no matching segments found

        for j in range(boundary_segnum):
            blj = boundary_line_and_tags[j][0]

            if ((numpy.allclose(blj[0][0], bounding_polygon[-1][0])) and \
                (numpy.allclose(blj[0][1], bounding_polygon[-1][1]))):

                found_match = True

                # Append all but the first point of blj to bounding_polygon

                old_lbp = len(bounding_polygon)
                extra_polygon = [blj[k] for k in range(1, len(blj))]
                bounding_polygon = bounding_polygon + extra_polygon

                # If we are on the last segment, then drop the final point
                # (for consistency with ANUGA format)

                if i == boundary_segnum - 2:
                    # The first and last points should match
                    x_agree = numpy.allclose(bounding_polygon[0][0], 
                                             bounding_polygon[-1][0])
                    y_agree = numpy.allclose(bounding_polygon[0][1], 
                                             bounding_polygon[-1][1])
                    if not (x_agree and y_agree):
                        msg = 'The first and last points of the bounding ' + \
                              'polygon are not identical'
                        raise ValueError(msg)

                    # Remove the last point
                    bounding_polygon = [bounding_polygon[k] for k in
                                        range(len(bounding_polygon) - 1)]
                    new_tags = range(old_lbp - 1, len(bounding_polygon))
                else:
                    new_tags = range(old_lbp - 1, len(bounding_polygon)
                                     - 1)

                # Either the key already exists and we append edge indices, OR
                # we add the key and the edge indices

                key = boundary_line_and_tags[j][1]
                if (key in boundary_tags.keys()):
                    boundary_tags[key] = boundary_tags[key] + new_tags
                else:
                    boundary_tags[key] = new_tags
                # Finish the j loop
                break
        # Check that a 'match' was found for every 'i' iteration
        if not found_match:
            raise Exception('Did not find a match')

    return (bounding_polygon, boundary_tags)
