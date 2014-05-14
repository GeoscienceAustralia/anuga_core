#!/usr/bin/env python

'''
A program to compare two SWW  files for "equality".

This program makes lots of assumptions about the structure of the SWW files, 
so if that structure changes, this program must change.
'''

import sys
import os
import os.path
import getopt
from anuga.file.netcdf import NetCDFFile
import numpy as num
from anuga.config import netcdf_mode_r


#####
# Various constants.
#####

# Global for the '-q' quiet flag
quiet = None

# default tolerances - same as numpy defaults
default_abs_tolerance = 1.0e-08
default_rel_tolerance = 1.0000000000000001e-05

# Global attributes that should exist and be same in both files
# Don't have to have all of these, and we don't care about others.
expect_global_attributes = ['smoothing', 'vertices_are_stored_uniquely',
                            'order', 'starttime',
                            'xllcorner', 'yllcorner',
                            'zone', 'false_easting', 'false_northing',
                            'datum', 'projection', 'units']

# dimensions expected, with expected values (None means unknown)
expected_dimensions = {'number_of_volumes': None,
                       'number_of_vertices': 3,
                       'numbers_in_range': 2,
                       'number_of_points': None,
                       'number_of_timesteps': None}

# Variables expected, with expected dimensions.
# Don't have to have all of these, and we don't care about others.
expected_variables = {'x': ('number_of_points',),
                      'y': ('number_of_points',),
                      'elevation': ('number_of_points',),
                      'elevation_range': ('numbers_in_range',),
                      'z': ('number_of_points',),
                      'volumes': ('number_of_volumes', 'number of vertices'),
                      'time': ('number_of_timesteps',),
                      'stage': ('number_of_timesteps', 'numbers_of_points',),
                      'stage_range': ('numbers_in_range',),
                      'xmomentum': ('number_of_timesteps', 'number_of_points'),
                      'xmomentum_range': ('numbers_in_range'),
                      'ymomentum': ('number_of_timesteps', 'number_of_points'),
                      'ymomentum_range': ('numbers_in_range')}

##
# @brief An exception to inform user of usage problems.
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


##
# @brief Compare two SWW files.
# @param files A tuple of two filenames.
# @param globals A list of global attribute names to compare.
# @param timesteps A list of timesteps to compare at.
# @param variables A list of variable names to compare.
# @return Returns if files 'equal', else raises RuntimeError.
def files_are_the_same(files, globals=None, timesteps=None, variables=None,
                       rel_tolerance=default_rel_tolerance,
                       abs_tolerance=default_abs_tolerance):
    # split out the filenames and check they exist
    (file1, file2) = files
    filename1 = os.path.basename(file1)
    filename2 = os.path.basename(file2)
#    width = max(len(filename1), len(filename2))
#    filename1 = filename1.rjust(width)
#    filename2 = filename2.rjust(width)

    error = False
    error_msg = ''

    try:
        fid1 = NetCDFFile(file1, netcdf_mode_r)
    except:
        error_msg += "\nFile '%s' can't be opened?\n" % file1
        error = True

    try:
        fid2 = NetCDFFile(file2, netcdf_mode_r)
    except:
        error_msg += "\nFile '%s' can't be opened?\n" % file2
        error = True
        fid1.close()

    if error:
        raise RuntimeError, error_msg

    if globals is None:
        globals = expect_global_attributes

    #####
    # First, check that files have the required structure
    #####

    # dimensions - only check expected dimensions
    for key in expected_dimensions:
        if key not in fid1.dimensions.keys():
            error_msg += ("\nFile %s doesn't contain dimension '%s'\n"
                          % (filename1, key))
            error = True
        if key not in fid2.dimensions.keys():
            error_msg += ("\nFile %s doesn't contain dimension '%s'\n"
                          % (filename2, key))
            error = True

    # now check that dimensions are the same length
    # NOTE: DOESN'T CHECK 'UNLIMITED' DIMENSIONS YET! (get None at the moment)
    for dim in expected_dimensions:
        dim1_shape = fid1.dimensions.get(dim, None)
        dim2_shape = fid2.dimensions.get(dim, None)
        if dim1_shape != dim2_shape:
            error_msg += ("\nFile %s has '%s' dimension of size %s,\n"
                          "file %s has that dimension of size %s\n"
                          % (filename1, dim, str(dim1_shape),
                             filename2, str(dim2_shape)))
            error = True

    # check that we have the required globals
    if globals:
        for glob in globals:
            if glob not in dir(fid1):
                error_msg += ("\nGlobal attribute '%s' isn't in file %s\n"
                              % (glob, filename1))
                error = True
            if glob not in dir(fid2):
                error_msg += ("\nGlobal attribute '%s' isn't in file %s\n"
                              % (glob, filename2))
                error = True
    else:
        # get list of global attributes
        glob_vars1 = []
        glob_vars2 = []
        for glob in expect_global_attributes:
            if glob in dir(fid1):
                glob_vars1.append(glob)
            if glob in dir(fid2):
                glob_vars2.append(glob)

        # now check attribute lists are same
        if glob_vars1 != glob_vars2:
            error_msg = ('\nFiles differ in global attributes:\n'
                         '%s: %s,\n'
                         '%s: %s\n' 
                         % (filename1, str(glob_vars1),
                            filename2, str(glob_vars2)))
            error = True

    # get variables to test
    if variables:
        for var in variables:
            if var not in fid1.variables.keys():
                error_msg += ("\nVariable '%s' isn't in file %s\n"
                              % (var, filename1))
                error = True
            if var not in fid2.variables.keys():
                error_msg += ("\nVariable '%s' isn't in file %s\n"
                              % (var, filename2))
                error = True
    else:
        # check that variables are as expected in both files
        var_names1 = []
        var_names2 = []
        for var_name in expected_variables:
            if fid1.variables.has_key(var_name):
                var_names1.append(var_name)
            if fid2.variables.has_key(var_name):
                var_names2.append(var_name)
    
        if var_names1 != var_names2:
            error_msg += ('\nVariables are not the same between files:\n'
                          '%s variables=%s,\n'
                          '%s variables=%s\n'
                          % (filename1, str(var_names1), filename2, str(var_names2)))
            error = True
        variables = var_names1

    # get size of time dimension
    num_timesteps1 = fid1.variables['time'].shape
    num_timesteps2 = fid2.variables['time'].shape
    if num_timesteps1 != num_timesteps2:
        error_msg += ('Files have different number of timesteps:\n'
                      '%s=%s,\n'
                      '%s=%s\n'
                      % (filename1, str(num_timesteps1),
                         filename2, str(num_timesteps2)))
        error = True

    num_timesteps = num_timesteps1[0]

    # variable shapes same?
    for var_name in variables:
        var1 = fid1.variables[var_name]
        var2 = fid2.variables[var_name]
        var1_shape = var1.shape
        var2_shape = var2.shape
        if var1_shape != var2_shape:
            error_msg += ('Files differ in variable %s shape:\n'
                          '%s: %s,\n'
                          '%s: %s\n'
                          % (var_name, filename1, str(var1_shape),
                                       filename2, str(var2_shape)))
            error = True

    if error:
        fid1.close()
        fid2.close()
        raise RuntimeError, error_msg

    #####
    # Now check that actual data values are the same
    #####

    error_msg = ''
    glob_vars_bad = {}
    data_vars_bad = {}

    # check values of global attributes
    for glob_name in globals:
        if getattr(fid1, glob_name) != getattr(fid2, glob_name):
            print("\nFiles differ in global '%s':\n"
                  "%s: '%s',\n"
                  "%s: '%s'"
                  % (glob_name, filename1, str(g1), filename2, str(g2)))
            glob_vars_bad[glob_name] = glob_vars_bad.get(glob_name, 0) + 1
            error = True

    # check data variables, be clever with time series data
    max_rel_difference = -1
    diff_count = 0
    for var_name in variables:
        var_dims = expected_variables[var_name]
        if (len(var_dims) > 1) and (var_dims[0] == 'number_of_timesteps'):
            # time series, check by timestep block
            for t in xrange(num_timesteps):
                var1 = num.array(fid1.variables[var_name][t,:])
                var2 = num.array(fid2.variables[var_name][t,:])
                if not num.allclose(var1, var2,
                                    rtol=rel_tolerance, atol=abs_tolerance):
                    error = True
                    for i in xrange(len(var1)):
                        if not num.allclose(var1[i], var2[i],
                                            rtol=rel_tolerance,
                                            atol=abs_tolerance):
                            abs_difference = num.abs(var1[i]-var2[i])
                            max_a_b = num.max(num.abs(var1[i]),
                                              num.abs(var2[i]))
                            rel_difference = num.abs(abs_difference/max_a_b)

                            if not quiet:
                                print('\nFiles differ in variable '
                                      '%s[%d,%d]:\n'
                                      '%s: %f\n'
                                      '%s: %f\n'
                                      'abs. difference=%e, rel. difference=%e\n'
                                      % (var_name, t, i,
                                         filename1, var1[i],
                                         filename2, var2[i],
                                         abs_difference,
                                         rel_difference))

                            if rel_difference > max_rel_difference:
                                max_rel_difference = rel_difference
                                max_rel_difference_abs = abs_difference
                                max_rel_difference_a = var1[i]
                                max_rel_difference_b = var2[i]

                            data_vars_bad[var_name] = data_vars_bad.get(var_name, 0) + 1
                            diff_count += 1
        else:
            # simple data, check whole thing at once
            var1 = num.array(fid1.variables[var_name][:])
            var2 = num.array(fid2.variables[var_name][:])
            if not num.allclose(var1, var2,
                                rtol=rel_tolerance, atol=abs_tolerance):
                for j in xrange(len(var1)):
                    if not num.allclose(var1[j], var2[j],
                                          rtol=rel_tolerance, atol=abs_tolerance):
                        abs_difference = num.abs(var1[j]-var2[j])
                        max_a_b = num.max(num.abs(var1[j]),
                                          num.abs(var2[j]))
                        rel_difference = num.abs(abs_difference/max_a_b)

                        if not quiet:
                            print('\nFiles differ in variable '
                                  '%s[%d]:\n'
                                  '%s: %f\n'
                                  '%s: %f\n'
                                  'abs. difference=%e, rel. difference=%e\n'
                                   % (var_name, j,  
                                      filename1, var1[j],
                                      filename2, var2[j],
                                      abs_difference,
                                      rel_difference))

                        if rel_difference > max_rel_difference:
                            max_rel_difference = rel_difference
                            max_rel_difference_abs = abs_difference
                            max_rel_difference_a = var1[j]
                            max_rel_difference_b = var2[j]

                        data_vars_bad[var_name] = data_vars_bad.get(var_name, 0) + 1
                        diff_count += 1
                error = True

    #####
    # close files and signal OK or ERROR
    #####

    fid1.close()
    fid2.close()

    if error:
        error_msg += ('\nNumber of data differences=%d\n'
                      'Maximum relative data difference=%e\n'
                      'associated absolute difference=%e\n'
                      "associated 'a' value=%e\n"
                      "associated 'b' value=%e\n"
                      % (diff_count, max_rel_difference, max_rel_difference_abs,
                         max_rel_difference_a, max_rel_difference_b))
        error_msg += ('\nglob_vars bad=%s\n' % str(glob_vars_bad))
        error_msg += ('\ndata_vars bad=%s\n' % str(data_vars_bad))
        raise RuntimeError, error_msg

    return


##
# @brief Return a usage string.
def usage():
    result = []
    a = result.append
    a('Usage: %s <options> <file1> <file2>\n' % ProgName)
    a('where <options> is zero or more of:\n')
    a('                   -h        print this help\n')
    a('                   -q        be quiet, print only summary of differences\n')
    a("                   -a <val>  set absolute threshold of 'equivalent'\n")
    a("                   -r <val>  set relative threshold of 'equivalent'\n")
    a('                   -g <arg>  check only global attributes specified\n')
    a('                             <arg> has the form <globname>[,<globname>[,...]]\n')
    a('                   -t <arg>  check only timesteps specified\n')
    a('                             <arg> has the form <starttime>[,<stoptime>[,<step>]]\n')
    a('                   -v <arg>  check only the named variables\n')
    a('                             <arg> has the form <varname>[,<varname>[,...]]\n')
    a('and <file1> and <file2> are two SWW files to compare.\n')
    a('\n')
    a('The program exit status is one of:\n')
    a('   0    the two files are equivalent\n')
    a('   else the files are not equivalent.')
    return ''.join(result)

##
# @brief Print a message to stderr.
def warn(msg):
    print >>sys.stderr, msg


##
# @brief 
# @param argv 
# @return The status code the program will exit with.
def main(argv=None):
    global quiet

    if argv is None:
        argv = sys.argv

    try:
        try:
            opts, args = getopt.getopt(argv[1:], 'hqa:g:r:t:v:',
                                       ['help', 'globals',
                                        'variables', 'timesteps'])
        except getopt.error, msg:
            raise Usage(msg)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    # process options
    globals = None
    timesteps = None
    variables = None
    quiet = False
    rel_tolerance = default_rel_tolerance
    abs_tolerance = default_abs_tolerance
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print usage()
            sys.exit(0)
        elif opt in ('-q', '--quiet'):
            quiet = True
        elif opt in ('-a', '--absolute'):
            abs_tolerance = float(arg)
        elif opt in ('-r', '--relative'):
            rel_tolerance = float(arg)
        elif opt in ('-g', '--globals'):
            globals = arg.split(',')
        elif opt in ('-t', '--timesteps'):
            timesteps = arg.split(',')
        elif opt in ('-v', '--variables'):
            variables = arg.split(',')

    # process arguments
    if len(args) != 2:
        msg = usage()
        print 'msg=%s' % msg
        raise Usage(msg)

    try:
        files_are_the_same(args, globals=globals,
                           timesteps=timesteps, variables=variables,
                           rel_tolerance=rel_tolerance,
                           abs_tolerance=abs_tolerance)
    except RuntimeError, msg:
         print msg
         return 10


if __name__ == "__main__":
    global ProgName

    ProgName = os.path.basename(sys.argv[0])

    sys.exit(main())

