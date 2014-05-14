'''
Required files are downloaded from the ANUGA servers if they are
out of date or missing.

Modified from the automated_validation_tests/patong/test_patong_scenarios.py code
'''

import sys
import os
import glob
import unittest
import time
import shutil

from anuga.utilities.system_tools \
     import get_web_file, untar_file, file_length, get_host_name
import anuga.utilities.log as log

log.log_filename = './validation.log'

# sourceforge download mirror hosts (must end with '/')
# try these in turn for each file
## NOTE: It would be more reliable if we could somehow 'poll' Sourceforge
##       for a list of mirrors instead of hard-coding a list here.  The only
##       way to do this at the moment is to 'screen-scrape' the data at
##       http://apps.sourceforge.net/trac/sourceforge/wiki/Mirrors
##       but that only gets the main web page for the entity, not the
##       Sourceforge download mirror server.
MIRRORS = ['http://jaist.dl.sourceforge.net/sourceforge/anuga/',          # jp
           'http://transact.dl.sourceforge.net/sourceforge/anuga/',       # au
           'http://voxel.dl.sourceforge.net/sourceforge/anuga/',          # us
           'http://superb-west.dl.sourceforge.net/sourceforge/anuga/',    # u
           'http://dfn.dl.sourceforge.net/sourceforge/anuga/'             # de
          ]

### for testing
##MIRRORS = ['http://10.7.64.243/patong_validation_data/']       # local linux box

# URL to hand-get data files, if required
DATA_FILES_URL = ('http://sourceforge.net/project/showfiles.php?'
                  'group_id=172848&package_id=319323&release_id=677531')

# sequence of mandatory local data objects
Mandatory_Data_Objects = ('data.tgz',)

# sequence of optional local data objects.
# these names must be of the form <scene>.sww.<type>.tgz
# as code below depends upon it.
Optional_Data_Objects = (
                         'patong.sww.TRIAL.tgz',
                         'patong.sww.BASIC.tgz',
                         'patong.sww.FINAL.tgz'
                        )
#Optional_Data_Objects = ('patong.sww.TRIAL.tgz',)

# Associated tolerances to be used in comparisons (would depend on discretisation errors)
epsilon = {'patong.sww.TRIAL.tgz': 1.0e-3,
           'patong.sww.BASIC.tgz': 2.0e-3,
           'patong.sww.FINAL.tgz': 1.0e-3}


# path to the local data directory
Local_Data_Directory = 'local_data'

# path to the remote data directory
Remote_Data_Directory = 'remote_data'

# name of stdout catch file for runmodel.py
RUNMODEL_STDOUT = 'runmodel.stdout'

# text at start of 'output dir' line in RUNMODEL_STDOUT file
OUTDIR_PREFIX = 'Output directory: '

# Name of SWW file produced by run_model.py
OUTPUT_SWW = 'patong.sww'


def setup():
    '''Prepare for the validation run.'''
   
    pass


def refresh_local_data(data_objects, target_dir, mirrors):
    '''Update local data objects from the server.

    data_objects:   list of files to refresh
    target_dir:     directory in which to put files
    mirrors:        list of mirror sites to use
    
    Each file has an associated *.digest file used to decide
    if the local file needs refreshing.
    
    Return True if all went well, else False.
    '''

    # decision function to decide if a file contains HTML
    def is_html(filename):
        '''Decide if given file contains HTML.'''
        
        fd = open(filename)
        data = fd.read(1024)
        fd.close()

        if 'DOCTYPE' in data:
            return True
        
        return False

    
    # local function to get remote file from one of mirrors
    def get_remote_from_mirrors(remote, local, auth, mirrors):
        '''Get 'remote' from one of 'mirrors', put in 'local'.'''

        # Get a unique date+time string to defeat caching.  The idea is to add
        # this to the end of any URL so proxy sees a different request.
        cache_defeat = '?' + time.strftime('%Y%m%d%H%M%S')

        # try each mirror when getting file
        for mirror in mirrors:
            log.debug('Fetching remote file %s from mirror %s'
                      % (remote, mirror))

            remote_url = mirror + remote + cache_defeat
            (result, auth) = get_web_file(remote_url, local, auth=auth)
            if result and is_html(local)==False:
                log.debug('Success fetching file %s' % remote)
                return (True, auth)
            log.debug('Failure fetching from %s' % mirror)
            auth = None

        log.debug('Failure fetching file %s' % remote)
        return (False, auth)            
                

    # local function to compare contents of two files
    def files_same(file_a, file_b):
        '''Compare two files to see if contents are the same.'''
        
        fd = open(file_a, 'r')
        data_a = fd.read()
        fd.close()

        fd = open(file_b, 'r')
        data_b = fd.read()
        fd.close()

        return data_a == data_b

        
    # local function to update one data object
    def refresh_object(obj, auth, mirrors):
        '''Update object 'obj' using authentication tuple 'auth'.
        
        Return (True, <updated_auth>) if all went well,
        else (False, <updated_auth>).
        '''

        # create local and remote file paths.
        obj_digest = obj + '.digest'
        
        remote_file = os.path.join(Remote_Data_Directory, obj)
        remote_digest = remote_file + '.digest'
        
        local_file = os.path.join(Local_Data_Directory, obj)
        local_digest = local_file + '.digest'
        
        # see if missing either digest or object .tgz
        if not os.path.exists(local_digest) or not os.path.exists(local_file):
            # no digest or no object, download both digest and object
            (res, auth) = get_remote_from_mirrors(obj_digest, local_digest, auth, mirrors)
            if res:
                (res, auth) = get_remote_from_mirrors(obj, local_file, auth, mirrors)
        else:
            # download object digest to remote data directory
            (res, auth) = get_remote_from_mirrors(obj_digest, remote_digest, auth, mirrors)
            if res:
                if not files_same(local_digest, remote_digest):
                    # digests differ, refresh object
                    shutil.move(remote_digest, local_digest)
                    (res, auth) = get_remote_from_mirrors(obj, local_file, auth, mirrors)

        return (res, auth)

    # create local data directory if required
    log.debug('Creating local directory: %s' % Local_Data_Directory)
    if not os.path.exists(Local_Data_Directory):
        os.mkdir(Local_Data_Directory)

    # clean out remote data copy directory
    log.debug('Cleaning remote directory: %s' % Remote_Data_Directory)
    shutil.rmtree(Remote_Data_Directory, ignore_errors=True)
    os.mkdir(Remote_Data_Directory)

    # success, refresh local files
    auth = None
    result = True
    for data_object in data_objects:
        log.info("Refreshing file '%s'" % data_object)
        log.debug('refresh_local_data: getting %s from mirrors, auth=%s'
                  % (data_object, str(auth)))
        (res, auth) = refresh_object(data_object, auth, mirrors)
        log.debug('refresh_local_data: returned (res,auth)=%s,%s'
                  % (str(res), str(auth)))
        if res == False:
            log.info('Refresh of file %s failed.' % data_object)
            result = False
            # don't use possibly bad 'auth' again,
            # some proxies lock out on repeated failures.
            auth = None

    if result:
        log.critical('Local data has been refreshed.')
    else:
        log.critical('Local data has been refreshed, with one or more errors.')
    log.critical()
    return result


def can_we_run():
    '''Decide if we can run with the files we have.
    
    Return True if we *can* run, else False.

    Tell user what is happening first, then untar files.
    '''

    log.critical('Checking if you have the required files to run:')

    # get max width of object name string
    max_width = 0
    for obj in Mandatory_Data_Objects:
        max_width = max(len(obj), max_width)
    for obj in Optional_Data_Objects:
        max_width = max(len(obj), max_width)

    # if we don't have *all* mandatory object, can't run
    have_mandatory_files = True
    for obj in Mandatory_Data_Objects:
        obj_path = os.path.join(Local_Data_Directory, obj)
        if os.path.exists(obj_path):
            log.info('\t%s  found' % obj.ljust(max_width))
        else:
            log.info('\t%s  MISSING AND REQUIRED!' % obj.ljust(max_width))
            have_mandatory_files = False

    # at least *one* of these must exist
    have_optional_files = False
    for obj in Optional_Data_Objects:
        obj_path = os.path.join(Local_Data_Directory, obj)
        if os.path.exists(obj_path):
            have_optional_files = True
            log.info('\t%s  found' % obj.ljust(max_width))
        else:
            log.info('\t%s  MISSING!' % obj.ljust(max_width))

    if not have_mandatory_files or not have_optional_files:
        log.critical('You must obtain the missing files before you can run '
                     'this validation.')
        return False

    log.critical('You have enough required files to run.')
    log.critical()

    return True


def set_environment():
    # modify environment so we use the local data
    new_inundationhome = os.path.join(Local_Data_Directory, '')
    os.environ['INUNDATIONHOME'] = new_inundationhome
    new_muxhome = os.path.join(Local_Data_Directory, 'data')
    os.environ['MUXHOME'] = new_muxhome


def run_simulation(vtype, sim_obj):
    '''Run a simulation.

    Returns True if all went well, else False.
    '''
    
    # untar the object
    tar_path = os.path.join(Local_Data_Directory, sim_obj)
    log.info('Untarring %s in directory %s ...'
             % (tar_path, Local_Data_Directory))
    untar_file(tar_path, target_dir=Local_Data_Directory)

    # modify project.py template
    log.debug("Creating '%s' version of project.py" % vtype)
    fd = open('project_template.py', 'r')
    project = fd.readlines()
    fd.close()

    new_project = []
    for line in project:
        new_project.append(line.replace('#!SETUP!#', vtype.lower()))
            
    fd = open('project.py', 'w')
    fd.write(''.join(new_project))
    fd.close()
    
    # import new project.py
    import project

    # run the simulation, produce SWW file
    log.info('Running the simulation ...')
    cmd = 'python run_model.py > %s' % RUNMODEL_STDOUT
    log.debug("run_simulation: doing '%s'" % cmd)
    res = os.system(cmd)
    log.debug("run_simulation: res=%d" % res)

    # 'unimport' project.py
    del project

    # check result
    if res != 0:
        log.critical('Simulation failed, check log')

    return res == 0

def check_that_output_is_as_expected(expected_sww, valid_sww, epsilon):
    '''Check that validation output is as required.'''

    # get path to expected SWW file
    log.critical('Checking that simulation results are as expected ...')
    local_sww = os.path.join(Local_Data_Directory, valid_sww)

    # get output directory from stdout capture file
    try:
        fd = open(RUNMODEL_STDOUT, 'r')
    except IOError, e:
        log.critical("Can't open catch file '%s': %s"
                     % (RUNMODEL_STDOUT, str(e)))
        return 1
    lines = fd.readlines()
    fd.close

    output_directory = None
    for line in lines:
        if line.startswith(OUTDIR_PREFIX):
            output_directory = line.replace(OUTDIR_PREFIX, '', 1)
            output_directory = output_directory.strip()
            break
    if output_directory is None:
        log.critical("Couldn't find line starting with '%s' in file '%s'"
                     % (OUTDIR_PREFIX, RUNMODEL_STDOUT))
        return 1

    log.debug('check_that_output_is_as_expected: output_directory=%s'
              % output_directory)
    
    # compare SWW files here and there
    new_output_sww = os.path.join(output_directory, expected_sww)
    #cmd = 'python cmpsww.py %s %s > cmpsww.stdout' % (local_sww, new_output_sww)
    cmd = 'python compare_model_timeseries.py %s %s %e > compare_model_timeseries.stdout' %\
          (local_sww, new_output_sww, epsilon)
    print '-------------------------------------'
    print cmd
    print '-------------------------------------'    
    
    log.debug("check_that_output_is_as_expected: doing '%s'" % cmd)
    res = os.system(cmd)
    log.debug("check_that_output_is_as_expected: res=%d" % res)
    log.critical()
    print 'Result', res
    if res == 0:
        log.info('Simulation results are as expected.')
    else:
        log.critical('Simulation results are NOT as expected.')
        fd = open('compare_model_timeseries.stdout', 'r')
        cmp_error = fd.readlines()
        fd.close()
        log.critical(''.join(cmp_error))


def teardown():
    '''Clean up after validation run.'''

    log.debug('teardown: called')
    
    # remove remote directory and stdout capture file
    #shutil.rmtree(Remote_Data_Directory, ignore_errors=True)
    #try:
    #    os.remove(RUNMODEL_STDOUT)
    #except OSError:
    #    pass
            

################################################################################
# Mainline - run the simulation, check output.
################################################################################

# set logging levels
log.console_logging_level = log.INFO
log.log_logging_level = log.DEBUG

log.debug('Machine we are running on is "%s"' % get_host_name())
setup()

# prepare user for what is about to happen
log.critical('''
This validation requires a working internet connection to refresh its files.
You may still run this validation without an internet connection if you have the
required files.

If you are behind a proxy server you will need to supply your proxy details
such as the proxy server address and your proxy username and password.  These
can be defined in one or more of the environment variables:
    HTTP_PROXY
    PROXY_USERNAME
    PROXY_PASSWORD
if you wish.  If not supplied in environment variables you will be prompted for
the information.
''')


# make sure local data is up to date
all_objects = Mandatory_Data_Objects + Optional_Data_Objects
if not refresh_local_data(all_objects, Local_Data_Directory, MIRRORS):
    if not can_we_run():
        log.critical("Can't refresh via the internet and you don't have the "
                     "required files.")
        log.critical('Terminating the validation.')
        log.critical('')
        log.critical('If you get the missing files from %s' % DATA_FILES_URL)
        log.critical('then you can try to run the validation again.  Put the '
                     'files into the directory')
        log.critical("%s." % Local_Data_Directory)
        sys.exit(10)

# now untar mandatory objects
for obj in Mandatory_Data_Objects:
    tar_path = os.path.join(Local_Data_Directory, obj)
    log.info('Untarring %s in directory %s ...'
             % (tar_path, Local_Data_Directory))
    untar_file(tar_path, target_dir=Local_Data_Directory)

# set required environment variables
set_environment()

## 6/6/2013 -- Removed this code, as purpose now is just to get data
## now run what simulations we can and check output is as expected

#for odo in Optional_Data_Objects:
#    start_time = time.time()
#
#    _, vtype, _ = odo.rsplit('.', 2)
#    vtype = vtype.lower()
#    log.critical('#' * 72)
#    log.critical('Running Patong "%s" validation ...' % vtype)
#    if run_simulation(vtype, odo):
#        # get SWW names expected and valid, check 'equal'
#        valid_sww, _ = odo.rsplit('.', 1)
#        expected_sww, _ = valid_sww.rsplit('.', 1)
#        check_that_output_is_as_expected(expected_sww, valid_sww, epsilon[odo])
#
#    stop_time = time.time()
#    log.critical('"%s" validation took %.1fs\n\n\n' % (vtype, stop_time - start_time))
#
## clean up
#log.critical('Tearing down ...')
#teardown()
