"""Track IP of data files in an entire directory tree.
See docstring for the public function IP_verified()
for details.
"""

from os import remove, walk, sep
from os.path import join, splitext

from anuga.utilities.xml_tools import xml2object, XML_element
from anuga.utilities.system_tools import compute_checksum


# Audit exceptions
class NotPublishable(Exception): pass
class FilenameMismatch(Exception): pass
class CRCMismatch(Exception): pass
class Invalid(Exception): pass
class WrongTags(Exception): pass

audit_exceptions = (NotPublishable,
                    FilenameMismatch,
                    CRCMismatch,
                    Invalid,
                    WrongTags)


def IP_verified(directory,
                extensions_to_ignore=None,
                directories_to_ignore=None,
                files_to_ignore=None,
                verbose=False):
    """Find and audit potential data files that might violate IP

    This is the public function to be used to ascertain that
    all data in the specified directory tree has been audited according
    to the GA data IP tracking process.

    if IP_verified is False:
        # Stop and take remedial action
        ...
    else:
        # Proceed boldly with confidence
        
    verbose controls standard output.
    If verbose is False, only diagnostics about failed audits will appear.
    All files that check OK will pass silently.

    Optional arguments extensions_to_ignore, directories_to_ignore, and
    files_to_ignore are lists of things to skip.

    Examples are:
    extensions_to_ignore = ['.py','.c','.h', '.f'] # Ignore source code
    files_to_ignore = ['README.txt']
    directories_to_ignore = ['.svn', 'misc']

    None is also OK for these parameters.
    
    """

    # Print header
    dirwidth = 72

    # Identify data files
    first_time = True
    all_files_accounted_for = True
    for dirpath, datafile in identify_datafiles(directory,
                                                extensions_to_ignore,
                                                directories_to_ignore,
                                                files_to_ignore):
        
        filename = join(dirpath, datafile)
        basename, ext = splitext(datafile)

        # Look for a XML license file with the .lic
        status = 'OK'
        try:
            fid = open(join(dirpath, basename + '.lic'))
        except IOError:
            status = 'NO LICENSE FILE'
            all_files_accounted_for = False
        else:
            try:
                license_file_is_valid(fid, datafile, dirpath,
                                      verbose=False)
            except audit_exceptions, e:
                all_files_accounted_for = False                                
                status = 'LICENSE FILE NOT VALID\n'
                status += 'REASON: %s\n' %e

                try:
                    doc = xml2object(fid)
                except:
                    status += 'XML file could not be read:'
                    fid.seek(0)
                    status += fid.read()                    
                else:    
                    status += str(doc)

            fid.close()
            
        if status != 'OK' or verbose is True:
            if first_time is True:
                # Print header
                print '---------------------------------------------'
                print 'Files that need to be assessed for IP issuses'.ljust(dirwidth), 'Status'
                print '---------------------------------------------'
                first_time = False

            print filename + ' (Checksum=%s): '\
                  %str(compute_checksum(filename)), status



    # Return result        
    return all_files_accounted_for



#------------------
# Private functions
#------------------
def identify_datafiles(root,
                       extensions_to_ignore=None,
                       directories_to_ignore=None,
                       files_to_ignore=None):
    """ Identify files that might contain data

    See function IP_verified() for details about optinoal parmeters
    """

    for dirpath, dirnames, filenames in walk(root):

        for ignore in directories_to_ignore:
            if ignore in dirnames:
                dirnames.remove(ignore)  # don't visit ignored directories


        for filename in filenames:


            # Ignore extensions that need no IP check
            ignore = False
            for ext in extensions_to_ignore:
                if filename.endswith(ext):
                    ignore = True

            if filename in files_to_ignore:
                ignore = True

            if ignore is False:
                yield dirpath, filename


def license_file_is_valid(fid, filename_to_verify,
                          dirpath='.', verbose=False):
    """Check that XML license file for given filename_to_verify is valid.

    Input:
        fid: Open file object for XML license file
        file_name_to_verify: The data filename that is being audited
        dir_path: Where the files live
        verbose: Optional verbosity
        

    Check for each datafile listed that

    * Datafile tags are there and match the one specified
    * Fields are non empty
    * Datafile exists
    * Checksum is correct
    * Datafile is flagged as publishable

    If anything is violated an appropriate exception is raised.
    If everything is honky dory the function will return True.
    """

    license_filename = fid.name

    doc = xml2object(fid)
    #print doc

    
    # Check that file is valid (e.g. all elements there)
    if not doc.has_key('ga_license_file'):
        msg = 'License file %s must have two elements' %license_filename
        msg += ' at the root level. They are\n'
        msg += '  <?xml version="1.0" encoding="iso-8859-1"?>\n'
        msg += '  <ga_license_file>\n'
        msg += 'The second element was found to be %s' %doc.keys()
        raise WrongTags, msg
    

    # Validate elements: metadata, datafile, datafile, ...
    # FIXME (Ole): I'd like this to verified by the parser
    # using a proper DTD template one day....
    # For not, let's check the main ones.
    elements = doc['ga_license_file']
    if not elements.has_key('metadata'):
        msg = 'Tag %s must have the element "metadata"'\
              %doc.keys()[0]
        msg += 'The element found was %s' %elements[0].nodeName
        raise WrongTags, msg

    if not elements.has_key('datafile'):
        msg = 'Tag %s must have the element "datafile"'\
              %doc.keys()[0]
        msg += 'The element found was %s' %elements[0].nodeName
        raise WrongTags, msg    

    for key in elements.keys():
        msg = 'Invalid tag: %s' %key
        if not key in ['metadata', 'datafile']:
            raise WrongTags, msg                    

    
    # Extract information for metadata section
    if verbose: print
    metadata = elements['metadata']

    author = metadata['author']
    if verbose: print 'Author:   ', author
    if author == '':
        msg = 'Missing author'
        raise Exception, msg                
    
    #svn_keywords = metadata['svn_keywords']
    #if verbose: print 'SVN keywords:   ', svn_keywords
    
        
    # Extract information for datafile sections
    datafile = elements['datafile']
    if isinstance(datafile, XML_element):
        datafile = [datafile]


    # Check that filename to verify is listed in license file
    found = False
    for data in datafile:    
        if data['filename'] == filename_to_verify:
            found = True
    if not found:
        msg = 'Specified filename to verify %s ' %filename_to_verify
        msg += 'did not appear in license file %s' %license_filename
        raise FilenameMismatch, msg                
            
        
    # Check contents
    for data in datafile:
        if verbose: print

        # Filename
        if data['filename'] == '':
            msg = 'Missing filename'
            raise FilenameMismatch, msg            
        else:
            filename = join(dirpath, data['filename'])
            if verbose: print 'Filename: "%s"' %filename
            try:
                fid = open(filename, 'r')
            except:
                msg = 'Specified filename %s could not be opened'\
                      %filename
                raise FilenameMismatch, msg

        # CRC
        reported_crc = data['checksum']
        if verbose: print 'Checksum: "%s"' %reported_crc
        
        file_crc = str(compute_checksum(filename))
        if reported_crc != file_crc:
            msg = 'Bad checksum (CRC).\n'
            msg += '  The CRC reported in license file "%s" is "%s"\n'\
                   %(license_filename, reported_crc)
            msg += '  The CRC computed from file "%s" is "%s"'\
                   %(filename, file_crc)
            raise CRCMismatch, msg
                
        # Accountable
        accountable = data['accountable']
        if verbose: print 'Accountable: "%s"' %accountable
        if accountable == '':
            msg = 'No accountable person specified'
            raise Exception, msg

        # Source
        source = data['source']
        if verbose: print 'Source: "%s"' %source
        if source == '':
            msg = 'No source specified'
            raise Exception, msg                

        # IP owner
        ip_owner = data['IP_owner']
        if verbose: print 'IP owner: "%s"' %ip_owner
        if ip_owner == '':
            msg = 'No IP owner specified'
            raise Exception, msg                                
                
        # IP info
        ip_info = data['IP_info']
        if verbose: print 'IP info: "%s"' %ip_info
        if ip_info == '':
            msg = 'No IP info specified'
            raise Exception, msg                                               

        # Publishable
        publishable = data['publishable']
        if verbose: print 'Publishable: "%s"' %publishable
        if publishable == '':
            msg = 'No publishable value specified'
            raise NotPublishable, msg
        
        if publishable.upper() != 'YES':
            msg = 'Data file %s is not flagged as publishable'\
                  %fid.name
            raise NotPublishable, msg



    # If we get this far, the license file is OK
    return True
