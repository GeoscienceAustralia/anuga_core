"""Track IP of data files included in this distribution. 
"""

from os import remove, walk, sep
from os.path import join, splitext

from anuga.utilities.xml_tools import xml2object, XML_element
from anuga.utilities.system_tools import compute_checksum

from data_audit_config import extensions_to_ignore, directories_to_ignore, files_to_ignore



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

def IP_verified(directory, verbose=False):
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
    

    """

    # Print header
    dirwidth = 72

    # Identify data files
    first_time = True
    all_files_accounted_for = True
    for dirpath, datafile in identify_datafiles(directory):
        
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
                license_file_is_valid(fid, dirpath, verbose=verbose)
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


def identify_datafiles(root):
    """ Identify files that might contain data
    """

    for dirpath, dirnames, filenames in walk(root):

        for ignore in directories_to_ignore:
            if ignore in dirnames:
                dirnames.remove(ignore)  # don't visit ignored directories

        #print 'Searching dir', dirpath
        

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


def license_file_is_valid(fid, dirpath='.', verbose=False):
    """Check that XML license file is valid
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
    
    #svn_keywords = metadata['svn_keywords']
    #if verbose: print 'SVN keywords:   ', svn_keywords
    
        
    # Extract information for datafile sections
    datafile = elements['datafile']
    if isinstance(datafile, XML_element):
        datafile = [datafile]

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
        publishable = data['publishable'].upper()
        if verbose: print 'Publishable: "%s"' %publishable        
        if publishable != 'YES':
            msg = 'Data file %s is not flagged as publishable'\
                  %fid.name
            raise NotPublishable, msg



    # If we get this far, the license file is OK
    return True
