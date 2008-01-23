"""Track IP of data files included in this distribution. 
"""

from os import remove, walk, sep
from os.path import join, splitext

from anuga.utilities.xml_tools import parse, pretty_print_tree, get_elements, get_text
from anuga.utilities.system_tools import compute_checksum

# Audit exceptions
class NotPublishable(Exception): pass
class FilenameMismatch(Exception): pass
class CRCMismatch(Exception): pass
class Invalid(Exception): pass
class WrongTags(Exception): pass

audit_exceptions = (NotPublishable, FilenameMismatch, CRCMismatch, Invalid, WrongTags)

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

    print '---------------------------------------------'
    print 'Files that need to be assessed for IP issues:'
    print '---------------------------------------------'

    # Print header
    dirwidth = 72
    print '---------------------------------------------'
    print 'File'.ljust(dirwidth), 'Status'
    print '---------------------------------------------'

    # Identify data files
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
                status = 'LICENSE FILE NOT VALID'
                status += 'REASON: %s' %e

                #doc = parse(fid)
                #pretty_print_tree(doc)
                fid.seek(0)
                status += fid.read()

            #else:        
            #    if verbose: print 'OK'

            fid.close()
            
        if status != 'OK' or verbose is True:
            print filename + ' (Checksum=%s): '\
                  %str(compute_checksum(filename)), status



    # Return result        
    return all_files_accounted_for


def identify_datafiles(root):
    """ Identify files that might contain data
    """

    # Ignore source code files
    extensions_to_ignore = ['.py','.c','.h', '.f'] #, '.gif', '.jpg', '.png']

    # Ignore generated stuff 
    extensions_to_ignore += ['.pyc', '.o', '.so', '~']
    extensions_to_ignore += ['.aux', '.log', '.idx', 'ilg', '.ind',
                             '.bbl', '.blg']

    # Ignore license files themselves
    extensions_to_ignore += ['.lic']    
    

    # Ignore certain other files
    files_to_ignore = ['README.txt']

    # Ignore directories
    directories_to_ignore = ['anuga_work', 'pymetis', 'obsolete_code',
                             'anuga_parallel', 'anuga_viewer',
                             'planning', 'coding_standards',
                             'experimentation',
                             '.svn', 'misc', '.metadata']

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
    doc = parse(fid)
    #print_tree(doc)

    # Check that file is valid (e.g. all elements there)
    # FIXME (Ole): Todo
    

    if doc.nodeName != '#document':
        msg = 'License file %s does not appear' %license_filename
        msg += 'to be a valid XML document'
        msg += 'The root node has name %s' %doc.nodeName
        msg += 'but it should be %s' %'#document'
        raise Invalid, msg        

    if len(doc.childNodes) != 1:
        msg = 'License file %s must have only one element' %license_filename
        msg += ' at the root level. It is\n '
        msg += '<ga_license_file>'
        raise Invalid, msg
    

    # Start looking at document in earnest
    root_node = doc.childNodes[0]
    if root_node.nodeName != 'ga_license_file':
        msg = 'License file %s must have two elements' %license_filename
        msg += ' at the root level. They are\n '
        msg += '<?xml version="1.0" encoding="iso-8859-1"?>\n'
        msg += '<ga_license_file>\n'
        msg += 'The second element was found to be %s' %root_node.nodeName
        raise WrongTags, msg
    

    # Validate elements: metadata, datafile, datafile, ...
    elements = get_elements(root_node.childNodes)
    if elements[0].nodeName != 'metadata':
        msg = 'The first element under %s must be "metadata"'\
              %root_node.nodeName
        msg += 'The element found was %s' %elements[0].nodeName
        raise WrongTags, msg

    for node in elements[1:]:
        if node.nodeName != 'datafile':
            msg = 'All elements, except the first, under %s must '\
                  %root_node.nodeName            
            msg += 'be "datafile"'
            msg += 'The element found was %s' %node.nodeName
            raise WrongTags, msg        

    if verbose: print    
    # Extract information for source section
    for node in get_elements(elements[0].childNodes):
        if node.nodeName == 'author':
            # Do something
            if verbose: print 'Author:   ', get_text(node.childNodes)

        if node.nodeName == 'svn_keywords':
            # Do nothing
            pass
        
    # Extract information for datafile sections
    for datanode in elements[1:]:
        if verbose: print
    
        for node in get_elements(datanode.childNodes):
            #print 'Node', node.nodeName, node.childNodes
            #continue
            
            if node.nodeName == 'filename':
                # FIXME Check correctness
                filename = join(dirpath, get_text(node.childNodes))
                if verbose: print 'Filename: "%s"' %filename
                try:
                    fid = open(filename, 'r')
                except:
                    msg = 'Specified filename %s could not be opened'\
                          %filename
                    raise FilenameMismatch, msg

            if node.nodeName == 'checksum':
                # FIXME (Ole): This relies on crc being preceded by filename
                reported_crc = get_text(node.childNodes)
                if verbose: print 'Checksum: "%s"' %reported_crc

                file_crc = str(compute_checksum(filename))

                if reported_crc != file_crc:
                    msg = 'Bad checksum (CRC).\n'
                    msg += '  The CRC reported in license file "%s" is "%s"\n'\
                          %(license_filename, reported_crc)
                    msg += '  The CRC computed from file "%s" is "%s"'\
                           %(filename, file_crc)
                    raise CRCMismatch, msg
                

            if node.nodeName == 'accountable':
                accountable = get_text(node.childNodes)
                if verbose: print 'Accountable: "%s"' %accountable
                if accountable == "":
                    msg = 'No accountable person specified'
                    raise Exception, msg

            if node.nodeName == 'source':
                source = get_text(node.childNodes)
                if verbose: print 'Source: "%s"' %source
                if source == "":
                    msg = 'No source specified'
                    raise Exception, msg                

            if node.nodeName == 'IP_owner':
                ip_owner = get_text(node.childNodes)
                if verbose: print 'IP owner: "%s"' %ip_owner
                if ip_owner == "":
                    msg = 'No IP owner specified'
                    raise Exception, msg                                
                

            if node.nodeName == 'IP_info':
                if verbose: print 'IP info: "%s"' %get_text(node.childNodes)  
                

            if node.nodeName == 'publishable':
                
                if verbose: print 'Publishable: %s' %fid.name                
                value = get_text(node.childNodes)
                if value.upper() != 'YES':
                    msg = 'Data file %s is not flagged as publishable'\
                          %fid.name
                    raise NotPublishable, msg



    # If we get this far, the license file is OK
    return True
