"""Track IP of data files included in this distribution. 
"""

from os import remove, walk, sep
from os.path import join, splitext

from anuga.utilities.xml_tools import parse, print_tree, get_elements, get_text

# Audit exceptions
class NotPublishable(Exception): pass
class Invalid(Exception): pass
class WrongTags(Exception): pass


def IP_verified(directory):
    """Find and audit potential data files that might violate IP

    This is the public function to be used to ascertain that
    all data in the specified directory tree has been audited according
    to the GA data IP tracking process.

    if IP_verified is False:
        # Stop and take remedial action
        ...
    else:
        # Proceed boldly with confidence
        

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
        print join(dirpath, datafile) + ': ',

        basename, ext = splitext(datafile)

        # Look for a XML license file with the .lic
        try:
            fid = open(join(dirpath, basename + '.lic'))
        except IOError:
            print 'NO LICENSE FILE'
            all_files_accounted_for = False
        else:
            if license_file_is_valid(fid):
                print 'OK'
            else:
                print 'LICENSE FILE NOT VALID'
                all_files_accounted_for = False
            fid.close()

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


def license_file_is_valid(fid):
    """Check that XML license file is valid
    """

    doc = parse(fid)
    #print_tree(doc)

    # Check that file is valid (e.g. all elements there)
    # FIXME (Ole): Todo
    

    if doc.nodeName != '#document':
        msg = 'License file %s does not appear' %fid.name
        msg += 'to be a valid XML document'
        msg += 'The root node has name %s' %doc.nodeName
        msg += 'but it should be %s' %'#document'
        raise Invalid, msg        

    if len(doc.childNodes) != 2:
        msg = 'License file %s must have two elements' %fid.name
        msg += ' at the root level. They are\n '
        msg += '<?xml version="1.0" encoding="iso-8859-1"?>\n'
        msg += '<ga_license_file>'
        raise Invalid, msg
    

    # Start looking at document in earnest
    root_node = doc.childNodes[1]
    if root_node.nodeName != 'ga_license_file':
        msg = 'License file %s must have two elements' %fid.name
        msg += ' at the root level. They are\n '
        msg += '<?xml version="1.0" encoding="iso-8859-1"?>\n'
        msg += '<ga_license_file>\n'
        msg += 'The second element was found to be %s' %root_node.nodeName
        raise WrongTags, msg
    

    # Validate elements: source, datafile, datafile, ...
    elements = get_elements(root_node.childNodes)
    if elements[0].nodeName != 'source':
        msg = 'The first element under %s must be "source"'\
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

    print 
    # Extract information for source section
    for node in get_elements(elements[0].childNodes):
        if node.nodeName == 'author':
            # Do something
            print 'Author is', get_text(node.childNodes)

        if node.nodeName == 'svn_keywords':
            # Do nothing
            pass
        
    # Extract information for datafile sections
    for datanode in elements[1:]:
        print    
    
        for node in get_elements(datanode.childNodes):
            #print 'Node', node.nodeName, node.childNodes
            #continue
            
            if node.nodeName == 'filename':
                # FIXME Check correctness
                print 'Filename is "%s"' %get_text(node.childNodes)

            if node.nodeName == 'accountable':
                print 'Accountable is "%s"' %get_text(node.childNodes)

            if node.nodeName == 'location':
                print 'Location is "%s"' %get_text(node.childNodes)

            if node.nodeName == 'IP_owner':
                print 'IP owner is "%s"' %get_text(node.childNodes)

            if node.nodeName == 'IP_info':
                print 'IP info is "%s"' %get_text(node.childNodes)                                
                

            if node.nodeName == 'publishable':
                value = get_text(node.childNodes)
                if value.upper() != 'YES':
                    msg = 'Data file %s is not flagged as publishable'\
                          %fid.name
                    print msg
                    #raise NotPublishable, msg
                else:
                    print 'Data file %s is flagged publishable' %fid.name                

    #FIXME (Ole): Use hash code for original datafile as an XML element
    # USE CRC32 in zlib or hash
    
    #for node in elements:
    #    print node
    #print 



    # Check that file is deemed publishable
    items = doc.getElementsByTagName('publishable')
    for i in items:
        print i
        #i.getAttribute()
