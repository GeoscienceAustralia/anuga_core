"""Basic XML utilities based on minidom - the built in Document Object Model
"""

import sys
from xml.dom import minidom, Node
#from xml.sax import make_parser, parse as validate, handler

def print_tree(n, indent=0):
    while n:
        #print 'nodeType', n.nodeType, Node.ELEMENT_NODE
        #if n.nodeType != Node.ELEMENT_NODE:
        #    break

        print ' '*indent,\
              'Node name: "%s",' %n.nodeName,\
              'Node type: "%s",' %n.nodeType,\
              'Node value: "%s"' %str(n.nodeValue).strip()
              
        
        print_tree(n.firstChild, indent+4)
        n = n.nextSibling


def parse(fid):
    """Parse XML file descriptor and return DOM object.
    """

    # FIXME (OLE): XML code should be validated against the DTD
    #validate(fid, handler)
    #doc = minidom.parse(fid, make_parser())

    doc = minidom.parse(fid)    
    return doc


def get_elements(nodelist):
    """Return list of nodes that are ELEMENT_NODE
    """

    element_list = []
    for node in nodelist:
        if node.nodeType == Node.ELEMENT_NODE:
            element_list.append(node)

    return element_list


def get_text(nodelist):
    """Return a concatenation of text fields from list of nodes
    """

    s = ''
    for node in nodelist:
        if node.nodeType == Node.TEXT_NODE:
            s += node.nodeValue + ', '

    if len(s)>0: s = s[:-2]
    return s
