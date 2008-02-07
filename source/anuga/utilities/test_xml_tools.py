#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose, Float
from tempfile import mkstemp, mktemp

import os

from xml_tools import *

class Test_xml_tools(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate_xml(self):
	"""Test that xml code is generated from XMLobject model
	"""

        elements = []
        for i in range(10):
            X = XML_element(tag='title %d' %i,
                            contents='example text %d' %i)
            elements.append(X)

            
        X1 = XML_element(tag='first element',
                         contents=XML_element(tag='some text',
                                              contents='hello world'))
        
        X2 = XML_element(tag='second element',
                         contents=XML_element(tag='texts',
                                              contents=elements))
        X3 = XML_element(tag='third element',
                         contents='42')        
        
            

        X = XML_element(tag='all',
                        contents=[X1, X2, X3])
        doc = XML_document(contents=X)

        #print doc.pretty_print()        
        #print doc
        assert doc['all']['second element']['texts']['title 4'].elements == 'example text 4'

        

    def test_xml2object(self):
        """Test that XML_document can be generated from file
        """

        tmp_fd , tmp_name = mkstemp(suffix='.xml', dir='.')
        fid = os.fdopen(tmp_fd, 'w')
        
	xml_string = """<?xml version="1.0" encoding="iso-8859-1"?>
  <ga_license_file>
    <metadata>
      <author>Ole Nielsen</author>
      <svn_keywords>
        <author>$Author$</author>  
        <date>$Date$</date>
        <revision>$Revision$</revision>
        <url>$URL:$</url>
        <id>$Id$</id>
      </svn_keywords>
    </metadata>
    <datafile>
      <filename>bathymetry.asc</filename>
      <checksum>%s</checksum>
      <publishable>Yes</publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a test</IP_info>
    </datafile>
  </ga_license_file>
""" %('1234')

	fid.write(xml_string)
	fid.close()

        fid = open(tmp_name)
	reference = fid.read()
        reflines = reference.split('\n')
	
        xmlobject = xml2object(fid, verbose=True)


        #print xmlobject.pretty_print()
        
        xmllines = str(xmlobject).split('\n')

        #for line in reflines:
        #    print line
        #print    
        #for line in xmllines:
        #    print line            

            
        assert len(reflines) == len(xmllines)

        for i, refline in enumerate(reflines):
            msg = '%s != %s' %(refline.strip(), xmllines[i].strip())
            assert refline.strip() == xmllines[i].strip(), msg

        # Check dictionary behaviour    
        for tag in xmlobject['ga_license_file'].keys():
            xmlobject['ga_license_file'][tag]

        assert xmlobject['ga_license_file']['datafile']['accountable'].elements == 'Jane Sexton'
        

        #print
        #print
        #print xmlobject['ga_license_file']['metadata']
        #print xmlobject['ga_license_file']['datafile']
        #print xmlobject['ga_license_file']['datafile']['accountable']        

        #print xmlobject['ga_license_file'].keys()

        #for tag in xmlobject['ga_license_file'].keys():
        #    print xmlobject['ga_license_file'][tag]
	    	
        # Clean up
	fid.close()

	

	
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_xml_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
        

