#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose, Float
from tempfile import NamedTemporaryFile
import os

from data_audit import *

class Test_data_audit(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def NOtest_license_file_is_not_valid(self):
	"""Basic test using an invalid XML file
	"""

        # FIXME(OLE): Needs work to ensure that the order of
        # problems is deterministic. Currently we check for checksum
        # but on some systems file or publishable may come first
        

	# Generate invalid example
	
        fid = NamedTemporaryFile(mode='w',
                                 suffix='.asc',
                                 dir='.')
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.flush()
	
	# Create associated license file
	basename, ext = os.path.splitext(fid.name)
	license_filename = basename + '.lic'
    
        #print fid.name, license_filename
	licfid = open(license_filename, 'w')
	xml_string = """<?xml version="1.0" encoding="iso-8859-1"?>

  <ga_license_file>
    <metadata>
      <author>Ole Nielsen</author>
      <svn_keywords>
        <author>$Author: ole $</author>  
        <date>$Date: 2008-01-21 18:58:15 +1100 (Mon, 21 Jan 2008) $</date>
        <revision>$Revision$</revision>
        <url>$URL: https://datamining.anu.edu.au/svn/ga/anuga_core/source/anuga/utilities/mainland_only.lic $</url>
        <id>$Id: mainland_only.lic 4963 2008-01-21 07:58:15Z ole $</id>
      </svn_keywords>
    </metadata>
    <datafile>
      <filename>mainland_only.csv</filename>
      <checksum>-1661725548</checksum>
      <publishable>No</publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a polygon comprising easting and northing locations 
      tracing parts of the coastline at Dampier WA as well as a rectangular area inland.
      This is used to specifically set the onshore initial condition in a tsunami scenario
      and here, it is used with a unit test in test_polygon.py.
      
      The coastline was derived from Maritime Boundaries which is a public dataset. However,
      rumour has it that some of it was digitised from a Landgate supplied image.
      
      The origin and license issues are still undecided</IP_info>
    </datafile>

  </ga_license_file>
"""
	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
	#print licfid.read()
	
        try:
            license_file_is_valid(licfid)
        except CRCMismatch:
	    pass
	else:
	    msg = 'Should have raised bad CRC exception' 
	    raise Exception, msg    	
	    	
        # Clean up
	licfid.close()
	fid.close()
        os.remove(license_filename)
	

    def NOtest_license_file_is_valid(self):
	"""Basic test using an valid XML file
	"""

        # FIXME(Ole): NOT FINISHED
	
	# Generate valid example
	
        fid = NamedTemporaryFile(mode='w',
                                 suffix='.asc',
                                 dir='.')
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.flush()
	
        # Strip leading dir (./)
	data_filename = os.path.split(fid.name)[1]
	
	print 'Name', data_filename
	
	# Create associated license file
	basename, ext = os.path.splitext(fid.name)
	license_filename = basename + '.lic'
    
	licfid = open(license_filename, 'w')
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
      <filename>%s</filename>
      <checksum>%s</checksum>
      <publishable>Yes</publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a test</IP_info>
    </datafile>

  </ga_license_file>
""" %(data_filename, '000')

	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
	#print licfid.read()
	
        #print fid.name, license_filename
        
	print os.listdir('.')
        license_file_is_valid(licfid, verbose=True)
	    	
        # Clean up
	licfid.close()
	fid.close()
        os.remove(license_filename)
	
		
	
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_data_audit, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)




