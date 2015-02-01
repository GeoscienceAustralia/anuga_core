#!/usr/bin/env python

import unittest
from tempfile import mkstemp
import os

from anuga.utilities.data_audit import *


class Test_data_audit(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_license_file_is_not_valid1(self):
	"""Basic test using an invalid XML file. This one
        should fail on bad CRC checksum
	"""

	# Generate invalid checksum example
	
        tmp_fd , tmp_name = mkstemp(suffix='.asc', dir='.')
        fid = os.fdopen(tmp_fd, 'w')
        
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()
	
	# Create associated license file
	basename, ext = os.path.splitext(tmp_name)
	license_filename = basename + '.lic'
    
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
      <filename>%s</filename>
      <checksum>-111111</checksum>
      <publishable>Yes</publishable>
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
""" %tmp_name
        
	licfid.write(xml_string)
	licfid.close()

	#licfid = open(license_filename)
	#print licfid.read()
        #licfid.close()

        try:
            license_file_is_valid(license_filename, tmp_name)
        except CRCMismatch:
	    pass
	else:
	    msg = 'Should have raised bad CRC exception' 
	    raise Exception, msg    	
	    	
        # Clean up
	#licfid.close()
        os.remove(license_filename)
        try:
            os.remove(tmp_name)        
        except:
            # FIXME(DSG) Windows seems to have a problem deleting this file
            # This is a work-a-round. It doesn't fix the root problem
            # It does delete the file though.
            fid = open(tmp_name, 'a')        
            string = 'Example data file'
            fid.write(string)
            fid.close()
            os.remove(tmp_name) 
        



    def test_license_file_is_not_valid2(self):
	"""Basic test using an invalid XML file. This one
        should fail on Not Publishable
	"""

	# Generate invalid checksum example
	
        tmp_fd , tmp_name = mkstemp(suffix='.asc', dir='.')
        fid = os.fdopen(tmp_fd, 'w')
        
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()
	
	# Create associated license file
	basename, ext = os.path.splitext(tmp_name)
	license_filename = basename + '.lic'
    
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
      <filename>%s</filename>
      <checksum>2810517858</checksum>
      <publishable>no</publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a polygon comprising easting and northing locations</IP_info>
    </datafile>

  </ga_license_file>
""" %tmp_name
        
	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
	#print licfid.read()


        try:
            license_file_is_valid(licfid, tmp_name)
        except NotPublishable:
	    pass
	else:
	    msg = 'Should have raised NotPublishable exception' 
	    raise Exception, msg    	
	    	
        # Clean up
	licfid.close()
        os.remove(license_filename)

	fid.close()        
        try:
            os.remove(tmp_name)        
        except:
            # FIXME(DSG) Windows seems to have a problem deleting this file
            # This is a work-a-round. It doesn't fix the root problem
            # It does delete the file though.
            fid = open(tmp_name, 'a')        
            string = 'Example data file'
            fid.write(string)
            fid.close()
            os.remove(tmp_name) 




    def test_license_file_is_not_valid3(self):
	"""Basic test using an invalid XML file. This one
        should fail on Filename Mismatch
	"""

	
        tmp_fd , tmp_name = mkstemp(suffix='.asc', dir='.')
        fid = os.fdopen(tmp_fd, 'w')
        
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()
	
	# Create associated license file
	basename, ext = os.path.splitext(tmp_name)
	license_filename = basename + '.lic'
    
	licfid = open(license_filename, 'w')
	xml_string = """<?xml version="1.0" encoding="iso-8859-1"?>

  <ga_license_file>
    <metadata>
      <author>Ole Nielsen</author>
      <svn_keywords>
        <author>$Author: ole $</author>  
        <date>$Date: 2008-01-21 18:58:15 +1100 (Mon, 21 Jan 2008) $</date>
        <revision>$Revision$</revision>
        <url>$URL:$</url>
        <id>$Id:$</id>
      </svn_keywords>
    </metadata>
    <datafile>
      <filename>%s</filename>
      <checksum>2810517858</checksum>
      <publishable>Yes</publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a polygon comprising easting and northing locations</IP_info>
    </datafile>

  </ga_license_file>
""" %(basename + '.no_exist')

        
	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
	#print licfid.read()


        try:
            license_file_is_valid(licfid, basename + '.no_exist')
        except FilenameMismatch:
	    pass
	else:
	    msg = 'Should have raised FilenameMismatch exception' 
	    raise Exception, msg    	
	    	
        # Clean up
	licfid.close()
	fid.close()
        os.remove(license_filename)
        os.remove(tmp_name)        




    def test_license_file_is_valid(self):
	"""Basic test using an valid XML file
	"""
	
	# Generate valid example
        tmp_fd , tmp_name = mkstemp(suffix='.asc', dir='.')
        fid = os.fdopen(tmp_fd, 'w')        

        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()
	
        # Strip leading dir (./)
	data_filename = os.path.split(tmp_name)[1]
	
	#print 'Name', data_filename
	
	# Create associated license file
	basename, ext = os.path.splitext(tmp_name)
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
""" %(data_filename, '2810517858')

	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
        license_file_is_valid(licfid, data_filename)
        licfid.close()
        
        # Clean up
        os.remove(license_filename)
        os.remove(tmp_name)        
	



    def test_valid_license_file_with_multiple_files(self):
	"""Test of XML file with more than one datafile element.
	"""
	
	# Generate example files
        tmp_fd , tmp_name = mkstemp(suffix='.asc', dir='.')
        fid = os.fdopen(tmp_fd, 'w')        
        string = 'Example data file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()

        # Derive filenames
	basename, ext = os.path.splitext(tmp_name)
	data_filename1 = basename + '.asc' 
        data_filename2 = basename + '.prj'
	license_filename = basename + '.lic'
        #print data_filename1, data_filename2, license_filename         

        # Write data to second data file
        fid = open(data_filename2, 'w')        
        string = 'Another example data file with text in it'
        fid.write(string)
        fid.close()        

        # Create license file
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
      <source>Generated on the fly</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a test</IP_info>
    </datafile>
    <datafile>
      <filename>%s</filename>
      <checksum>%s</checksum>
      <publishable>Yes</publishable>
      <accountable>Ole Nielsen</accountable>
      <source>Generated on the fly</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is another test</IP_info>
    </datafile>    
  </ga_license_file>
""" %(data_filename1, '2810517858', data_filename2, '2972536556')

	licfid.write(xml_string)
	licfid.close()

	licfid = open(license_filename)
        license_file_is_valid(licfid, data_filename1)
        license_file_is_valid(licfid, data_filename2)        
        licfid.close()
	    	
        # Clean up
        os.remove(license_filename)
        os.remove(data_filename1)
        os.remove(data_filename2)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_data_audit, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

