#!/usr/bin/env python


import unittest
from tempfile import mkstemp, mktemp

import os

from anuga.utilities.xml_tools import *


class Test_xml_tools(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_generate_xml(self):
        """Test that xml code is generated from XMLobject model
        """

        X1 = XML_element(tag='first element',
                         value=XML_element(tag='some text',
                                           value='hello world'))

        text_elements = []
        for i in range(10):
            X = XML_element(tag='title %d' % i,
                            value='example text %d' % i)
            text_elements.append(X)

        X2 = XML_element(tag='second element',
                         value=XML_element(tag='texts',
                                           value=text_elements))
        X3 = XML_element(tag='third element',
                         value='42')

        doc = XML_element(value=[X1, X2, X3])

        # print doc.pretty_print()
        # print doc

        assert doc['second element']['texts']['title 4'] == 'example text 4'

        assert doc.has_key('first element')

    def test_xml2object(self):
        """Test that XML_document can be generated from file
        """

        tmp_name = mktemp(suffix='.xml')
        fid = open(tmp_name, 'w')

        xml_string = """<?xml version="1.0" encoding="iso-8859-1"?>
  <ga_license_file>
    <metadata>
      <author>Ole Nielsen</author>
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
""" % ('1234')

        fid.write(xml_string)
        fid.close()

        fid = open(tmp_name)
        reference = fid.read()
        reflines = reference.split('\n')

        xmlobject = xml2object(fid, verbose=True)

        # print xmlobject.pretty_print()

        xmllines = str(xmlobject).split('\n')

        # for line in reflines:
        #    print line
        # print
        # for line in xmllines:
        #    print line

        assert len(reflines) == len(xmllines)

        for i, refline in enumerate(reflines):
            msg = '%s != %s' % (refline.strip(), xmllines[i].strip())
            assert refline.strip() == xmllines[i].strip(), msg

        # Check dictionary behaviour
        for tag in list(xmlobject['ga_license_file'].keys()):
            xmlobject['ga_license_file'][tag]

        assert xmlobject['ga_license_file']['datafile']['accountable'] == 'Jane Sexton'

        # print
        # print xmlobject['ga_license_file']['datafile']
        # print xmlobject['ga_license_file']['metadata']
        # print xmlobject['ga_license_file']['datafile']
        # print xmlobject['ga_license_file']['datafile']['accountable']
        # print xmlobject['ga_license_file']['datafile'].keys()

        # for tag in xmlobject['ga_license_file'].keys():
        #    print xmlobject['ga_license_file'][tag]

        # Clean up
        fid.close()
        os.remove(tmp_name)

    def test_xml2object_empty_fields(self):
        """Test that XML_document can be generated from file
        This on tests that empty fields are treated as ''
        """

        tmp_name = mktemp(suffix='.xml')
        fid = open(tmp_name, 'w')

        xml_string = """<?xml version="1.0" encoding="iso-8859-1"?>
  <ga_license_file>
    <metadata>
      <author>Ole Nielsen</author>
    </metadata>
    <datafile>
      <filename>bathymetry.asc</filename>
      <checksum>%s</checksum>
      <publishable></publishable>
      <accountable>Jane Sexton</accountable>
      <source>Unknown</source>
      <IP_owner>Geoscience Australia</IP_owner>
      <IP_info>This is a test</IP_info>
    </datafile>
  </ga_license_file>
""" % ('1234')

        fid.write(xml_string)
        fid.close()

        fid = open(tmp_name)
        reference = fid.read()
        reflines = reference.split('\n')

        xmlobject = xml2object(fid, verbose=True)

        # print xmlobject.pretty_print()

        xmllines = str(xmlobject).split('\n')

        # for line in reflines:
        #    print line
        # print
        # for line in xmllines:
        #    print line

        ##assert len(reflines) == len(xmllines)
        x = xmlobject['ga_license_file']['datafile']['publishable']
        msg = 'Got %s, should have been an empty string' % x
        assert x == '', msg

        for i, refline in enumerate(reflines):
            msg = '%s != %s' % (refline.strip(), xmllines[i].strip())
            assert refline.strip() == xmllines[i].strip(), msg

        # Check dictionary behaviour
        for tag in list(xmlobject['ga_license_file'].keys()):
            xmlobject['ga_license_file'][tag]

        assert xmlobject['ga_license_file']['datafile']['accountable'] == 'Jane Sexton'

        # print
        # print xmlobject['ga_license_file']['datafile']
        # print xmlobject['ga_license_file']['metadata']
        # print xmlobject['ga_license_file']['datafile']
        # print xmlobject['ga_license_file']['datafile']['accountable']
        # print xmlobject['ga_license_file']['datafile'].keys()

        # for tag in xmlobject['ga_license_file'].keys():
        #    print xmlobject['ga_license_file'][tag]

        # Clean up
        fid.close()
        os.remove(tmp_name)

    def test_generate_and_read_back(self):
        """Test that xml code generated from XMLobject model
        can be read back.
        """

        X1 = XML_element(tag='first_element',
                         value=XML_element(tag='some_text',
                                           value='hello world'))

        text_elements = []
        for i in range(10):
            X = XML_element(tag='title_%d' % i,
                            value='example text %d' % i)
            text_elements.append(X)

        X2 = XML_element(tag='second_element',
                         value=XML_element(tag='texts',
                                           value=text_elements))
        X3 = XML_element(tag='third_element',
                         value='42')

        # Need to have one main element according to minidom
        main = XML_element(tag='all', value=[X1, X2, X3])
        xmldoc = XML_element(value=main)
        # print xmldoc

        tmp_name = mktemp(suffix='.xml')
        fid = open(tmp_name, 'w')
        fid.write(str(xmldoc))
        fid.close()

        # Now read it back
        xmlobject = xml2object(tmp_name, verbose=True)

        assert str(xmldoc) == str(xmlobject)

        os.remove(tmp_name)

    def test_duplicate_tags(self):
        """Test handling of duplicate tags.
        """

        X1 = XML_element(tag='datafile',
                         value=XML_element(tag='some_text',
                                           value='hello world'))

        X2 = XML_element(tag='second_element',
                         value=XML_element(tag='texts',
                                           value='egg and spam'))
        X3 = XML_element(tag='datafile',
                         value='42')

        # Need to have one main element according to minidom
        main = XML_element(tag='all', value=[X1, X2, X3])
        xmldoc = XML_element(value=main)
        # print xmldoc

        tmp_name = mktemp(suffix='.xml')
        fid = open(tmp_name, 'w')
        fid.write(str(xmldoc))
        fid.close()

        # Now read it back
        xmlobject = xml2object(tmp_name, verbose=True)
        # print xmlobject

        assert str(xmldoc) == str(xmlobject)

        assert xmlobject['all'].has_key('datafile')

        assert len(xmlobject['all']['datafile']) == 2
        # print xmlobject['all']['datafile']

        os.remove(tmp_name)

class Test_xml_tools_extra(unittest.TestCase):
    """Cover previously uncovered lines in xml_tools.py."""

    def _make_simple_xml_file(self):
        """Helper: write a simple XML file and return its path."""
        import tempfile
        content = '<root><item>hello</item><item>world</item></root>'
        fd, path = tempfile.mkstemp(suffix='.xml')
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        return path

    def test_print_tree(self):
        """print_tree covers lines 8-15."""
        from io import StringIO
        path = self._make_simple_xml_file()
        try:
            with open(path) as f:
                doc = parse(f)
            import sys
            # Just call it — output goes to stdout
            print_tree(doc.childNodes[0])
        finally:
            os.unlink(path)

    def test_pretty_print_tree(self):
        """pretty_print_tree covers line 19."""
        path = self._make_simple_xml_file()
        try:
            with open(path) as f:
                doc = parse(f)
            pretty_print_tree(doc)
        finally:
            os.unlink(path)

    def test_get_elements(self):
        """get_elements covers lines 39-44."""
        path = self._make_simple_xml_file()
        try:
            with open(path) as f:
                doc = parse(f)
            elems = get_elements(doc.childNodes)
            self.assertGreater(len(elems), 0)
        finally:
            os.unlink(path)

    def test_get_text(self):
        """get_text covers lines 51-58."""
        path = self._make_simple_xml_file()
        try:
            with open(path) as f:
                doc = parse(f)
            root = doc.childNodes[0]
            text = get_text(root.childNodes)
            # text is concatenation of text nodes
            self.assertIsInstance(text, str)
        finally:
            os.unlink(path)

    def test_xml_element_add_radd_repr(self):
        """__add__, __radd__, __repr__ cover lines 105, 108, 111."""
        e = XML_element(tag='x', value='42')
        s1 = e + ' suffix'
        self.assertIn('42', s1)
        s2 = 'prefix ' + e
        self.assertIn('42', s2)
        self.assertIn('42', repr(e))

    def test_getitem_none(self):
        """__getitem__ returns None when key not found (line 158)."""
        e = XML_element(tag='root',
                        value=[XML_element(tag='a', value='1')])
        result = e['missing']
        self.assertIsNone(result)

    def test_getitem_multiple(self):
        """__getitem__ returns list when multiple matches (lines 161→exit)."""
        e = XML_element(tag='root',
                        value=[XML_element(tag='item', value='1'),
                               XML_element(tag='item', value='2')])
        result = e['item']
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_pretty_print_method(self):
        """XML_element.pretty_print covers lines 179-189."""
        e = XML_element(tag='root',
                        value=[XML_element(tag='child', value='hello')])
        s = e.pretty_print()
        self.assertIn('root', s)
        self.assertIn('child', s)

    def test_xml2object_parse_error(self):
        """xml2object with bad XML raises (lines 219-224)."""
        import tempfile
        fd, path = tempfile.mkstemp(suffix='.xml')
        with os.fdopen(fd, 'w') as f:
            f.write('not valid xml <<<')
        try:
            with self.assertRaises(Exception) as cm:
                xml2object(path)
            self.assertIn('could not be parsed', str(cm.exception))
        finally:
            os.unlink(path)


################################################################################


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_xml_tools)
    runner = unittest.TextTestRunner()
    runner.run(suite)
