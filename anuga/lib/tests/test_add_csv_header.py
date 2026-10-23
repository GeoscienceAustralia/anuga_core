"""Tests for anuga.lib.add_csv_header."""

import os
import tempfile
import unittest

from anuga.lib.add_csv_header import add_csv_header


class Test_add_csv_header(unittest.TestCase):

    def _write_csv(self, content):
        """Write content to a temp file, return its path."""
        fd, fname = tempfile.mkstemp(suffix='.csv')
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        return fname

    def tearDown(self):
        pass

    def test_add_header_normal(self):
        """Normal case: adds header line to a CSV file."""
        fname = self._write_csv('1,2,3\n4,5,6\n')
        try:
            add_csv_header(fname, ['a', 'b', 'c'])
            with open(fname) as f:
                lines = f.readlines()
            self.assertEqual(lines[0].strip(), 'a,b,c')
            self.assertEqual(lines[1].strip(), '1,2,3')
        finally:
            os.unlink(fname)

    def test_add_header_be_green(self):
        """be_green=True path: line-by-line reading."""
        fname = self._write_csv('1,2,3\n4,5,6\n')
        try:
            add_csv_header(fname, ['x', 'y', 'z'], be_green=True)
            with open(fname) as f:
                lines = f.readlines()
            self.assertEqual(lines[0].strip(), 'x,y,z')
            self.assertEqual(lines[1].strip(), '1,2,3')
        finally:
            os.unlink(fname)

    def test_missing_file_raises(self):
        """Non-existent file raises Exception."""
        with self.assertRaises(Exception):
            add_csv_header('/no/such/file.csv', ['a', 'b'])

    def test_empty_header_raises(self):
        """Empty header_list raises Exception."""
        fname = self._write_csv('1,2\n3,4\n')
        try:
            with self.assertRaises(Exception):
                add_csv_header(fname, [])
        finally:
            os.unlink(fname)

    def test_none_header_raises(self):
        """None header_list raises Exception."""
        fname = self._write_csv('1,2\n3,4\n')
        try:
            with self.assertRaises(Exception):
                add_csv_header(fname, None)
        finally:
            os.unlink(fname)

    def test_column_mismatch_raises(self):
        """Mismatched column count raises Exception (normal path)."""
        fname = self._write_csv('1,2,3\n4,5,6\n')
        try:
            with self.assertRaises(Exception):
                add_csv_header(fname, ['a', 'b'])  # file has 3 cols, header has 2
        finally:
            os.unlink(fname)

    def test_column_mismatch_be_green_raises(self):
        """Mismatched column count raises Exception (be_green path)."""
        fname = self._write_csv('1,2,3\n4,5,6\n')
        try:
            with self.assertRaises(Exception):
                add_csv_header(fname, ['a', 'b'], be_green=True)
        finally:
            os.unlink(fname)


if __name__ == '__main__':
    unittest.main()
