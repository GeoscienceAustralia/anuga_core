#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose

from system_tools import *

class Test_system_tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_user_name(self):
        user = get_user_name()

        # print user
        assert isinstance(user, basestring), 'User name should be a string'

    def test_host_name(self):
        host = get_host_name()

        # print host
        assert isinstance(host, basestring), 'User name should be a string'        
        
        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)




