"""Dummy class for use with test_caching.py
"""


from builtins import str
from builtins import object
class Dummy:
    def __init__(self, value, another):
        self.value = value
        self.another = another
    
    def __repr__(self):
        return str(self.value) + ', ' + str(self.another)
        
    def __eq__(self, other):
        return (self.value == other.value and
                self.another == other.another)
    

# Define class Dummy_memorytest before any tests are run
# to make sure it has a different memory address
# to the one defined in test 'test_objects_are_created_memory'
class Dummy_memorytest:
    def __init__(self, value, another):
        self.value = value      
    
