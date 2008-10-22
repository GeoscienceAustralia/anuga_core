"""Dummy class for use with test_caching.py
"""


class Dummy:
    def __init__(self, value, another):
        self.value = value
        self.another = another
    
    def __repr__(self):
        return str(self.value) + ', ' + str(self.another)
    

# Define class Dummy_memorytest before any tests are run
# to make sure it has a different memory address
# to the one defined in test 'test_objects_are_created_memory'
class Dummy_memorytest:
    def __init__(self, value, another):
        self.value = value      
    
