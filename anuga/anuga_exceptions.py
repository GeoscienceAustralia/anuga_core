"""Exceptions used by ANUGA
"""

import exceptions

class TitleError(exceptions.IOError):
    """ Incorrect header in a file. """
    pass

class ParsingError(exceptions.IOError):
    """ Could not parse a file. """
    pass
    
class ShapeError(exceptions.IOError):
    """ Pathological shape in data. """
    pass

class ANUGAError(exceptions.Exception):
    """ Generic ANUGA error. """
    #def __init__(self, args=None):
    #self.args = args
    pass

class DataMissingValuesError(exceptions.Exception):
    """ Missing values in file. """
    pass
    
class DataFileNotOpenError(exceptions.Exception):
    """ File is not open. """
    pass
    
class DataTimeError(exceptions.Exception):
    """ Pathological time data. """
    pass
    
class DataDomainError(exceptions.Exception):
    """ Pathological domain. """
    pass
    
class NewQuantity(exceptions.Exception):
    """ Quantity used but not defined. """
    pass
    
class TitleValueError(exceptions.Exception):
    """ Title of data column in file has wrong value. """
    pass
