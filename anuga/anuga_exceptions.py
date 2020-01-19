"""Exceptions used by ANUGA
"""

# Handle PendingDeprecationWarning causing an ImportError if using Python 3
try:
    from exceptions import IOError
    from exceptions import Exception
except ImportError:
    pass

class TitleError(IOError):
    """ Incorrect header in a file. """
    pass

class ParsingError(IOError):
    """ Could not parse a file. """
    pass

class ShapeError(IOError):
    """ Pathological shape in data. """
    pass

class ANUGAError(Exception):
    """ Generic ANUGA error. """
    #def __init__(self, args=None):
    #self.args = args
    pass

class DataMissingValuesError(Exception):
    """ Missing values in file. """
    pass

class DataFileNotOpenError(Exception):
    """ File is not open. """
    pass

class DataTimeError(Exception):
    """ Pathological time data. """
    pass

class DataDomainError(Exception):
    """ Pathological domain. """
    pass

class NewQuantity(Exception):
    """ Quantity used but not defined. """
    pass

class TitleValueError(Exception):
    """ Title of data column in file has wrong value. """
    pass
