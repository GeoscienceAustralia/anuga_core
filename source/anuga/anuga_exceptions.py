"""Exceptions used by ANUGA
"""

import exceptions

class TitleError(exceptions.IOError): pass
class ParsingError(exceptions.IOError): pass
class ShapeError(exceptions.IOError): pass

class ANUGAError(Exception):
    def __init__(self, args=None):
        self.args = args

class DataMissingValuesError(exceptions.Exception): pass
class DataFileNotOpenError(exceptions.Exception): pass
class DataTimeError(exceptions.Exception): pass
class DataDomainError(exceptions.Exception): pass
class NewQuantity(exceptions.Exception): pass
class TitleValueError(exceptions.Exception): pass
