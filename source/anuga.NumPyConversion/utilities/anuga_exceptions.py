"""Exceptions used by ANUGA
"""

import exceptions
class TitleError(exceptions.IOError): pass
class ParsingError(exceptions.IOError): pass
class ShapeError(exceptions.IOError): pass


class ANUGAError(Exception):
    def __init__(self, args=None):
        self.args = args
        
