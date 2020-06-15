## module error
""" err(string).
    Prints 'string' and terminates program.

    Taken from the book Numerical Methods in Engineering with Python
    by J. Kiusalaas    
"""
# FIXME (Ole): Use exception handling instead of this
import sys
def err(string):
    print(string)
    raw_input('Press return to exit')
    sys.exit()
