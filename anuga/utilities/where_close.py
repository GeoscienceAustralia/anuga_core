#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Single-function module.

   See function docstring for description.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: where_close.py,v 1.2 2004/04/28 00:28:15 jlin Exp $
#
# Modification History:
# - 19 Mar 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed reasonable tests.
#
# Notes:
# - Written for Python 2.2.2.
# - Function is based on code from the MA module by Paul F. Dubois.
#   Some snippets of code in this function are copied directly from 
#   lines in that module.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python where_close.py".
# - See import statements throughout for packages/modules required.
#
# Copyright (c) 2004 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/gemath/doc/;

# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA.

# You can contact Johnny Lin at his email address or at the University
# of Chicago, Department of the Geophysical Sciences, 5734 S. Ellis
# Ave., Chicago, IL 60637, USA.

#=======================================================================



#--------------------------- General Function --------------------------

def where_close(x, y, rtol=1.e-5, atol=1.e-8):
    """Mask of where x and y are element-wise "equal" to each other.

    Returns a long integer array with elements equal to 1 where x and 
    y are "equal", and 0 otherwise.  If x or y are floating point, 
    "equal" means where abs(x-y) <= atol + rtol * abs(y).  This is 
    essentially the same algorithm used in the numeric function 
    allclose.  If x and y are integer, "equal" means strict equality.  
    Shape and size of output is the same as x and y; if one is an 
    array and the other is scalar, shape and size of the output is the 
    same as the array.  Output is a numeric array, unless both inputs 
    are scalar in which the output is a Python integer scalar.

    Positional Input Arguments:
    * x:  Scalar or numeric array, Python list/tuple of any size and 
      shape.  Floating or integer type.
    * y:  Scalar or numeric array, Python list/tuple of any size and 
      shape.  Floating or integer type.

    Keyword Input Arguments:
    * rtol:   "Relative" tolerance.  Default is 1.e-5.  Used in the
              comparison between x and y only if the two are floating 
              point.
    * atol:   "Absolute" tolerance.  Default is 1.e-8.  Used in the
              comparison between x and y only if the two are floating 
              point.

    Examples:
    >>> import numpy as N
    >>> from where_close import where_close
    >>> x = [20.,  -32., -1., 2.            , 5., 29.]
    >>> y = [20.1, -31., -1., 2.000000000001, 3., 28.99]
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['0', '0', '1', '1', '0', '0']

    >>> x = N.array([1,  5,  7, -2, 10])
    >>> y = N.array([1, -5, 17, -2,  0])
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['1', '0', '0', '1', '0']
    """
    import numpy as N
    abs = N.absolute


    #- Make sure input is numeric type:

    xN = N.array(x)
    yN = N.array(y)


    #- Safe compare if floating.  Strict compare if integer.  Any other
    #  type returns an error:

    if (xN.dtype.char in N.typecodes['Float']) or \
       (yN.dtype.char in N.typecodes['Float']):
        return N.less_equal(abs(xN-yN), atol+rtol*abs(yN))

    elif (xN.dtype.char in N.typecodes['Integer']) and \
         (yN.dtype.char in N.typecodes['Integer']):
        return N.equal(xN, yN)

    else:
        raise ValueError("where_close:  Inputs must be Float or Integer")




#-------------------------- Main:  Test Module -------------------------

#- Define additional examples for doctest to use:

__test__ = { 'Additional Examples':
    """
    >>> from where_close import where_close
    >>> import numpy as N
    >>> x = [20.,  -32., -1., 2.            , 5., 29.]
    >>> y = [20.1, -31., -1., 2.000000000001, 3., 28.99]
    >>> x = N.reshape(x, (2,3))
    >>> y = N.reshape(y, (2,3))
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[0,i] for i in range(3)]
    ['0', '0', '1']
    >>> ['%.1g' % ind[1,i] for i in range(3)]
    ['1', '0', '0']
    >>> ind.shape
    (2, 3)
    >>> ind.dtype.char
    '?'
    >>> type(ind)
    <type 'numpy.ndarray'>

    >>> x = [20.,  -32., -1., 2.            , 5., 29.]
    >>> y = [20.1, -31., -1., 2.000000000001, 3.]
    >>> ind = where_close(x, y)
    Traceback (most recent call last):
        ...
    ValueError: frames are not aligned

    >>> x = [20, -32, -1]
    >>> y = [20, -32]
    >>> ind = where_close(x, y)
    Traceback (most recent call last):
        ...
    ValueError: frames are not aligned

    >>> x = [20,  -32, -1.0, 2, 5, 29]
    >>> y = 2.
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['0', '0', '0', '1', '0', '0']
    >>> x = -32
    >>> y = N.array([20.,  -32., -1., 2., 5., 29.])
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['0', '1', '0', '0', '0', '0']
    >>> x = -32
    >>> y = -33.
    >>> ind = where_close(x, y)
    >>> print ind
    0
    >>> type(ind)
    <type 'numpy.bool_'>
    >>> x = -33
    >>> y = -33.
    >>> ind = where_close(x, y)
    >>> print ind
    1
    >>> x = -33
    >>> y = -33
    >>> ind = where_close(x, y)
    >>> print ind
    1
    >>> x = -3.
    >>> y = -3.
    >>> ind = where_close(x, y)
    >>> print ind
    1
    """ }


#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module.

    Tests the examples in all the module documentation 
    strings, plus __test__.

    Note:  To help ensure that module testing of this file works, the
    parent directory to the current directory is added to sys.path.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====

