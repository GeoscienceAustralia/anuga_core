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
#   $Id: interp.py,v 1.2 2004/03/23 04:28:16 jlin Exp $
#
# Modification History:
# - 22 Mar 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.2.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python interp.py".
# - See import statements throughout for non-"built-in" packages and
#   modules required.
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




#----------------------- Overall Module Imports ------------------------

#- Set module version to package version:




#------------------------ Non-Private Function -------------------------

from builtins import range
def interp(y, x, xinterp, missing=1e+20):
    """Simple linear interpolation for ordinate with missing values.


    Vectors x and y are the data describing a piecewise linear function.
    Function returns the interpolated values of the ordinate function 
    at abscissa values in xinterp.  Values of xinterp outside the range 
    of x are returned as missing.  Any elements in the output that uses
    missing values in y for the interpolation are set to missing.


    Positional Input Arguments:
    * y:  Ordinate values of data.  Rank 1 numeric vector.  Required.
      Can have missing values.  Floating or integer type.

    * x:  Abscissa values of data.  Rank 1 numeric vector.  Required.
      Can have no missing values.  Must be monotonically ascending.  
      Floating or integer type.

    * xinterp:  Abscissa values to calculate interpolated ordinate 
      values at.  Rank 1 numeric vector or numeric scalar.  Required.  
      Can have no missing values.  Can be in any order.  Floating or 
      integer type.


    Keyword Input Arguments:
    * missing:  If input has missing values, this is the missing value 
      value.  Scalar.  Floating or integer type.  Default is 1e+20.


    Output Result:
    * Interpolated ordinate values at xinterp.  Rank 1 numeric vector 
      of same length as xinterp (if xinterp is a numeric scalar, 
      output is also a numeric scalar).  Missing values are set to the 
      value of argument missing.  Type is Float, even if argument 
      missing and inputs are all integer.


    References:
    * Lin, J. W.-B.:  "Simple Interpolation."
      Python/CDAT for Earth Scientists: Tips and Examples.
      http://www.johnny-lin.com/cdat_tips/tips_math/interp.html


    Example with no missing values (gives same output as function
    arrayfns.interp):

    >>> from interp import interp
    >>> import numpy as N
    >>> x = N.array([1., 2., 3., 4., 5.])
    >>> y = N.array([3., 6., 2.,-5.,-3.])
    >>> xint = N.array([3.4, 2.3])
    >>> yint = interp(y, x, xint, missing=1e+20)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['-0.8', '4.8']

    Example with missing values:

    >>> x = N.array([1.,    2., 3.,  4.,  5.])
    >>> y = N.array([3., 1e+20, 2., -5., -3.])
    >>> xint = N.array([3.4, 2.3])
    >>> yint = interp(y, x, xint, missing=1e+20)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['-0.8', '1e+20']

    Example with values out of range of the data:

    >>> x = N.array([1.,   2.1, 3.,  4., 5.1])
    >>> y = N.array([3., 1e+20, 2., -5., -3.])
    >>> xint = N.array([3.4, -2.3, 6.])
    >>> yint = interp(y, x, xint, missing=1e+20)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['-0.8', '1e+20', '1e+20']
    """
    import arrayfns
    import numpy.ma as MA
    import numpy as N
    from .where_close import where_close


    #- Check inputs for possible errors:

    if (N.rank(y) != 1) or (N.rank(x) != 1):
        raise ValueError("interp:  Input(s) not a vector")
    if N.rank(xinterp) > 1:
        raise ValueError("interp:  xinterp not a vector or scalar")
    if x[-1] <= x[0]:
        raise ValueError("interp:  x not monotonically increasing")


    #- Establish constants and define xint, a rank 1 version of
    #  xinterp to be used for the rest of the function:

    if N.rank(xinterp) == 0:
        xint = N.reshape(xinterp, (1,))
    else:
        xint = xinterp

    num_xint = N.size(xint)


    #- Mask as missing values of xint that are outside of the range
    #  of x:

    yint_outrange_mask = N.logical_or( N.less(xint, x[0]) \
                                     , N.greater(xint, x[-1]) )


    #- Mask of elements with missing values in y, if there are any
    #  missing values in y.  If xint equals a value in x, missing 
    #  values mask for that xint is the same as the corresponding 
    #  value in x; and mask elements in xint which fall in an interval 
    #  (whose upper bound index is top_idx) where one of the endpoints 
    #  is missing:

    y_miss_mask    = where_close(y, missing)
    yint_miss_mask = N.zeros(num_xint)

    if MA.maximum(y_miss_mask) == 1:

        for i in range(num_xint):
            if yint_outrange_mask[i] == 0:
                x_eq_xint = where_close(x, xint[i])
                if MA.maximum(x_eq_xint) == 1:
                    yint_miss_mask[i] = y_miss_mask[N.nonzero(x_eq_xint)]
                else:
                    top_idx = N.nonzero(N.greater(x, xint[i]))[0]
                    yint_miss_mask[i] = y_miss_mask[top_idx] or \
                                        y_miss_mask[top_idx-1]


    #- Return interpolated values, set to missing values as 
    #  appropriate, and making a scalar if xinterp is a scalar:

    yint = arrayfns.interp(y, x, xint)
    N.putmask( yint, N.logical_or(yint_miss_mask, yint_outrange_mask) \
             , missing)
    if N.rank(xinterp) == 0:  yint = yint[0]

    return yint




#-------------------------- Main:  Test Module -------------------------

#- Define additional examples for doctest to use:

__test__ = {'Additional Examples':
    """
    (1) General error catching:

    >>> from interp import interp
    >>> import numpy as N
    >>> x = N.array([1.,    2., 3.,  4.,  5.,  6.])
    >>> y = N.array([3., 1e+20, 2., -5., -3., -4.])
    >>> x = N.reshape(x, (2,3))
    >>> y = N.reshape(y, (2,3))
    >>> xint = N.array([3.4, 2.3])
    >>> yint = interp(y, x, xint, missing=1e+20)
    Traceback (most recent call last):
        ...
    ValueError: interp:  Input(s) not a vector

    >>> x = N.array([1.,    2., 3.,  4.,  5.,  6.])
    >>> y = N.array([3., 1e+20, 2., -5., -3., -4.])
    >>> xint = N.array([[3.4, 2.3],[3.4, 2.3]])
    >>> yint = interp(y, x, xint, missing=1e+20)
    Traceback (most recent call last):
        ...
    ValueError: interp:  xinterp not a vector or scalar

    >>> x = N.array([1.,    2., 3.,  4.,  5.,  0.])
    >>> y = N.array([3., 1e+20, 2., -5., -3., -4.])
    >>> xint = N.array([3.4, 2.3])
    >>> yint = interp(y, x, xint, missing=1e+20)
    Traceback (most recent call last):
        ...
    ValueError: interp:  x not monotonically increasing

    >>> x = N.array([1.,    2., 3.,  4.,  5.,  6.])
    >>> y = N.array([3., None, 2., -5., -3., -4.])
    >>> xint = N.array([3.4, 2.3, 2., 5., 3., 1.])
    >>> yint = interp(y, x, xint, missing=None)
    Traceback (most recent call last):
        ...
    ValueError: where_close:  Inputs must be Float or Integer

    (2) Values right on the border of intervals:

    >>> x = N.array([1.,    2., 3.,  4.,  5.,  6.])
    >>> y = N.array([3., 1e+20, 2., -5., -3., -4.])
    >>> xint = N.array([3.4, 2.3, 2., 5., 3., 1.])
    >>> yint = interp(y, x, xint, missing=1e+20)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['-0.8', '1e+20', '1e+20', '-3', '2', '3']

    (3) Single element vector input:

    >>> yint = interp(y, x, N.array([6.]), missing=1e+20)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['-4']

    (4) Scalar xint:

    >>> x = N.array([1.,    2., 3.,  4.,  5.,  6.])
    >>> y = N.array([3., 1e+20, 2., -5., -3., -4.])
    >>> yint = interp(y, x, N.array(6.), missing=1e+20)
    >>> yint
    -4.0
    >>> N.rank(yint)
    0

    (5) Integer values:

    >>> x = N.arange(6)
    >>> y = N.arange(6)
    >>> xint = N.array([3.4, 2.3])
    >>> yint = interp(y, x, xint, missing=-9999999)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['3.4', '2.3']
    >>> yint.dtype.char
    'd'
    >>> x = N.arange(6)
    >>> y = N.arange(6)
    >>> xint = N.array([3, 2])
    >>> yint = interp(y, x, xint, missing=-9999999)
    >>> ['%.7g' % yint[i] for i in range(len(yint))]
    ['3', '2']
    >>> yint.dtype.char
    'd'
    """}


#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module.

    Tests the examples in all the module documentation strings, plus
    __test__.

    Note:  To help ensure that module testing of this file works, the
    parent directory to the current directory is added to sys.path.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====

