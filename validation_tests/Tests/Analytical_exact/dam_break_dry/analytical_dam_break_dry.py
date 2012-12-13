#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       untitled.py
#       
#       Copyright 2011 Stephen Roberts <steve@babbage>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

    
def wu(h2, h0=1.0, h1=10.0):
    """ The intermediate state in dambreak problem
    must satisfy wu(h2) = 0.0
    
    Taken from formula (3.14) from Zoppou and Roberts SWW notes
    """

    import math

    h21 = h2/h1
    h01 = h0/h1

    #print h0,h1,h01
    

   # Calculate
    #h21**3 -9.0*(h21**2)*h01 +16.0*(h21**(1.5))*h01 -h21*h01*(h01+8.0) +h01**3



    sqrth21 = math.sqrt(h21)

    #calc_1 = h21**3 -9.0*(h21**2)*h01 +16.0*(h21**(1.5))*h01 -h21*h01*(h01+8.0) +h01**3

    calc_2 = h21*(sqrth21*(sqrth21*(h21-9*h01)+16.0*h01)- h01*(h01+8.0)) + h01**3


    #print calc_1, calc_2
    
    return calc_2


def calc_h2(h0=1.0, h1=10.0):

    assert h1 >  0.0
    assert h0 <=  h1
    assert h0 >= 0.0


    h01 = h0/h1

    if h01 >= 1 - 1.0e-4:
        return 0.5*(h0+h1)


    from scipy.optimize import brentq
    return brentq(wu, h0, h1, args=(h0,h1))





def dam_break(x,t,h0=1.0,h1 = 10.0):
    import math
    from anuga import g

    h2 = calc_h2(h0,h1)
    u2 = 2.0*(math.sqrt(g*h1) - math.sqrt(g*h2))
    #u2 = 2.0*(h1-h2)*g/(math.sqrt(g*h1) + math.sqrt(g*h2))

    try:
        s = u2*h2/(h2 - h0)
    except ZeroDivisionError:
        s = math.sqrt(g*h2)

    c1 = math.sqrt(g*h1)
    c2 = math.sqrt(g*h2)

    # use numpy.select
    if x < -t*c1:
        h = h1
        u = 0.0
    elif x < t*(u2 - c2):
        h = 1.0/g*(2.0/3.0 *c1 - 1.0/3.0*x/t)**2
        u = 2.0/3.0*(c1 + x/t)
    elif x < s*t:
        h = h2
        u = u2
    else:
        h = h0
        u = 0

    return h,u






def vec_dam_break(x,t,h0=1.0,h1 = 10.0):
    import math
    import numpy
    from anuga import g

    msg = 'Argument x should be a numpy array'
    assert isinstance(x,numpy.ndarray), msg

    h2 = calc_h2(h0,h1)
    u2 = 2.0*(math.sqrt(g*h1) - math.sqrt(g*h2))

    try:
        s = u2*h2/(h2 - h0)
    except ZeroDivisionError:
        s = math.sqrt(g*h2)

    c1 = math.sqrt(g*h1)
    c2 = math.sqrt(g*h2)


    condlist = [ x < -t*c1, x < t*(u2-c2), x < s*t]
    hchoicelist = [h1, 1.0/g*(2.0/3.0 *c1 - 1.0/3.0*x/t)**2, h2 ]
    uchoicelist = [0.0, 2.0/3.0*(c1 + x/t) , u2 ]


    h = numpy.select(condlist,hchoicelist,default=h0)
    u = numpy.select(condlist,uchoicelist,default=0.0)

    return h,u



if __name__ == "__main__":
    import numpy
    import pylab

    h1 = 10.0
    h0 = 1.0
    
    n = 100
    x = -1.0 + numpy.arange(0,2*n+1)/float(n)
    t = 0.075

    try:
        h,u = vec_dam_break(2.0,t,h0=h0)
    except:
        pass
    else:
        raise Exception('should have raised an exception')
    
    h,u = vec_dam_break(x,t,h0=h0)

    m = n/4
    hstar, ustar = dam_break(x[m],t,h0=h0)

    assert hstar == h[m]
    assert ustar == u[m]

    pylab.plot(x,h)
    #pylab.plot(x,u)

    pylab.show()
