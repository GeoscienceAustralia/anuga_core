#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$16/06/2010 4:29:59 PM$"



def  minmod(a,b):
    from numpy import abs, where

    phi = where((abs(a) < abs(b)) & (a*b >= 0.0), a, 0.0)
    phi = where((abs(b) < abs(a)) & (a*b >= 0.0), b, phi)

    return phi

def  minmod_kurganov(a,b,c):
    from numpy import sign, abs, minimum, where

    return where( (sign(a)*sign(b) >= 0.0) & (sign(a)*sign(c)>0.0),
        sign(a)*minimum(minimum(abs(a),abs(b)),abs(c)), 0.0 )



def  maxmod(a,b):

    from numpy import abs, where

    phi =  where((abs(a) > abs(b)) & (a*b > 0.0), a, 0.0)
    phi =  where((abs(b) > abs(a)) & (a*b > 0.0), b, phi)

    return phi


def vanleer(a,b):

    from numpy import fabs, where

    return where((fabs(a) + fabs(b) >= 1.0e-12), (a*fabs(b)+fabs(a)*b)/(fabs(a)+fabs(b)+ 1.0e-50), 0.0)


def vanalbada(a,b):

    from numpy import abs, where

    return where((a*a + b*b >= 1.0e-32), (a*a*b+a*b*b)/(a*a+b*b), 0.0)



def  minmod_old(beta_p,beta_m):
    if (abs(beta_p) < abs(beta_m)) & (beta_p*beta_m > 0.0):
        phi = beta_p
    elif (abs(beta_m) < abs(beta_p)) & (beta_p*beta_m > 0.0):
        phi = beta_m
    else:
        phi = 0.0
    return phi


def vanleer_old(a,b):
    if abs(a)+abs(b) > 1e-12:
        return (a*abs(b)+abs(a)*b)/(abs(a)+abs(b))
    else:
        return 0.0


def vanalbada_old(a,b):
    if a*a+b*b > 1e-12:
        return (a*a*b+a*b*b)/(a*a+b*b)
    else:
        return 0.0


def  maxmod_old(a,b):
    if (abs(a) > abs(b)) & (a*b > 0.0):
        phi = a
    elif (abs(b) > abs(a)) & (a*b > 0.0):
        phi = b
    else:
        phi = 0.0
    return phi

def  minmod_kurganov_old(a,b,c):
    from numpy import sign
    if sign(a)==sign(b)==sign(c):
        return sign(a)*min(abs(a),abs(b),abs(c))
    else:
        return 0.0
