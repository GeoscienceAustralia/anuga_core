import exceptions
class VectorShapeError(exceptions.Exception): pass
class ConvergenceError(exceptions.Exception): pass

import numpy as num

import anuga.utilities.log as log

class Stats:

    def __init__(self):

        self.iter = None
        self.rTr = None
        self.dt = None

def conjugate_gradient(A, b, x0=None, imax=10000, tol=1.0e-8, atol=1.0e-14,
                        iprint=None, output_stats=False):
    """
    Try to solve linear equation Ax = b using
    conjugate gradient method

    If b is an array, solve it as if it was a set of vectors, solving each
    vector.
    """
    
    if x0 is None:
        x0 = num.zeros(b.shape, dtype=num.float)
    else:
        x0 = num.array(x0, dtype=num.float)

    b  = num.array(b, dtype=num.float)


    if len(b.shape) != 1:
       
        for i in range(b.shape[1]):
            x0[:, i], stats = _conjugate_gradient(A, b[:, i], x0[:, i],
                                           imax, tol, atol, iprint)
    else:
        x0 , stats = _conjugate_gradient(A, b, x0, imax, tol, atol, iprint)

    if output_stats:
        return x0, stats
    else:
        return x0

    
def _conjugate_gradient(A, b, x0, 
                        imax=10000, tol=1.0e-8, atol=1.0e-10, iprint=None):
    """
   Try to solve linear equation Ax = b using
   conjugate gradient method

   Input
   A: matrix or function which applies a matrix, assumed symmetric
      A can be either dense or sparse or a function
      (__mul__ just needs to be defined)
   b: right hand side
   x0: inital guess (default the 0 vector)
   imax: max number of iterations
   tol: tolerance used for residual

   Output
   x: approximate solution
   """

    b  = num.array(b, dtype=num.float)
    if len(b.shape) != 1:
        raise VectorShapeError, 'input vector should consist of only one column'

    if x0 is None:
        x0 = num.zeros(b.shape, dtype=num.float)
    else:
        x0 = num.array(x0, dtype=num.float)


    if iprint == None  or iprint == 0:
        iprint = imax

    dx = 0.0
    
    i = 1
    x = x0
    r = b - A * x
    d = r
    rTr = num.dot(r, r)
    rTr0 = rTr
    
    #FIXME Let the iterations stop if starting with a small residual
    while (i < imax and rTr > tol ** 2 * rTr0 and rTr > atol ** 2):
        q = A * d
        alpha = rTr / num.dot(d, q)
        xold = x
        x = x + alpha * d

        dx = num.linalg.norm(x-xold)
        #if dx < atol :
        #    break
            
        if i % 50:
            r = b - A * x
        else:
            r = r - alpha * q
        rTrOld = rTr
        rTr = num.dot(r, r)
        bt = rTr / rTrOld

        d = r + bt * d
        i = i + 1
        if i % iprint == 0:
            log.info('i = %g rTr = %15.8e dx = %15.8e' % (i, rTr, dx))

        if i == imax:
            log.warning('max number of iterations attained')
            msg = 'Conjugate gradient solver did not converge: rTr==%20.15e' % rTr
            raise ConvergenceError, msg

    #print x

    stats = Stats()
    stats.iter = i
    stats.rTr = rTr
    stats.dx = dx

    return x, stats

