import numpy as num
from .cg_ext import jacobi_precon_c
from .cg_ext import cg_solve_c_precon
from .cg_ext import cg_solve_c
from anuga.utilities.sparse import Sparse, Sparse_CSR
import anuga.utilities.log as log

class VectorShapeError(Exception):
    pass


class ConvergenceError(Exception):
    pass


class PreconditionerError(Exception):
    pass


# Setup for C conjugate gradient solver


class Stats(object):

    def __init__(self):

        self.iter = None
        self.rTr = None
        self.dx = None
        self.rTr0 = None
        self.x = None
        self.x0 = None

    def __str__(self):
        msg = ' iter %.5g rTr %.5g x %.5g dx %.5g rTr0 %.5g x0 %.5g' \
              % (self.iter, self.rTr, self.x, self.dx, self.rTr0, self.x0)
        return msg

# Note Padarn 26/11/12: This function has been modified to include an
# additional argument 'use_c_cg' solve, which instead of using the current
# python implementation calls a c implementation of the cg algorithm. This
# has not been tested when trying to perform the cg routine on multiple
# quantities, but should work. No stats are currently returned by the c
# function.
# Note Padarn 26/11/12: Further note that to use the c routine, the matrix
# A must currently be in the sparse_csr format implemented in anuga.util.sparse

# Test that matrix is in correct format if c routine is being called


def conjugate_gradient(A, b, x0=None, imax=10000, tol=1.0e-8, atol=1.0e-14,
                       iprint=None, output_stats=False, use_c_cg=False, precon='None'):
    """
    Try to solve linear equation Ax = b using
    conjugate gradient method

    If b is an array, solve it as if it was a set of vectors, solving each
    vector.
    """

    if use_c_cg:
        from anuga.utilities.sparse import Sparse_CSR
        msg = ('c implementation of conjugate gradient requires that matrix A\
                be of type %s') % (str(Sparse_CSR))
        assert isinstance(A, Sparse_CSR), msg

    if x0 is None:
        x0 = num.zeros(b.shape, dtype=float)
    else:
        x0 = num.array(x0, dtype=float)

    b = num.array(b, dtype=float)

    err = 0

    # preconditioner
    # Padarn Note: currently a fairly lazy implementation, needs fixing
    M = None
    if precon == 'Jacobi':

        M = num.zeros(b.shape[0])
        jacobi_precon_c(A, M)
        x0 = b.copy()

        if len(b.shape) != 1:

            for i in range(b.shape[1]):

                if not use_c_cg:
                    x0[:, i], stats = _conjugate_gradient_preconditioned(A, b[:, i], x0[:, i], M,
                                                                         imax, tol, atol, iprint, Type="Jacobi")
                else:
                    # need to copy into new array to ensure contiguous access
                    xnew = x0[:, i].copy()
                    err = cg_solve_c_precon(
                        A, xnew, b[:, i].copy(), imax, tol, atol, b.shape[1], M)
                    x0[:, i] = xnew
        else:

            if not use_c_cg:
                x0, stats = _conjugate_gradient_preconditioned(
                    A, b, x0, M, imax, tol, atol, iprint, Type="Jacobi")
            else:
                err = cg_solve_c_precon(A, x0, b, imax, tol, atol, 1, M)

    else:

        if len(b.shape) != 1:

            for i in range(b.shape[1]):

                if not use_c_cg:
                    x0[:, i], stats = _conjugate_gradient(A, b[:, i], x0[:, i],
                                                          imax, tol, atol, iprint)
                else:
                    # need to copy into new array to ensure contiguous access
                    xnew = x0[:, i].copy()
                    err = cg_solve_c(
                        A, xnew, b[:, i].copy(), imax, tol, atol, b.shape[1])
                    x0[:, i] = xnew
        else:

            if not use_c_cg:
                x0, stats = _conjugate_gradient(
                    A, b, x0, imax, tol, atol, iprint)
            else:
                x0 = b.copy()
                err = cg_solve_c(A, x0, b, imax, tol, atol, 1)

    if err == -1:

        log.warning('max number of iterations attained from c cg')
        msg = 'Conjugate gradient solver did not converge'
        raise ConvergenceError(msg)

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

    stats = Stats()

    b = num.array(b, dtype=float)
    if len(b.shape) != 1:
        raise VectorShapeError(
            'input vector should consist of only one column')

    if x0 is None:
        x0 = num.zeros(b.shape, dtype=float)
    else:
        x0 = num.array(x0, dtype=float)

    stats.x0 = num.linalg.norm(x0)

    if iprint is None or iprint == 0:
        iprint = imax

    dx = 0.0

    i = 1
    x = x0
    r = b - A * x
    d = r
    rTr = num.dot(r, r)
    rTr0 = rTr

    stats.rTr0 = rTr0

    # FIXME Let the iterations stop if starting with a small residual
    while (i < imax and rTr > tol ** 2 * rTr0 and rTr > atol ** 2):
        q = A * d
        alpha = rTr / num.dot(d, q)
        xold = x
        x = x + alpha * d

        dx = num.linalg.norm(x-xold)

        # if dx < atol :
        #    break

        # Padarn Note 26/11/12: This modification to the algorithm seems
        # unnecessary, but also seem to have been implemented incorrectly -
        # it was set to perform the more expensive r = b - A * x routine in
        # 49/50 iterations. Suggest this being either removed completely or
        # changed to 'if i%50==0' (or equvialent).
        # if i % 50:
        if False:
            r = b - A * x
        else:
            r = r - alpha * q
        rTrOld = rTr
        rTr = num.dot(r, r)
        bt =  rTr / rTrOld

        d = r + bt * d
        i = i + 1

        if i % iprint == 0:
            log.info('i = %g rTr = %15.8e dx = %15.8e' % (i, rTr, dx))

        if i == imax:
            log.warning('max number of iterations attained')
            msg = 'Conjugate gradient solver did not converge: rTr==%20.15e' % rTr
            raise ConvergenceError(msg)

    stats.x = num.linalg.norm(x)
    stats.iter = i
    stats.rTr = rTr
    stats.dx = dx

    return x, stats


def _conjugate_gradient_preconditioned(A, b, x0, M,
                                       imax=10000, tol=1.0e-8, atol=1.0e-10, iprint=None, Type='None'):
    """
   Try to solve linear equation Ax = b using
   preconditioned conjugate gradient method

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

    # Padarn note: This is temporary while the Jacboi preconditioner is the only
    # one avaliable.
    D = []
    if not Type == 'Jacobi':
        log.warning(
            'Only the Jacobi Preconditioner is impletment cg_solve python')
        msg = 'Only the Jacobi Preconditioner is impletment in cg_solve python'
        raise PreconditionerError(msg)
    else:
        D = Sparse(A.M, A.M)
        for i in range(A.M):
            D[i, i] = 1 / M[i]
        D = Sparse_CSR(D)

    stats = Stats()

    b = num.array(b, dtype=float)
    if len(b.shape) != 1:
        raise VectorShapeError(
            'input vector should consist of only one column')

    if x0 is None:
        x0 = num.zeros(b.shape, dtype=float)
    else:
        x0 = num.array(x0, dtype=float)

    stats.x0 = num.linalg.norm(x0)

    if iprint is None or iprint == 0:
        iprint = imax

    dx = 0.0

    i = 1
    x = x0
    r = b - A * x
    z = D * r
    d = r
    rTr = num.dot(r, z)
    rTr0 = rTr

    stats.rTr0 = rTr0

    # FIXME Let the iterations stop if starting with a small residual
    while (i < imax and rTr > tol ** 2 * rTr0 and rTr > atol ** 2):
        q = A * d
        alpha = rTr / num.dot(d, q)
        xold = x
        x = x + alpha * d

        dx = num.linalg.norm(x-xold)

        # if dx < atol :
        #    break

        # FIXME: Padarn Note 26/11/12: This modification to the algorithm seems
        # unnecessary, but also seem to have been implemented incorrectly -
        # it was set to perform the more expensive r = b - A * x routine in
        # 49/50 iterations. Suggest this being either removed completely or
        # changed to 'if i%50==0' (or equvialent).
        # if i % 50:
        if False:
            r = b - A * x
        else:
            r = r - alpha * q
        rTrOld = rTr
        z = D * r
        rTr = num.dot(r, z)
        bt = rTr / rTrOld

        d = z + bt * d
        i = i + 1
        if i % iprint == 0:
            log.info('i = %g rTr = %15.8e dx = %15.8e' % (i, rTr, dx))

        if i == imax:
            log.warning('max number of iterations attained')
            msg = 'Conjugate gradient solver did not converge: rTr==%20.15e' % rTr
            raise ConvergenceError(msg)

    stats.x = num.linalg.norm(x)
    stats.iter = i
    stats.rTr = rTr
    stats.dx = dx

    return x, stats
