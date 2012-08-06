"""Fundamental routines for computing the Mandelbrot set

   Ole Nielsen, SUT 2003
"""

def balance(N, P, p):
    """Compute p'th interval when N is distributed over P bins.
    """

    from math import floor

    L = int(floor(float(N)/P))
    K = N - P*L
    if p < K:
        Nlo = p*L + p
        Nhi = Nlo + L + 1
    else:
        Nlo = p*L + K
        Nhi = Nlo + L

    return Nlo, Nhi


def calculate_point(c, kmax):
    """Python version for calculating on point of the set
       This is slow and for reference purposes only.
       Use version from mandel_ext instead.
    """

    z = complex(0,0)

    count = 0
    while count < kmax and abs(z) <= 2:
        z = z*z + c
        count += 1

    return count    


def calculate_region(real_min, real_max, imag_min, imag_max, kmax, M, N,
                     Mlo = 0, Mhi = None, Nlo = 0, Nhi = None):
    """Calculate the mandelbrot set in the given region with resolution M by N
       If Mlo, Mhi or Nlo, Nhi are specified computed only given subinterval.
    """

    from numpy import zeros
    from mandel_ext import calculate_point  #Fast C implementation

    if Mhi is None: Mhi = M
    if Nhi is None: Nhi = N    

    real_step = (real_max-real_min)/M
    imag_step = (imag_max-imag_min)/N

    A = zeros((M, N), dtype='i')   # Create M x N matrix

    for i in range(Mlo, Mhi):
        for j in range(Nlo, Nhi):
            c = complex(real_min + i*real_step, imag_min + j*imag_step)
            A[i,j] = calculate_point(c, kmax)

    return A



def calculate_region_cyclic(real_min, real_max, imag_min, imag_max, kmax,
                            M, N, p=0, P=1, row = 1):
    """Calculate rows p+nP, n in N of the mandelbrot set in the given region
    with resolution M by N

    This is the most efficient way of partitioning the work.
    """


    from numpy import zeros
    from mandel_ext import calculate_point  #Fast C implementation

    real_step = (real_max-real_min)/M
    imag_step = (imag_max-imag_min)/N

    A = zeros((M, N), dtype='i')   # Create M x N matrix

    if row:
        for i in range(M):
            if i%P == p:
                for j in range(N):
                    c = complex(real_min + i*real_step, imag_min + j*imag_step)
                    A[i,j] = calculate_point(c, kmax)
    else:
        for j in range(N):        
            if j%P == p:
                for i in range(M):
                    c = complex(real_min + i*real_step, imag_min + j*imag_step)
                    A[i,j] = calculate_point(c, kmax)
    return A

    
