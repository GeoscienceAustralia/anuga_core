"""Fundamentals for computing the Mandelbrot set
"""

def calculate_point(c, kmax):
    """Calculate one point of the Mandelbrot set by iterating the 
    governing equation
    
        z_{k+1} = z_k^2 + c,  k = 0, 1, ... kmax-1
        z_0 = 0 + 0i        
   
    for as long as k < kmax  and |z_k| <= 2
    
    Inputs:
      c:    A complex number for which the iteration is computed
      kmax: The maximal number of iterations
      
    Output:
      count: The value of k after the iteration has completed. 
    """

    z = complex(0,0) #Create complex number with initial value
    k = 0            #Initialise iteration counter

    while k < kmax and abs(z) <= 2:
        z = z*z + c
        k = k+1

    return k



def calculate_region(real_min, real_max, imag_min, imag_max, kmax, M, N):
    """Calculate the mandelbrot set in a given rectangular subset of the complex plane.
    Inputs:
       real_min: Left boundary
       real_max: Right boundary
       imag_min: Lower boundary
       imag_max: Upper boundary       
       kmax: Maximal iteration count
       M: Number of points along the real axis
       N: Number of points along the imaginary axis
    Output:
       Matrix A (MxN) containing iteration counts for each point   
    """
    
    from Numeric import zeros
    from mandel_ext import calculate_point #Faster C-implementation

    #Compute resolution
    real_step = (real_max-real_min)/M  
    imag_step = (imag_max-imag_min)/N
    
    
    A = zeros((M, N))   # Create M x N matrix
    
    #Compute Mandelbrot iteration for each point in the rectangular subset
    #and store iteration counts in matrix A                   
    for i in range(M):
        for j in range(N):
            c = complex(real_min + i*real_step, imag_min + j*imag_step)
            A[i,j] = calculate_point(c, kmax)
            
    return A
