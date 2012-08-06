"""Plotting routine for use with the mandelbrot set

   Ole Nielsen, SUT 2003
"""


def plot(A, kmax = None):
    """Plot matrix A as an RGB image using the Python Imaging Library and Tkinter

       A is converted to an RGB image using PIL and saved to disk
       Then it is displayed using PhotoImage

       A must be square, integer valued, two dimensional and non-negative

       If kmax is omitted it will be set to the max value of A
       
       Ole Nielsen, SUT 2003
    """
    
    #from Tkinter import Frame, Canvas, TOP, NW, PhotoImage
    from numpy import transpose
    from Image import new             # PIL    
    import time

    t0 = time.time()

    # User definable parameters

    imtype = 'ppm'      # Image format (e.g 'ppm', 'bmp', 'tiff')
    filename ='mandel'  # Base filename
    
    exponent = 0.998    # Exponent for morphing colors, good with kmax = 2**15
    rgbmax = 2**24      # Normalisation constant for RGB
    

    # Input check
    assert len(A.shape) == 2, 'Matrix must be 2 dimensional'
    assert A.shape[0] == A.shape[1], 'Matrix must be square'
    msg = 'A must contain integers, I got %c' %A.dtype.char
    assert A.dtype.char in 'iIbBhHl', msg
    assert min(A.flat)>=0, 'A must be non-negative'

    if kmax is None:
        kmax = max(A.flat)

    # Convert values from A into RGB values (0 to 255) in each band
    N = A.shape[0]
    A = transpose(A).astype('d') # Cast as double         

    im = new("RGB", A.shape)

    
    L = []        
    try:
        from mandelplot_ext import normalise_and_convert
        normalise_and_convert(A, L, kmax, rgbmax, exponent)
        
    except:
        print 'WARNING: Could not import C extension from mandelplot_ext'
                
    
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):    
                
                c = A[i,j]/kmax

                if c == 1: c = 0       #Map convergent point (kmax) to black (0)
                c = c**exponent        #Morph slightly
                
                c = int(c * rgbmax)    #Normalise to 256 levels per channel
                
                red   = c / 256 / 256
                green = (c / 256) % 256
                blue  = c % 256
    
                L.append( (red, green, blue) )
            
            
    
    # Save image to file        
    im.putdata(L)
    im.save(filename + '.' + imtype, imtype)
    print 'Computed plot in %.2f seconds: ' %(time.time()-t0)

    # Display image on screen
    #answer = raw_input('Show image [Y/N][Y]?')
    #if answer.lower() in ['n', 'no']:
    #   import sys
    #   sys.exit()

    # Try to display using a image viewer
    import os
    os.system('eog %s' %(filename + '.' + imtype))
    #os.system('xv %s' %(filename + '.' + imtype))    



 
