"""Abstract parallel interface - suitable for sequential programs

Use pypar for parallism if installed.
Otherwise define a rudimentary interface for sequential execution.
"""

try:
    import pypar 
except:
    print 'WARNING: Could not import pypar - defining sequential interface'
    def size(): return 1
    def rank(): return 0

    def get_processor_name():
        import os
        try:
            hostname = os.environ['HOST']
        except:
            try:  
                hostname = os.environ['HOSTNAME']  
            except:
                hostname = 'Unknown'  

        return hostname
      
    def abort():
        import sys
        sys.exit()

    def finalize(): pass
  
    def barrier(): pass  

    def time():
        import time
        return time.time()

    #def send(*args, **kwargs):
    #    pass

    #def receive(*args, **kwargs):
    #    pass    

else:
    from pypar import *
