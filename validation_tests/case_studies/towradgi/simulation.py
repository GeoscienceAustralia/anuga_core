""" 
Towradgi Creek 17 August 1998 Storm Event Calibration
By Petar Milevski, some revisions by Gareth Davies

Setting up a simulation class
"""

#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------


import anuga
import time
import numpy
import os





from anuga import myid, barrier, finalize, numprocs
from anuga import sequential_distribute_load
from anuga import sequential_distribute_dump
from anuga_parallel.sequential_distribute import sequential_distribute_load_pickle_file
from os.path import join

#from project import *



class Simulation(object):
    
    def __init__(self, outname = 'domain', partition_dir = 'PARTITIONS', verbose=True):
        self.verbose = verbose
        self.outname = outname
        self.partition_dir = partition_dir
    
    def setup_original_domain(self, np=None):
        """
        Create sequential domain and partition
        
        """
        verbose = self.verbose
        outname = self.outname
        partition_dir = self.partition_dir
        if np is None:
            np = numprocs
        
        # Only create the sequential domain on processor 0
        if myid == 0:
            
            if verbose: print 'CREATING PARTITIONED DOMAIN'
            
            # Let's see if we have already pickled this domain
            pickle_name = outname+'_P%g_%g.pickle'% (1,0)
            pickle_name = join(partition_dir,pickle_name)
    
            if os.path.exists(pickle_name):
                if verbose: print 'Saved domain seems to already exist'
            else:
                from setup_domain import setup_domain
                domain = setup_domain(verbose=self.verbose)
                
                if verbose: print 'Saving Domain'
                sequential_distribute_dump(domain, 1, partition_dir=partition_dir, verbose=verbose)    

            
            # Now partition the domain

            
            if verbose: print 'Load in saved sequential pickled domain'
            domain = sequential_distribute_load_pickle_file(pickle_name, np=1, verbose = verbose)
         
            par_pickle_name = outname+'_P%g_%g.pickle'% (np,0)
            par_pickle_name = join(partition_dir,par_pickle_name)
            if os.path.exists(par_pickle_name):
                if verbose: print 'Saved partitioned domain seems to already exist'
            else:
                if verbose: print 'Dump partitioned domains'
                sequential_distribute_dump(domain, np, partition_dir=partition_dir, verbose=verbose) 

        barrier()
        
        #===============================================================================
        # Setup parallel domain
        #===============================================================================
        if myid == 0 and verbose: print 'LOADING PARTITIONED DOMAIN'
        
        self.domain = sequential_distribute_load(filename=join(partition_dir,outname), verbose = self.verbose)
        #print self.domain.get_name()
        
        
    def setup_rainfall(self):
        """
        Create rainfall functions associated with polygons
        """
        
        if myid == 0 and self.verbose: print 'CREATING RAINFALL FUNCTIONS'
        from setup_rainfall import setup_rainfall
        
        setup_rainfall(self.domain) 

    def setup_structures(self):
        """
        Create structures such as culverts and bridges
        """
        
        if myid == 0 and self.verbose: print 'CREATING STRUCTURES'
        from setup_structures import setup_structures
        
        setup_structures(self.domain)
   
    def setup_boundaries(self):
        """
        Setup Boundary Conditions
        """
        
        if myid == 0 and self.verbose: print 'SETUP BOUNDARY CONDITIONS'
        from setup_boundaries import setup_boundaries
        
        setup_boundaries(self.domain)
        




 




