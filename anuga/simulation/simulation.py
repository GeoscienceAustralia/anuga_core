"""
Setting up a simulation class
"""


#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
from builtins import range
from builtins import object
import anuga
import time
import numpy
import os


from anuga import myid, barrier, finalize, numprocs
from anuga.parallel.sequential_distribute import sequential_distribute_load
from anuga.parallel.sequential_distribute import sequential_distribute_dump
from anuga.parallel.sequential_distribute import sequential_distribute_load_pickle_file
from os.path import join


class Simulation(object):

    def __init__(self,
                 argument_adder=None,
                 from_commandline=True,
                 setup_boundaries=None,
                 setup_domain=None,
                 setup_rainfall=None,
                 setup_structures=None,
                 **kwargs):

        args = parse_args_and_parameters(argument_adder, from_commandline, **kwargs)

        self.verbose = args.verbose
        self.outname = args.outname
        self.partition_dir = args.partition_dir
        self.checkpoint_dir = args.checkpoint_dir
        self.alg = args.alg
        self.args = args
        self.checkpointing = args.checkpointing
        self.checkpoint_time = args.checkpoint_time

        self.setup_boundaries = setup_boundaries
        self.setup_domain = setup_domain
        self.setup_rainfall = setup_rainfall
        self.setup_structures = setup_structures

        if self.checkpointing:
            # try to read in from checkpoint file
            from anuga import load_checkpoint_file
            try:
                if myid == 0 and self.verbose:
                    print('TRYING TO OPEN CHECKPOINT FILES')
                self.domain = load_checkpoint_file(domain_name = self.outname, checkpoint_dir = self.checkpoint_dir)
                if myid == 0 and self.verbose:
                    print('OPENNED CHECKPOINT FILE at time = {}'.format(self.domain.get_time()))
            except:
                self.initialize_simulation()

            self.domain.set_checkpointing(checkpoint_time = self.checkpoint_time)
        else:
            self.initialize_simulation()


    def initialize_simulation(self):

        if myid == 0 and self.verbose:
            print('INITIALIZE SIMULATION')

        self._setup_original_domain()
        self._setup_structures()
        self._setup_rainfall()
        self._setup_boundaries()


    def run(self, yieldstep = None, finaltime =None):

        if yieldstep is None:
            yieldstep = self.args.yieldstep

        if finaltime is None:
            finaltime = self.args.finaltime

        if myid == 0 and self.verbose:
            print('EVOLVE(yieldstep = {}, finaltime = {})'.format(yieldstep,finaltime))


        domain = self.domain
        #import time
        t0 = time.time()

        for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):#= 83700.):

            if myid == 0 and self.verbose:
                domain.write_time()

        barrier()
        for p in range(numprocs):
            if myid == p:
                print('Processor %g ' %myid)
                print('That took %.2f seconds' %(time.time()-t0))
                print('Communication time %.2f seconds'%domain.communication_time)
                print('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
                print('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
            else:
                pass

            barrier()

#         if myid == 0 and self.verbose:
#             print 'Number of processors %g ' %numprocs
#             print 'That took %.2f seconds' %(time.time()-t0)
#             print 'Communication time %.2f seconds'%domain.communication_time
#             print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
#             print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time


        # --------------------------------------------------
        # Merge the individual sww files into one file
        # But don't delete the sub domain sww files
        # --------------------------------------------------
        domain.sww_merge()

        
        finalize()

    def _setup_original_domain(self, np=None):
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

            if verbose: print('CREATING PARTITIONED DOMAIN')

            # Let's see if we have already pickled this domain
            pickle_name = outname+'_P%g_%g.pickle'% (1,0)
            pickle_name = join(partition_dir,pickle_name)

            if os.path.exists(pickle_name):
                if verbose: print('Saved domain seems to already exist')
            else:
                domain = self.setup_domain(self)

                if verbose: print('Saving Domain')
                sequential_distribute_dump(domain, 1, partition_dir=partition_dir, verbose=verbose)


            # Now partition the domain


            if verbose: print('Load in saved sequential pickled domain')
            domain = sequential_distribute_load_pickle_file(pickle_name, np=1, verbose = verbose)

            par_pickle_name = outname+'_P%g_%g.pickle'% (np,0)
            par_pickle_name = join(partition_dir,par_pickle_name)
            if os.path.exists(par_pickle_name):
                if verbose: print('Saved partitioned domain seems to already exist')
            else:
                if verbose: print('Dump partitioned domains')
                sequential_distribute_dump(domain, np, partition_dir=partition_dir, verbose=verbose)

        barrier()

        #===============================================================================
        # Setup parallel domain
        #===============================================================================
        if myid == 0 and verbose: print('LOADING PARTITIONED DOMAIN')

        self.domain = sequential_distribute_load(filename=join(partition_dir,outname), verbose = self.verbose)
        #print self.domain.get_name()


    def _setup_rainfall(self):
        """
        Create rainfall functions associated with polygons
        """

        if myid == 0 and self.verbose: print('CREATING RAINFALL FUNCTIONS')

        if self.setup_rainfall is None:
            pass
        else:
            self.setup_rainfall(self)

    def _setup_structures(self):
        """
        Create structures such as culverts and bridges
        """

        if myid == 0 and self.verbose: print('CREATING STRUCTURES')

        if self.setup_structures is None:
            pass
        else:
            self.setup_structures(self)

    def _setup_boundaries(self):
        """
        Setup Boundary Conditions
        """

        if myid == 0 and self.verbose: print('SETUP BOUNDARY CONDITIONS')

        if self.setup_boundaries is None:
            pass
        else:
            self.setup_boundaries(self)





def parse_args(argument_adder, from_commandline=False, **kwargs):
    """
    Return Namespace of arguments.

    If from_commandline, the arguments are parsed from sys.argv;
    otherwise, defaults are returned.

    argument_adder is a function which modifies an argparse.ArgumentParser
    to give it simulation-specific arguments.

    Optional kwargs allow overriding default argument values.

    """
    # set up default arguments common to all simulations
    parser = anuga.create_standard_parser()

    # add any benchmark-specific arguments
    argument_adder(parser)

    # override default argument values
    parser.set_defaults(**kwargs)

    if from_commandline:
        return parser.parse_args()
    else:
        return parser.parse_args([])


def parse_args_and_parameters(argument_adder=None, from_commandline=False, **kwargs):
    """
    Return Namespace of arguments.

    If from_commandline, the arguments are parsed from sys.argv;
    otherwise, defaults are returned.

    argument_adder is a function which modifies an argparse.ArgumentParser
    to give it simulation-specific arguments.

    Optional kwargs allow overriding default argument values.

    Combine with local variables from project.py
    """



    parser = anuga.create_standard_parser()

    try:
        import project

        # Make variables in project.py into a dictionary
        project_dict = vars(project)
        del project_dict['__builtins__']
        del project_dict['__doc__']
        del project_dict['__file__']
        del project_dict['__name__']
        del project_dict['__package__']
        del project_dict['join']
    except:
        project_dict = {}

    # add any benchmark-specific arguments
    if argument_adder:
        argument_adder(parser)

    # override default argument values from project.py
    parser.set_defaults(**project_dict)

    # override default argument values
    parser.set_defaults(**kwargs)

    if from_commandline:
        return parser.parse_args()
    else:
        return parser.parse_args([])
