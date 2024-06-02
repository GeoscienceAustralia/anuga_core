

from builtins import object
from anuga.utilities.system_tools import log_to_file


class Operator(object):
    """Operator - generic structure for a fractional operator
    
    This is the base class for all fractional step operators
    """ 

    counter = 0

    def __init__(self,
                 domain,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        
        self.domain = domain
        self.domain.set_fractional_step_operator(self)


        # useful aliases
        self.stage_c = self.domain.quantities['stage'].centroid_values
        self.height_c = self.domain.quantities['height'].centroid_values
        self.xmom_c  = self.domain.quantities['xmomentum'].centroid_values
        self.ymom_c  = self.domain.quantities['ymomentum'].centroid_values
        self.elev_c  = self.domain.quantities['elevation'].centroid_values
        self.coord_c = self.domain.centroid_coordinates
        self.areas = self.domain.areas


        if domain.numproc > 1:
            msg = 'Not implemented to run in parallel'
            assert self.parallel_safe(), msg

        if description is None:
            self.description = ' '
        else:
            self.description = description


        self.set_label(label)


        self.verbose = verbose

        # Keep count of inlet operator
        Operator.counter += 1

        self.set_logging(logging)


    def __call__(self):

        #timestep = self.domain.get_timestep()
        raise Exception('Need to implement __call__ for your operator')
                    
    def get_timestep(self):

        return self.domain.get_timestep()


    def get_time(self):

        return self.domain.get_time()

    def parallel_safe(self):
        """By default an operator is not parallel safe
        """
        return False

    def statistics(self):

        message = 'You need to implement operator statistics for your operator'
        return message
    
    def timestepping_statistics(self):

        message  = 'You need to implement timestepping statistics for your operator'
        return message


    def print_statistics(self):

        print(self.statistics())

    def print_timestepping_statistics(self):

        print(self.timestepping_statistics())


    def log_timestepping_statistics(self):

        from anuga.utilities.system_tools import log_to_file
        if self.logging:
            log_to_file(self.log_filename, self.timestepping_statistics())

    def set_label(self, label=None):

        if label is None:
            self.label = "operator_%g" % Operator.counter
        else:
            self.label = label + '_%g' % Operator.counter


    def set_logging(self, flag=True):

        self.logging = flag
        


    def activate_logging(self):

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging:
            self.log_filename = self.label + '.log'
            log_to_file(self.log_filename, self.statistics(), mode='w')
            #log_to_file(self.log_filename, 'time,Q')

            #log_to_file(self.log_filename, self.culvert_type)




