"""
Collect stage max stage info


"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



import numpy as num
from anuga.geometry.polygon import inside_polygon
from anuga.operators.base_operator import Operator
from anuga import Quantity


class Collect_max_stage_operator(Operator):
    """
    Simple operator to collect the max stage during a run

    """

    def __init__(self,
                 domain,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        #------------------------------------------
        # Setup a quantity to store max_stage
        #------------------------------------------
        self.max_stage = Quantity(domain, name = 'max_stage', register=True)
        self.max_stage.set_values(-1.0e+100)

        #------------------------------------------
        # Aliases for stage quantity
        #------------------------------------------
        self.stage  = domain.quantities['stage']

        

    def __call__(self):
        """
        Calculate max_stage at each timestep
        """

        self.max_stage.maximum(self.stage)


    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Collect_max_stage operator'
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Collecting_max_stage'
        return message


    def save_centroid_data_to_csv(self, filename=None):

        self.max_stage.save_centroid_data_to_csv(filename)


    def plot_quantity(self, filename=None, show=True):

        self.max_stage.plot_quantity(filename=filename, show=show)




