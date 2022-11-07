"""Classes for implementing damage curves and calculating financial damage

   Duncan Gray, Ole Nielsen, Jane Sexton, Nick Bartzis
   Geoscience Australia, 2006
"""

import os
from math import sqrt
from scipy.interpolate import interp1d
scipy_available = True


from random import choice

import numpy as num


try:  
    import kinds  
except ImportError:  
    # Hand-built mockup of the things we need from the kinds package, since it
    # was recently removed from the standard numeric distro.  Some users may  
    # not have it by default.  
    class _bunch(object):  
        pass  
         
    class _kinds(_bunch):  
        default_float_kind = _bunch()  
        default_float_kind.MIN = 2.2250738585072014e-308  #smallest +ve number
        default_float_kind.MAX = 1.7976931348623157e+308  
     
    kinds = _kinds()
    

from anuga.utilities.numerical_tools import ensure_numeric
from .exposure import Exposure
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.utilities.numerical_tools import NAN
from anuga.config import epsilon
import anuga.utilities.log as log

depth_epsilon = epsilon

# Change these if the ouput from nexix changes
SHORE_DIST_LABEL = 'SHORE_DIST'
WALL_TYPE_LABEL = 'WALL_TYPE'
STR_VALUE_LABEL = 'STR_VALUE'
CONT_VALUE_LABEL = 'CONT_VALUE'

def inundation_damage(sww_base_name, exposure_files_in,
                      exposure_file_out_marker=None,
                      ground_floor_height=0.3,
                      overwrite=False, verbose=True,
                                 use_cache = True):
    """
    This is the main function for calculating tsunami damage due to
    inundation.  It gets the location of structures from the exposure
    file and gets the inundation of these structures from the
    sww file.

    It then calculates the damage loss.

    Note, structures outside of the sww file get the minimum inundation
    (-ground_floor_height).
    
    These calculations are done over all the sww files with the sww_base_name
    in the specified directory.

    exposure_files_in - a file or a list of files to input from
    exposure_file_out_marker -  this string will be added to the input file
                                name to get the output file name
    """
    if isinstance(exposure_files_in, str):
        exposure_files_in = [exposure_files_in]


    for exposure_file_in in exposure_files_in:
        csv = Exposure(exposure_file_in,
                           title_check_list=[SHORE_DIST_LABEL,WALL_TYPE_LABEL,
                                             STR_VALUE_LABEL,CONT_VALUE_LABEL])
        geospatial = csv.get_location()
        geospatial = ensure_absolute(geospatial)
        max_depths, max_momentums = calc_max_depth_and_momentum(sww_base_name,
                        geospatial,
                        ground_floor_height=ground_floor_height,
                        verbose=verbose,
                        use_cache=use_cache)
        edm = EventDamageModel(max_depths,
                               csv.get_column(SHORE_DIST_LABEL),
                               csv.get_column(WALL_TYPE_LABEL),
                               csv.get_column(STR_VALUE_LABEL),
                               csv.get_column(CONT_VALUE_LABEL)
                               )
        results_dic = edm.calc_damage_and_costs(verbose_csv=True,
                                                verbose=verbose)
        for title, value in results_dic.items():
            csv.set_column(title, value, overwrite=overwrite)
    
        # Save info back to csv file
        if exposure_file_out_marker is None:
            exposure_file_out = exposure_file_in
        else:
            # split off extension, in such a way to deal with more than one '.' in the name of file
            split_name = exposure_file_in.split('.')
            exposure_file_out =  '.'.join(split_name[:-1]) + exposure_file_out_marker + \
                                '.' + split_name[-1]
        csv.save(exposure_file_out)
        if verbose: log.critical('Augmented building file written to %s'
                                 % exposure_file_out)
    
def add_depth_and_momentum2csv(sww_base_name, exposure_file_in,
                      exposure_file_out=None,
                      overwrite=False, verbose=True,
                                 use_cache = True):
    """
    Calculate the maximum depth and momemtum in an sww file, for locations
    specified in an csv exposure file.
    
    These calculations are done over all the sww files with the sww_base_name
    in the specified directory.
    """

    csv = Exposure(exposure_file_in)
    geospatial = csv.get_location()
    max_depths, max_momentums = calc_max_depth_and_momentum(sww_base_name,
                                                          geospatial,
                                                          verbose=verbose,
                                                          use_cache=use_cache)
    csv.set_column("MAX INUNDATION DEPTH (m)",max_depths, overwrite=overwrite)
    csv.set_column("MOMENTUM (m^2/s) ",max_momentums, overwrite=overwrite)
    csv.save(exposure_file_out)
    
def calc_max_depth_and_momentum(sww_base_name, points,
                                ground_floor_height=0.0,
                                verbose=True,
                                 use_cache = True):
    """
    Calculate the maximum inundation height above ground floor for a list
    of locations.

    The inundation value is in the range -ground_floor_height to
    overflow errors.

    These calculations are done over all the sww files with the sww_base_name
    in the specified directory.
    """

    quantities =  ['stage', 'elevation', 'xmomentum', 'ymomentum']
    points = ensure_absolute(points)
    point_count = len(points)

    # initialise the max lists
    max_depths = [-ground_floor_height]*point_count
    max_momentums = [-ground_floor_height]*point_count
    
    # How many sww files are there?
    dir, base = os.path.split(sww_base_name)
    if base[-4:] == '.sww':
        base = base[:-4]
    if dir == "": dir = "." # Unix compatibility
    dir_ls = os.listdir(dir)
    interate_over = [x for x in dir_ls if base in x and x[-4:] == '.sww']
    if len(interate_over) == 0:
        msg = 'No files of the base name %s.'\
              %(sww_base_name)
        raise IOError(msg)
    from os import sep

    for this_sww_file in interate_over:
        callable_sww = file_function(dir+sep+this_sww_file,
                                     quantities=quantities,
                                     interpolation_points=points,
                                     verbose=verbose,
                                     use_cache=use_cache)

        for point_i, point in enumerate(points):
            for time in callable_sww.get_time():
                quantity_values = callable_sww(time,point_i)
                w = quantity_values[0]
                z = quantity_values[1]
                uh = quantity_values[2] 
                vh = quantity_values[3]

                #print w,z,uh,vh
                if w == NAN or z == NAN or uh == NAN or vh == NAN:
                    continue
                    
                #  -ground_floor_height is the minimum value.
                depth = w - z - ground_floor_height
              
                if depth > max_depths[point_i]:
                    max_depths[point_i] = depth
                
                momentum = sqrt(uh*uh + vh*vh)
                if momentum > max_momentums[point_i]:
                    max_momentums[point_i] = momentum


    return max_depths, max_momentums

class EventDamageModel(object):
    """
    Object for working out the damage and cost

    """
    STRUCT_LOSS_TITLE = "STRUCT_LOSS_$"#"Structure Loss ($)"
    CONTENTS_LOSS_TITLE = "CONTENTS_LOSS_$"#"Contents Loss ($)"
    CONTENTS_DAMAGE_TITLE = "CONTENTS_DAMAGE_fraction"#"Contents damaged (fraction)"
    STRUCT_DAMAGE_TITLE = "STRUCT_DAMAGE_fraction" #"Structure damaged (fraction)"
    COLLAPSE_CSV_INFO_TITLE = "COLLAPSE_CSV_INFO"#"Calculation notes"
    MAX_DEPTH_TITLE = "MAX_DEPTH_m" #"Inundation height above ground floor (m)"
    STRUCT_COLLAPSED_TITLE = "STRUCT_COLLAPSED"#"collapsed structure if 1"
    STRUCT_INUNDATED_TITLE = "STRUCT_INUNDATED"#"inundated structure if 1"
    double_brick_damage_array = num.array([#[-kinds.default_float_kind.MAX, 0.0],
                                           [-1000.0, 0.0],
                                           [0.0-depth_epsilon, 0.0],
                                           [0.0,0.016],
                                           [0.1,0.150],
                                           [0.3,0.425],
                                           [0.5,0.449],
                                           [1.0,0.572],
                                           [1.5,0.582],
                                           [2.0,0.587],
                                           [2.5,0.647],
                                           [1000.0, 64.7]
                                           #[kinds.default_float_kind.MAX,64.7]
                                           ])

    if scipy_available:
        double_brick_damage_curve = interp1d(double_brick_damage_array[:,0],double_brick_damage_array[:,1])
    else:
        double_brick_damage_curve = InterpolatingFunction( \
             (num.ravel(double_brick_damage_array[:,0:1]),),
              num.ravel(double_brick_damage_array[:,1:]))
    
    brick_veeer_damage_array = num.array([#[-kinds.default_float_kind.MAX, 0.0],
                                          [-1000.0,0.0],
                                          [0.0-depth_epsilon, 0.0],
                                          [0.0,0.016],
                                          [0.1,0.169],
                                          [0.3,0.445],
                                          [0.5,0.472],
                                          [1.0,0.618],
                                          [1.5,0.629],
                                          [2.0,0.633],
                                          [2.5,0.694],
                                          [1000.0,69.4]
                                          #[kinds.default_float_kind.MAX,69.4]
                                          ])

    if scipy_available:
        brick_veeer_damage_curve = interp1d(brick_veeer_damage_array[:,0],brick_veeer_damage_array[:,1])
    else:
        brick_veeer_damage_curve = InterpolatingFunction( \
                                 (num.ravel(brick_veeer_damage_array[:,0:1]),),
                                  num.ravel(brick_veeer_damage_array[:,1:]))

    

    struct_damage_curve = {'Double Brick':double_brick_damage_curve,
                           'Brick Veneer':brick_veeer_damage_curve}
    default_struct_damage_curve = brick_veeer_damage_curve

    contents_damage_array = num.array([#[-kinds.default_float_kind.MAX, 0.0],
                                       [-1000.0,0.0],
                                       [0.0-depth_epsilon, 0.0],
                                       [0.0,0.013],
                                       [0.1,0.102],
                                       [0.3,0.381],
                                       [0.5,0.500],
                                       [1.0,0.970],
                                       [1.5,0.976],
                                       [2.0,0.986],
                                       [1000.0,98.6]
                                       #[kinds.default_float_kind.MAX,98.6]
                                       ])


    if scipy_available:
        contents_damage_curve = interp1d(contents_damage_array[:,0],contents_damage_array[:,1])
    else:
        contents_damage_curve = InterpolatingFunction( \
             (num.ravel(contents_damage_array[:,0:1]),),
              num.ravel(contents_damage_array[:,1:]))



    #building collapse probability
    # inundation depth above ground floor, m
    depth_upper_limits = [depth_epsilon, 1.0, 2.0, 3.0, 5.0,
                          kinds.default_float_kind.MAX]
    # shore mistance, m
    shore_upper_limits = [125,200,250, kinds.default_float_kind.MAX]
    # Building collapse probability
    collapse_probability = [[0.0, 0.0, 0.0, 0.0], #Code below assumes 0.0
                            [0.05, 0.02, 0.01, 0.0],
                            [0.6, 0.3, 0.1, 0.05],
                            [0.8, 0.4, 0.25, 0.15],
                            [0.95, 0.7, 0.5, 0.3],
                            [0.99, 0.9, 0.65, 0.45]]
    
    def __init__(self, max_depths, shore_distances, walls,
                 struct_costs, content_costs):
        """
        max depth is Inundation height above ground floor (m), so
                  the ground floor has been taken into account.
        """
        self.max_depths = [float(x) for x in max_depths]
        self.shore_distances = [float(x) for x in shore_distances]
        self.walls = walls
        self.struct_costs = [float(x) for x in struct_costs]
        self.content_costs = [float(x) for x in content_costs]

        self.structure_count = len(self.max_depths)
        #Fixme expand
        assert self.structure_count == len(self.shore_distances)
        assert  self.structure_count == len(self.walls)
        assert  self.structure_count == len(self.struct_costs)
        assert  self.structure_count == len(self.content_costs)
        #assert  self.structure_count == len(self.)

    def calc_damage_and_costs(self, verbose_csv=False, verbose=False):
        """
        This is an overall method to calculate the % damage and collapsed
        structures and then the $ loss.
        """
        self.calc_damage_percentages()
        collapse_probability = self.calc_collapse_probability()
        self._calc_collapse_structures(collapse_probability,
                                  verbose_csv=verbose_csv)
        self.calc_cost()
        results_dict = {self.STRUCT_LOSS_TITLE:self.struct_loss
                        ,self.STRUCT_DAMAGE_TITLE:self.struct_damage
                        ,self.CONTENTS_LOSS_TITLE:self.contents_loss
                        ,self.CONTENTS_DAMAGE_TITLE:self.contents_damage
                        ,self.MAX_DEPTH_TITLE:self.max_depths
                        ,self.STRUCT_COLLAPSED_TITLE:self.struct_collapsed
                        ,self.STRUCT_INUNDATED_TITLE:self.struct_inundated
                        }
        if verbose_csv:
           results_dict[self.COLLAPSE_CSV_INFO_TITLE] = self.collapse_csv_info
        return results_dict
            
    def calc_damage_percentages(self):
        """
        Using stage curves calc the damage to structures and contents
        """

        # the data being created
        struct_damage = num.zeros(self.structure_count, float)
        contents_damage = num.zeros(self.structure_count, float)
        self.struct_inundated = ['']* self.structure_count

        for i,max_depth,shore_distance,wall in zip(
                                                   list(range(self.structure_count)),
                                                   self.max_depths,
                                                   self.shore_distances,
                                                   self.walls):
            ## WARNING SKIP IF DEPTH < 0.0 
            if 0.0 > max_depth:
                continue
            
            # The definition of inundated is if the max_depth is > 0.0
            self.struct_inundated[i] = 1.0
            
            #calc structural damage %
            damage_curve = self.struct_damage_curve.get(wall,
                                              self.default_struct_damage_curve)
            struct_damage[i] = damage_curve(max_depth)
            contents_damage[i] = self.contents_damage_curve(max_depth)
           
        self.struct_damage = struct_damage
        self.contents_damage = contents_damage
           
        
    def calc_cost(self):
        """
        Once the damage has been calculated, determine the $ cost.
        """
        # ensure_numeric does not cut it.
        self.struct_loss = self.struct_damage * \
                           ensure_numeric(self.struct_costs)
        self.contents_loss = self.contents_damage * \
                           ensure_numeric(self.content_costs)
        
    def calc_collapse_probability(self):
        """
        return a dict of which structures have x probability of collapse.
             key is collapse probability
             value is list of struct indexes with key probability of collapse 
        """
        # I could've done this is the calc_damage_percentages and
        # Just had one loop.
        # But for ease of testing and bug finding I'm seperating the loops.
        # I'm make the outer loop for both of them the same though,
        # so this loop can easily be folded into the other loop.
        
        # dict of which structures have x probability of collapse.
        # key of collapse probability
        # value of list of struct indexes 
        struct_coll_prob = {}
        
        for i,max_depth,shore_distance,wall in zip(
                                                   list(range(self.structure_count)),
                                                   self.max_depths,
                                                   self.shore_distances,
                                                   self.walls):
            # WARNING ASSUMING THE FIRST BIN OF DEPTHS GIVE A ZERO PROBABILITY
            depth_upper_limits = self.depth_upper_limits
            shore_upper_limits = self.shore_upper_limits
            collapse_probability = self.collapse_probability
            if max_depth <= depth_upper_limits[0]:
                continue
            start = 1
            for i_depth, depth_limit in enumerate(depth_upper_limits[start:]):
                #Have to change i_depth so it indexes into the lists correctly
                i_depth += start
                if max_depth <= depth_limit:
                    for i_shore, shore_limit in enumerate(shore_upper_limits):
                        if shore_distance <= shore_limit:
                            coll_prob = collapse_probability[i_depth][i_shore]
                            if 0.0 == collapse_probability[i_depth][i_shore]:
                                break
                            struct_coll_prob.setdefault(coll_prob,[]).append(i)
                            break
                    break
                            
        return struct_coll_prob
    
    def _calc_collapse_structures(self, collapse_probability,
                                  verbose_csv=False):
        """
        Given the collapse probabilities, throw the dice
        and collapse some houses
        """
        
        self.struct_collapsed = [''] * self.structure_count
        if verbose_csv:
            self.collapse_csv_info = [''] * self.structure_count
        #for a given 'bin', work out how many houses will collapse
        for probability, house_indexes in collapse_probability.items():
            collapse_count = round(len(house_indexes) * probability)
            
            if verbose_csv:
                for i in house_indexes:
                    # This could be sped up I think
                    self.collapse_csv_info[i] = str(probability) + ' prob.( ' \
                           + str(int(collapse_count)) + ' collapsed out of ' \
                           + str(len(house_indexes)) + ')'
            for _ in range(int(collapse_count)):
                house_index = choice(house_indexes)
                self.struct_damage[house_index] = 1.0
                self.contents_damage[house_index] = 1.0
                house_indexes.remove(house_index)
                self.struct_collapsed[house_index] = 1
            
            # Warning, the collapse_probability list now lists 
            # houses that did not collapse, (though not all of them)
            
#############################################################################
if __name__ == "__main__":
    pass 
