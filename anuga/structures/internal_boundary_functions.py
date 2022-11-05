
import numpy
import scipy
from scipy.interpolate import interp1d


class hecras_internal_boundary_function(object):

    """ Read internal boundary curves for a bridge / culvert from hecras and
        convert to a function which can be input as a structure in ANUGA

       If using hec-ras for only part of the structure (e.g. under the bridge),
        be careful to get the tables from an adjusted hecras model where only
        flow under the bridge is allowed (and ineffective flow areas etc are
        adjusted to reflect this). Also, the curves will change as manning's n
        is adjusted so you might want to make a number of curves with different
        n values to make calibration easier

        Note also that when used in an ANUGA structure, this can still give a
        different result to hecras. 
            * One reason is that the water surface elevation in hecras is
              assumed constant over a cross-section, whereas substantial
              cross-channel variations can occur in ANUGA, especialy near a
              bridge opening. The enquiry points used by ANUGA are a single
              location, so won't 'average-out' the cross-sectional variation
            * ANUGA doesn't include side-wall friction


    """

    def __init__(self, internal_boundary_curves_file, skip_header_rows=4,
                 skip_columns=1, allow_sign_reversal=False, verbose=True, 
                 vertical_datum_offset = 0.):
        """ Use a csv file containing the htab-curves from hecras to create an
            interpolation function for the structure. 

            The header / initial column probably need to be skipped (hecras
            format)

            It is assumed that the first 2 columns are the free-flow curve,
            and the others are Q vs headwater rating curves for various
            tailwater values

            @param internal_boundary_curves_file A csv file with the format
                    copied from hecras's internal boundary curves table for
                    the structure. This is obtained from (in hecras 5.00) the
                    menu via View-Hydraulic Property Tables, then click on
                    the table tab, copy the table, paste into a spreadsheet
                    and then save as csv.
            @param skip_header_rows Number of header rows in
                   internal_boundary_curves_file (4 in the examples I have seen)
            @param skip_columns Number of columns to ignore (on the
                   left-hand-side of the table), 1 in the examples I have seen
            @param allow_sign_reversal If False, then if the _call_ method is invoked
                    with the headwater < tailwater, it produces an exception.
                    If True, then in the same situation the code reverses the 
                    headwater and tailwater, and also reverses the sign of the
                    resulting discharge. 
            @param verbose True/False more messages

        """

        if verbose:
            print('########################################')
            print('HECRAS INTERNAL BOUNDARY FUNCTION')
            print('THIS IS EXPERIMENTAL')
            print('SUBJECT TO CHANGE WITHOUT NOTICE')
            print('########################################')


        internal_boundary_curves = numpy.genfromtxt(
            internal_boundary_curves_file, delimiter=',',
            skip_header=skip_header_rows)

        internal_boundary_curves = internal_boundary_curves[:, skip_columns:]

        # Do we use the table for flows from upstream to downstream, and the reverse?
        # If so, then reverse flow has a negative sign
        self.allow_sign_reversal = allow_sign_reversal

        # Adjust HW/TW curves by vertical datum offset
        for i in range(1, internal_boundary_curves.shape[1], 2):
            internal_boundary_curves[:,i] = internal_boundary_curves[:,i] +\
                 vertical_datum_offset

        # The first 2 columns consist of the free overflow curve (Q, HW). This
        # is apparently used when, for a given tail-water, the head-water is
        # not on the corresponding curve.
        free_flow_data = remove_repeated_curve_points(
            internal_boundary_curves[:, 1],
            internal_boundary_curves[:, 0])
        self.free_flow_curve = interp1d(
            free_flow_data[:, 0], free_flow_data[:, 1], kind='linear')
        self.free_flow_data = free_flow_data

        self.free_flow_hw_range = numpy.array(
            [internal_boundary_curves[:, 1].min(), 
             internal_boundary_curves[:, 1].max()])

        # Aside from the first 2 columns, the file consists of column-stacked
        # Q,HW curves defining the discharge Q when TW = HW[0]. Hence Q[0] = 0
        # always.
        #
        # We store these as a list of interpolation functions for each TW, with
        # other arrays which give the TW, and corresponding maximum HW
        #

        ncol_ibc = internal_boundary_curves.shape[1]
        self.nonfree_flow_tw = \
            internal_boundary_curves[0, list(range(3, ncol_ibc + 1, 2))]

        # Ensure the nonfree_flow_tw is monotonic increasing
        assert numpy.all(self.nonfree_flow_tw[0:-1] < self.nonfree_flow_tw[1:])

        self.nonfree_flow_tw_range = numpy.array(
            [self.nonfree_flow_tw.min(), self.nonfree_flow_tw.max()])

        # Populate in loop
        self.hw_max_given_tw = self.nonfree_flow_tw * numpy.nan
        curve_counter = 0
        self.nonfree_flow_curves = list()
        self.nonfree_flow_data = list()
        for i in range(2, internal_boundary_curves.shape[1], 2):
            Q = internal_boundary_curves[:, i]
            HW = internal_boundary_curves[:, i + 1]
            TW = HW[0]

            assert Q[0] == 0.

            HW_max = numpy.nanmax(HW)
            Q_max = numpy.nanmax(Q)
            self.hw_max_given_tw[curve_counter] = HW_max

            if(HW_max == HW[0]):
                # Curve with no range in HW
                # Can happen for the extreme curve
                if i == internal_boundary_curves.shape[1] - 2:
                    if verbose:
                        print('Skipping final rating curve with no HW range')

                    self.nonfree_flow_tw = self.nonfree_flow_tw[0:-1]
                    self.nonfree_flow_tw_range[1] = self.nonfree_flow_tw[-1]
                    self.hw_max_given_tw = self.hw_max_given_tw[0:-1]
                    continue
                else:
                    print(i, internal_boundary_curves.shape[1], HW_max)
                    raise Exception('Rating curve with no HW range')

            curve_data = remove_repeated_curve_points(HW, Q)

            self.nonfree_flow_data.append(curve_data)

            self.nonfree_flow_curves.append(
                interp1d(curve_data[:, 0], curve_data[:, 1], kind='linear'))
            curve_counter += 1

        self.internal_boundary_curves = internal_boundary_curves
        self.name = internal_boundary_curves_file

        return



    def __call__(self, hw_in, tw_in):
        """
            Compute the discharge from the headwater / tailwater

            The basic approach is:
                1) Identify the 2 HW-Q curves corresponding to TW levels just
                   above/below tw. Say TW_l is the lower TW level, and TW_u is
                   the upper one
                2) Evaluate Q_l, Q_u from each of the l and u curves,
                   interpolating at "(hw - tw) + TW_l" and "(hw-tw) + TW_u"
                   respectively
                    - For each curve, if the adjusted 'hw' value exceeds the hw
                      range of the curve, then use the free overflow value.
                    - NOTE: This means that we may have to evaluate the function
                            at values above hw -- which can lead to errors if the
                            table does not cover enough hw range.
                3) Take a weighted average of each Q, with weights 
                   inversely proportional to (tw - TW_l) and (TW_u-tw)

            This weighting method is nice since it preserves steady states.
            Note the upper/lower curves are evaluated with the 'correct head
            loss', rather than the correct hw value (but incorrect tw value).
            This seems a smoother approach.

            If self.allow_sign_reversal is True, then if hw < tw, we swap
            hw and tw, and return the resultant Q with sign reversed

            @param hw_in the input headwater
            @param tw_in the input tailwater

        """

        # Usually hw >= tw. If not, see if we should reverse the sign of Q
        if ((hw_in < tw_in) and self.allow_sign_reversal):
            tw = 1.0*hw_in
            hw = 1.0*tw_in
            sign_multiplier = -1.0
        else:
            tw = 1.0*tw_in
            hw = 1.0*hw_in
            sign_multiplier = 1.0

        # Logical checks

        if hw < tw:
            msg = 'HW: ' + str(hw) + ' < TW: ' + str(tw) + ' in ' + self.name
            raise Exception(msg)

        # Quick exit
        if hw < self.free_flow_hw_range[0]:
            return 0.0

        if hw > self.free_flow_hw_range[1]:
            msg = 'HW: ' + str(hw) + ' is outside free_flow_hw_range ' +\
                str(self.free_flow_hw_range) + ' in ' + self.name
            raise Exception(msg)

        if tw > self.nonfree_flow_tw.max():
            msg = 'TW: ' + str(tw) + ' exceeds nonfree_flow_tw_max ' +\
                str(self.nonfree_flow_tw.max()) + ' in ' + self.name
            raise Exception(msg)

        # Compute discharge
        if tw < self.nonfree_flow_tw[0]:
            # Use free flow curve
            # This could induce a discontinuity in the flow as
            # tw crosses the minimum
            Q = self.free_flow_curve(hw)

        else:
            # Try to use nonfree flow curves
            tw_lower_index = (self.nonfree_flow_tw <= tw).sum() - 1
            max_allowed_tw_index = len(self.nonfree_flow_tw) - 1

            # Interpolate along the 'lower-tw' curve, with hw adjusted
            # to get the head-loss right. This will preserve stationary states
            lower_tw = self.nonfree_flow_tw[tw_lower_index]
            lower_hw = (hw - tw) + lower_tw
            # Interpolation weight
            w0 = tw - lower_tw

            if tw_lower_index < max_allowed_tw_index:
                # Get upper curve variables
                tw_upper_index = tw_lower_index + 1
                upper_tw = self.nonfree_flow_tw[tw_upper_index]
                upper_hw = (hw - tw) + upper_tw
                w1 = upper_tw - tw
            else:
                # tw_lower_index corresponds to the highest curve
                # In that case just use the lower curve
                tw_upper_index = tw_lower_index
                upper_tw = lower_tw
                upper_hw = lower_hw
                w1 = 1.
                w0 = 0.

            assert ((w0 >= 0.) & (w1 >= 0.)), str(w0) + ' ' + str(w1)

            # Check whether lower_hw or upper_hw exceed the allowed rating
            # curve range
            max_possible_lower_hw = \
                max(self.hw_max_given_tw[tw_lower_index],
                    self.free_flow_hw_range[1])
            if lower_hw > max_possible_lower_hw:
                msg = 'lower_hw ' + str(lower_hw) + \
                    ' < max_possible_lower_hw: ' +\
                    str(max_possible_lower_hw) + \
                    ' at structure ' + self.name + ' \n'+\
                    ' This may exceed hw: ' + str(hw) +\
                    ' but is required for our interpolation method.' +\
                    ' Fix by extending the range of the rating curves'
                raise Exception(msg)

            max_possible_upper_hw = \
                max(self.hw_max_given_tw[tw_upper_index],
                    self.free_flow_hw_range[1])
            if(upper_hw > max_possible_upper_hw):
                msg = 'upper_hw ' + str(upper_hw) + \
                    ' < max_possible_upper_hw: ' +\
                    str(max_possible_upper_hw) +\
                    ' at structure ' + self.name + ' \n'+\
                    ' This may exceed hw: ' + str(hw) +\
                    ' but is required for our interpolation method.' +\
                    ' Fix by extending the range of the rating curves'

                raise Exception(msg)

            # Get lower curve value
            # NOTE: This could introduce a (probably small) discontinuity
            # in the interpolation if the tabulated curves do not exactly agree
            # (which seems to be the case for hecras output)
            if lower_hw <= self.hw_max_given_tw[tw_lower_index]:
                lower_curve_Q = self.nonfree_flow_curves[
                    tw_lower_index](lower_hw)
            else:
                # Use free flow curve. The 'real' Q is hw
                lower_curve_Q = self.free_flow_curve(hw)

            # Get upper curve value
            if upper_hw < self.hw_max_given_tw[tw_upper_index]:
                upper_curve_Q = self.nonfree_flow_curves[
                    tw_upper_index](upper_hw)
            else:
                # Use free flow curve.
                upper_curve_Q = self.free_flow_curve(hw)

            Q = (w0 * upper_curve_Q + w1 * lower_curve_Q)/ (w0 + w1)

        #print 'Q: ', Q , ' HW: ', hw, ' TW:', tw

        return(Q*sign_multiplier)


    def grid_function(self, interactive_plot=True):
        """ Compute Q for each valid HW / TW combination
            Optionally plot it.
            Return a list with [HW_values, TW_values, Q].

        """
        from matplotlib import pyplot

        HW_min = self.free_flow_hw_range[0]
        HW_max = self.free_flow_hw_range[1]
        HW_range = scipy.linspace(HW_min, HW_max, num=101)
        TW_range = scipy.linspace(HW_min, self.nonfree_flow_tw.max(), num=101)

        Q_store = numpy.array((101 ** 2) * [0.])
        Q_store = Q_store.reshape((101, 101))

        for h, hw in enumerate(HW_range):
            for t, tw in enumerate(TW_range):
                if tw > hw:
                    continue
                try:
                    Q_store[h, t] = self(hw, tw)
                except:
                    Q_store[h, t] = numpy.nan

        if interactive_plot:
            pyplot.ion()
            pyplot.contourf(HW_range, TW_range, Q_store, 100)
            pyplot.colorbar()
            pyplot.title('Q given headwater and tailwater: ' + self.name)
            pyplot.xlabel('Tailwater (m)')
            pyplot.ylabel('Headwater (m)')

        # Be nice to use the grid as a quick lookup function. 
        # But it fails when nan is present
        # self.gridded_interpolation_function = scipy.interpolate.RectBivariateSpline(HW_range, TW_range, Q_store)

        return [HW_range, TW_range, Q_store]


def remove_repeated_curve_points(HW, Q):
    """Utility function to get rid of round-off-induced points with repeated
    HW values from the HW-Q interpolation curves

    Return a 2-column array with HW, Q but repeated HW values removed. We
    retain the smallest Q for each HW, except at the highest HW value where the
    largest Q is used. 

    """

    HW_max = numpy.nanmax(HW)
    Q_max = numpy.nanmax(Q)

    points2keep = list()
    for j in range(Q.shape[0]):
        # Check for nan
        if Q[j] != Q[j]:
            continue

        # Because of rounding, sometimes a single curve contains
        # multiple HW
        # Exclude curve points which repeat earlier HW
        # Don't to this for the first or last points on the curve
        if ((j > 0) & (HW[j] == HW[j - 1]) & (HW[j] != HW_max)):
            continue

        if((HW[j] == HW_max) & (Q[j] < Q_max)):
            continue

        points2keep.append([HW[j], Q[j]])

    return numpy.array(points2keep)



class pumping_station_function(object):
    """ Transfer water from one site to another at a given rate,
        based on the pump capacities and observed headwater/tailwater.

        This class returns a callable object which can compute the rate
        of pumping. This can be passed to an internal_boundary_operator.

        The locations of the latter are defined based on the enquiry points
        of the structure. The positive pump direction (positive Q) is from
        headwater to tailwater. 
    """

    def __init__(self, domain, pump_capacity, hw_to_start_pumping, hw_to_stop_pumping,
                 initial_pump_rate=0., pump_rate_of_increase = 1.0e+100, 
                 pump_rate_of_decrease = 1.0e+100, verbose=True):
        """
            @param domain ANUGA domain
            @param pump_capacity (m^3/s)
            @param hw_to_start_pumping Turn pumps on if hw exceeds this (m)
            @param hw_to_stop_pumping Turn pumps off if hw below this (m)
            @param initial_pump_rate rate of pumps at start of simulation  (m^3/s)
            @param pump_rate_of_increase Accelleration of pump rate when turning on (m^3/s/s)
            @param pump_rate_of_decrease Decelleration of pump rate when turning off (m^3/s/s)
            @param verbose 
        """

        self.pump_capacity = pump_capacity
        self.hw_to_start_pumping = hw_to_start_pumping
        self.hw_to_stop_pumping = hw_to_stop_pumping

        self.pump_rate_of_increase = pump_rate_of_increase
        self.pump_rate_of_decrease = pump_rate_of_decrease

        self.domain=domain
        self.last_time_called = domain.get_time()
        self.time = domain.get_time()

        if hw_to_start_pumping < hw_to_stop_pumping:
            raise Exception('hw_to_start_pumping should be >= hw_to_stop_pumping')

        if initial_pump_rate > pump_capacity:
            raise Exception('Initial pump rate is > pump capacity')
        
        if ((self.pump_rate_of_increase < 0.) | (self.pump_rate_of_decrease < 0.)):
            raise Exception('Pump rates of increase / decrease MUST be non-negative')

        if ( (pump_capacity < 0.) | (initial_pump_rate < 0.) ):
            raise Exception('Pump rates cannot be negative')

        self.pump_rate = initial_pump_rate

        if verbose:
            print('########################################')
            print('PUMPING STATION FUNCTION')
            print('THIS IS EXPERIMENTAL')
            print('SUBJECT TO CHANGE WITHOUT NOTICE')
            print('########################################')


    def __call__(self, hw_in, tw_in):
        """
            @param hw_in The stage (or energy) at the headwater site (m)
            @param tw_in The stage (or energy) at the tailwater site (m)
        """

        # Compute the time since last called, so we can increase / decrease the pump rate if needed
        self.time = self.domain.get_time()
        if self.time > self.last_time_called:
            dt = self.time - self.last_time_called
            self.last_time_called = self.time
        else:
            dt = 0.
            if self.time != self.last_time_called:
                print(self.time)
                print(self.last_time_called)
                print(self.time - self.last_time_called)
                raise Exception('Impossible timestepping, ask Gareth')
          
        # Increase / decrease the pump rate if needed 
        if hw_in < self.hw_to_stop_pumping:
            self.pump_rate = max(0., self.pump_rate - dt*self.pump_rate_of_decrease)
        elif hw_in > self.hw_to_start_pumping:
            self.pump_rate = min(self.pump_capacity, self.pump_rate + dt*self.pump_rate_of_increase)

        return self.pump_rate

