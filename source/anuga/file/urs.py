##
# @brief 
class Read_urs:
    """
    Read the info in URS mux files.

    for the quantities here's a correlation between the file names and
    what they mean;
    z-mux is height above sea level, m
    e-mux is velocity is Eastern direction, m/s
    n-mux is velocity is Northern direction, m/s
    """

    ##
    # @brief Initialize this instance of Urs_points.
    # @param urs_file Path to the underlying data file.
    def __init__(self, urs_file):
        self.iterated = False
        columns = 3                         # long, lat , depth
        mux_file = open(urs_file, 'rb')

        # Number of points/stations
        (self.points_num,) = unpack('i', mux_file.read(4))

        # nt, int - Number of time steps
        (self.time_step_count,) = unpack('i', mux_file.read(4))
        #dt, float - time step, seconds
        (self.time_step,) = unpack('f', mux_file.read(4))
        msg = "Bad data in the urs file."
        if self.points_num < 0:
            mux_file.close()
            raise ANUGAError, msg
        if self.time_step_count < 0:
            mux_file.close()
            raise ANUGAError, msg
        if self.time_step < 0:
            mux_file.close()
            raise ANUGAError, msg

        # The depth is in meters, and it is the distance from the ocean
        # to the sea bottom.
        lonlatdep = p_array.array('f')
        lonlatdep.read(mux_file, columns * self.points_num)
        lonlatdep = num.array(lonlatdep, dtype=num.float)
        lonlatdep = num.reshape(lonlatdep, (self.points_num, columns))
        self.lonlatdep = lonlatdep

        self.mux_file = mux_file
        # check this array

    ##
    # @brief Allow iteration over quantity data wrt time.
    def __iter__(self):
        """
        iterate over quantity data which is with respect to time.

        Note: You can only iterate once over an object

        returns quantity infomation for each time slice
        """

        msg =  "You can only interate once over a urs file."
        assert not self.iterated, msg

        self.iter_time_step = 0
        self.iterated = True

        return self

    ##
    # @brief 
    def next(self):
        if self.time_step_count == self.iter_time_step:
            self.close()
            raise StopIteration

        #Read in a time slice from mux file
        hz_p_array = p_array.array('f')
        hz_p_array.read(self.mux_file, self.points_num)
        hz_p = num.array(hz_p_array, dtype=num.float)
        self.iter_time_step += 1

        return hz_p

    ##
    # @brief Close the mux file.
    def close(self):
        self.mux_file.close()
   
