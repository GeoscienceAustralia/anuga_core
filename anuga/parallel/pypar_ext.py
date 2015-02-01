# =============================================================================
# pypar_ext.py - Extension to Parallel Python using MPI
# Copyright (C) 2012 Stephen G Roberts
#              (ANU)
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#
#
# Contact address: Stephen.Roberts@anu.edu.au
#
# Version: See pypar_ext.__version__
# =============================================================================

"""Module pypar_ext.py -Extension to pypar

Need to import pypar first to init mpi.

Public functions:



isend() -- Asyncronous send (arrays)
receive() --  Asyncronous receive (arrays)
allreduce() -- wrapper for MPI_Allreduce (array)


See doc strings of individual functions for detailed documentation.
"""

from numpy import zeros, reshape, product
#from __metadata__ import __version__, __date__, __author__


# Constants
#
max_tag = 32767      # Max tag value (MPI_TAG_UB didn't work and returned 0)
control_tag = 13849  # Reserved tag used to identify control information
default_tag = 1      # Tag used as default if not specified

control_sep = ':'          # Separator for fields in control info (NOT ',')
control_data_max_size = 64 # Maximal size of string holding control data


#---------------------------------------------------------------------------
# Communication functions
#--------------------------------------------------------------------------

def isend(x, destination, use_buffer=False, vanilla=False,
         tag=default_tag, bypass=False):
    """Wrapper for easy MPI send.
       Send x to destination.

       Automatically determine appropriate protocol
       and call corresponding send function.
       Also passes type and size information on as preceding message to
       simplify the receive call.

       The variable x can be any (picklable) type, but
       numpy variables and text strings will most efficient.
       Setting vanilla = 1 forces vanilla mode for any type.

       If use_buffer is True, workspace x will be used for the return value.
       In this case the corresponding receive call must specify a buffer.
       Otherwise a new workspace will be created by receive.

       If bypass is True, all admin and error checks
       get bypassed to reduce the latency. Should only
       be used for sending numpy arrays and should be matched
       with a bypass in the corresponding receive command.

    """
    import types, string

    if bypass:
        isend_array(x, destination, tag)
        return

    # Input check
    errmsg = 'Destination id (%s) must be an integer.' % destination
    assert type(destination) == types.IntType, errmsg

    errmsg = 'Tag %d is reserved by pypar - please use another.' % control_tag
    assert tag != control_tag, errmsg

    # Create metadata about object to be sent
    control_info, x = create_control_info(x, vanilla, return_object=True)
    protocol = control_info[0]


    # Possibly transmit control data
    if use_buffer is False:
        send_control_info(control_info, destination)


    # Transmit payload data
    if protocol == 'array':
        send_array(x, destination, tag)
    elif protocol in ['string', 'vanilla']:
        send_string(x, destination, tag)
    else:
        raise 'Unknown protocol: %s' %protocol


def ireceive(source, buffer=None, vanilla=False, tag=default_tag,
            return_status=False, bypass=False):
    """receive - blocking MPI receive

       Receive data from source.

       Optional parameters:
         buffer: Use specified buffer for received data (faster). Default None.
         vanilla: Specify to enforce vanilla protocol for any type.
                  Default False

         tag: Only received messages tagged as specified. Default default_tag
         return_status: Return Status object along with result. Default False.

       If no buffer is specified, receive will try to receive a
       preceding message containing protocol, type, size and shape and
       then create a suitable buffer.

       If buffer is specified the corresponding send must specify
       use_buffer = True.
       The variable buffer can be any (picklable) type, but
       numpy variables and text strings will most efficient.

       Appropriate protocol will be automatically determined
       and corresponding receive function called.


       If bypass is True, all admin and error checks
       get bypassed to reduce the latency. Should only
       be used for receiving numpy arrays and should
       be matched with a bypass in the corresponding send command.
       Also buffer must be specified.
    """

    if bypass:
        # errmsg = 'bypass mode must be used with specified buffer'
        # assert buffer is not None, msg
        stat = ireceive_array(buffer, source, tag)
    else:

        import types

        # Input check
        errmsg = 'Source id (%s) must be an integer.' %source
        assert type(source) == types.IntType, errmsg

        errmsg = 'Tag %d is reserved by pypar - please use another.'\
            % control_tag
        assert tag != control_tag, errmsg


        # Either receive or create metadata about object to receive
        if buffer is None:
            control_info, source = receive_control_info(source,
                                                        return_source=True)
            protocol, typecode, size, shape = control_info
        else:
            protocol, typecode, size, shape = create_control_info(buffer,
                                                                  vanilla)


        # Receive payload data
        if protocol == 'array':
            if buffer is None:
                buffer = zeros(size,typecode)
                buffer = reshape(buffer, shape)

            stat = receive_array(buffer, source, tag)

        elif protocol == 'string':
            if buffer is None:
                buffer = ' '*size

            stat = receive_string(buffer, source, tag)

        elif protocol == 'vanilla':
            from cPickle import dumps, loads, UnpicklingError
            if buffer is None:
                s = ' '*size
            else:
                s = dumps(buffer, protocol=2)
                s = s + ' '*int(0.1*len(s)) #safety

            stat = receive_string(s, source, tag)
            try:
                buffer = loads(s)   #Replace buffer with received result
            except UnpicklingError, err:
                raise UnpicklingError(str(err) + " - '%s'" % s)
        else:
            raise 'Unknown protocol: %s' %protocol

    # Return received data and possibly the status object
    if return_status:
        return buffer, Status(stat)
    else:
        return buffer




def allreduce(x, op, buffer=None, vanilla=0, bypass=False):
    """Allreduce elements in x to buffer (of the same size as x)
       applying operation op elementwise.

       If bypass is True, all admin and error checks
       get bypassed to reduce the latency.
       The buffer must be specified explicitly in this case.
    """

    if bypass:
        allreduce_array(x, buffer, op)
        return


    import types
    from pypar import size
    numproc = size()         # Needed to determine buffer size



    # Create metadata about object
    protocol, typecode, size, shape = create_control_info(x)

    # Allreduce
    if protocol == 'array':
        if buffer is None:
            buffer = zeros(size*numproc, typecode)

            # Modify shape along axis=0 to match size
            shape = list(shape)
            shape[0] *= numproc
            buffer = reshape(buffer, shape)


        msg = 'Data array and buffer must have same type '
        msg = 'in allreduce. I got types "%s" and "%s"' % (x.dtype.char,
                                                        buffer.dtype.char)
        assert x.dtype.char == buffer.dtype.char, msg
        allreduce_array(x, buffer, op)


    elif (protocol == 'vanilla' or protocol == 'string'):
        raise 'Protocol: %s unsupported for allreduce' % protocol
    else:
        raise 'Unknown protocol: %s' % protocol

    return buffer




#---------------------------------------------------------
# INTERNAL FUNCTIONS
#---------------------------------------------------------

class Status:
    """ MPI Status block returned by receive if
        specified with parameter return_status=True
    """

    def __init__(self, status_tuple):
        self.source = status_tuple[0]  # Id of sender
        self.tag = status_tuple[1]     # Tag of received message
        self.error = status_tuple[2]   # MPI Error code
        self.length = status_tuple[3]  # Number of elements transmitted
        self.size = status_tuple[4]    # Size of one element

    def __repr__(self):
        return 'Pypar Status Object:\n  source=%d\n  tag=%d\n  error=%d\n  length=%d\n  size=%d\n' %(self.source, self.tag, self.error, self.length, self.size)

    def bytes(self):
        """Number of bytes transmitted (excl control info)
        """
        return self.length * self.size



def create_control_info(x, vanilla=0, return_object=False):
    """Determine which protocol to use for communication:
       (numpy) arrays, strings, or vanilla based x's type.

       There are three protocols:
       'array':   numpy arrays of type 'i', 'l', 'f', 'd', 'F' or 'D' can be
                  communicated  with mpiext.send_array and mpiext.receive_array.
       'string':  Text strings can be communicated with mpiext.send_string and
                  mpiext.receive_string.
       'vanilla': All other types can be communicated as string representations
                  provided that the objects
                  can be serialised using pickle (or cPickle).
                  The latter mode is less efficient than the
                  first two but it can handle general structures.

       Rules:
       If keyword argument vanilla == 1, vanilla is chosen regardless of
       x's type.
       Otherwise if x is a string, the string protocol is chosen
       If x is an array, the 'array' protocol is chosen provided that x has one
       of the admissible typecodes.

       The optional argument return_object asks to return object as well.
       This is useful in case it gets modified as in the case of general structures
       using the vanilla protocol.
    """

    import types

    # Default values
    protocol = 'vanilla'
    typecode = ' '
    size = 0
    shape = ()

    # Determine protocol in case
    if not vanilla:
        if type(x) == types.StringType:
            protocol = 'string'
            typecode = 'c'
            size = len(x)
        elif type(x).__name__ == 'ndarray': # numpy isn't imported yet
            try:
                import numpy
            except:
                print "WARNING (pypar.py): numpy module could not be imported,",
                print "reverting to vanilla mode"
                protocol = 'vanilla'
            else:
                typecode = x.dtype.char
                if typecode in ['i', 'l', 'f', 'd', 'F', 'D']:
                    protocol = 'array'
                    shape = x.shape
                    size = product(shape)
                else:
                    print "WARNING (pypar.py): numpy object type %s is not supported."\
                          %(x.dtype.char)
                    print "Only types 'i', 'l', 'f', 'd', 'F', 'D' are supported,",
                    print "Reverting to vanilla mode."
                    protocol = 'vanilla'

    # Pickle general structures using the vanilla protocol
    if protocol == 'vanilla':
        from cPickle import dumps
        x = dumps(x, protocol=2)
        size = len(x) # Let count be length of pickled object

    # Return
    if return_object:
        return [protocol, typecode, size, shape], x
    else:
        return [protocol, typecode, size, shape]



#----------------------------------------------

def send_control_info(control_info, destination):
    """Send control info to destination
    """
    import string

    # Convert to strings
    control_info = [str(c) for c in control_info]

    control_msg = string.join(control_info,control_sep)
    if len(control_msg) > control_data_max_size:
        errmsg = 'Length of control_info exceeds specified maximium (%d)'\
                 %control_data_max_size
        errmsg += ' - Please increase it (in pypar.py)'
        raise errmsg

    send_string(control_msg, destination, control_tag)


def receive_control_info(source, return_source=False):
    """Receive control info from source

    The optional argument (due to Jim Bosch) also returns the actual source node
    which can be used to require that the data message come from the same node.
    """

    # FIXME (Ole): Perhaps we should include actual source in the control info?

    import string

    msg = ' '*control_data_max_size

    stat = receive_string(msg, source, control_tag)
    # No need to create status object here - it is reserved
    # for payload communications only

    msg = msg[:stat[3]] # Trim buffer to actual received length (needed?)

    control_info = msg.split(control_sep)

    assert len(control_info) == 4, 'len(control_info) = %d' %len(control_info)
    control_info[2] = eval(control_info[2]) # Convert back to int
    control_info[3] = eval(control_info[3]) # Convert back to tuple

    if return_source:
        return control_info, int(stat[0])
    else:
        return control_info


#----------------------------------------------------------------------------
# Initialise module
#----------------------------------------------------------------------------


# Take care of situation where module is part of package
import sys, os, string, os.path
dirname = os.path.dirname(string.replace(__name__,'.',os.sep)).strip()

if not dirname:
    dirname = '.'

if dirname[-1] != os.sep:
    dirname += os.sep



# Import MPI extension
#
# Verify existence of mpiext.so.

try:
    import mpiextras
except:
    errmsg = 'ERROR: C extension mpiextras could not be imported.\n'
    errmsg += 'Please compile mpiextras.c manually.\n'
    #raise Exception, errmsg
    error = 1
    print errmsg
else:

    # Determine if MPI program is allowed to run sequentially on current platform
    # Attempting to check this automatically may case some systems to hang.

    if sys.platform in ['linux2', 'sunos5', 'win32', 'darwin']:
        # Linux (LAM,MPICH) or Sun (MPICH)
        error = 0  #Sequential execution of MPI is allowed
    else:
        # Platform: Alpha 'osf1V5'
        cmdstring = '"import mpiext, sys; mpiext.init(sys.argv); mpiext.finalize()"'
        #s = 'cd %s; python -c %s' %(dirname, cmdstring)
        s = 'python -c %s >/dev/null 2>/dev/null' %cmdstring
        error = os.system(s)

        # The check is performed in a separate shell.
        # Reason: The Alpha server, LAM/Linux or the Sun cannot recover from a
        # try:
        #   mpiext.init(sys.argv)

        # However, on LAM/Linux, this test causes system to hang.
        # Verified (OMN 12/12/2)
        # If lamboot is started, the system, will hang when init is called
        # again further down in this file.
        # If lamboot is not loaded error will be nozero as it should.
        # I don't know how to deal with this
        #
        #Comparisons of two strategies using LAM
        #
        # Strategy 1: Assume seq execution is OK (i.e. set error = 0)
        # Strategy 2: Try to test if mpi can be initialised (in a separate shell)
        #
        #
        # Strategy 1 (currently used)
        #                    | Lam booted  | Lam not booted
        #-----------------------------------------------------
        #
        # Sequential exec    |  OK         | Not OK
        # Parallel exec      |  OK         | Not OK
        #
        #
        # Strategy 2
        #                    | Lam booted  | Lam not booted
        #-----------------------------------------------------
        #
        # Sequential exec    |  Hangs      | Not OK
        # Parallel exec      |  Hangs      | OK
        #



# Initialise MPI
#
# Attempt to initialise mpiext.so
# If this fails, define a rudimentary interface suitable for
# sequential execution.

if error:
    print "WARNING: MPI library could not be initialised"

else:
    from mpiextras import \
         isend_array, \
         ireceive_array, \
         allreduce_array

    # Work around bug in OpenMPI (December 2009):
    # https://bugs.launchpad.net/ubuntu/+source/petsc4py/+bug/232036

    from ctypes import *
    CDLL('libmpi.so', RTLD_GLOBAL)
    # End work around








