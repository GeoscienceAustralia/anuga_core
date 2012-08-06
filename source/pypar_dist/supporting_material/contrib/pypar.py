# =============================================================================
# pypar.py - Parallel Python using MPI
# Copyright (C) 2001, 2002 Ole M. Nielsen, Gian Paolo Ciceri          
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
# Contact addresses: Ole.Nielsen@anu.edu.au, gp.ciceri@acm.org         
#
# version 1.2.1, 16 February 2002                                      
#   Status block, MPI_ANY_TAG, MPI_ANY_SOURCE exported                 
# Version 1.2, 15 February 2002                                       	   
#   Scatter added by Gian Paolo Ciceri                                   
# Version 1.1, 14 February 2002                                       	   
#   Bcast added by Gian Paolo Ciceri                                   
# Version 1.0.2, 10 February 2002                                          
#   Modified by Gian Paulo Ciceri to allow pypar run under Python 2.2  
# Version 1.0.1, 8 February 2002                                       
#   Modified to install on SUN enterprise systems under Mpich          
# Version 1.0, 7 February 2002                                         
#   First public release for Python 2.1 (OMN)                          
# =============================================================================

"""Module pypar.py - Parallel Python using MPI

Public functions:

size() -- Number of processors
rank() -- Id of current processor
Get_processor_name() -- Return host name of current node

send() -- Blocking send (all types)
receive() -- Blocking receive (all types)
raw_send() -- Blocking send (Numeric arrays and strings)
raw_receive() -- Blocking receive (Numeric arrays and strings)
bcast() -- Broadcast
Wtime() -- MPI wall time
Barrier() -- Synchronisation point. Makes processors wait until all processors
             have reached this point.
Abort() -- Terminate all processes. 
Finalize() -- Cleanup MPI. No parallelism can take place after this point. 


See doc strings of individual functions for detailed documentation.
"""

import os, sys

# -----------------------------------------------------------------------------
# Options directory with default values - to be set by user
#

options = { 
  'vanilla_bufsize': None   # Buffersize (chars) used for vanilla sends
}

# Constants
#
control_tag = 32            # Tag used to identify control information
control_data_max_size = 256 # Maximal size of string holding control data


#------------------------------------------------------------------------
# MPI Status block (can be queried by user after each receive
#------------------------------------------------------------------------

class Status:
  def __init__(self):
    pass
  
  def set_values(self, status_tuple):
    self.source = status_tuple[0] 
    self.tag = status_tuple[1] 
    self.error = status_tuple[2] 
    self.length = status_tuple[3]
    
status = Status() #Initialise status object      

#---------------------------------------------------------------------------
# Communication functions
#--------------------------------------------------------------------------


def raw_send(x, destination, tag=0, vanilla=0):
  """Wrapper for raw MPI send.
     Send x to destination with tag.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
     
     The variable x can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.

  """

  protocol = get_control_info(x, vanilla)[0]
  if protocol == 'array':
    send_array(x, destination, tag)  
  elif protocol == 'string':
    send_string(x, destination, tag)            
  else:  
    send_vanilla(x, destination, tag)

      
def raw_receive(x, source, tag=0, vanilla=0):
  """Wrapper for raw MPI receive.
     Receive something of same size as x from source with tag.
     
     Automatically determine appropriate protocol
     and call corresponding receive function.
     
     The variable x can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.
  """

  
  protocol = get_control_info(x, vanilla)[0]
  if protocol == 'array':
    err, stat = receive_array(x, source, tag) 
    if not err:
      status.set_values(stat)
    else:
      raise 'receive_array failed with error code %d' %err  
  elif protocol == 'string':
    err, stat = receive_string(x, source, tag)
    if not err:
      status.set_values(stat)
    else:
      raise 'receive_string failed with error code %d' %err  
  else:  
    x = receive_vanilla(x, source, tag)

  return x


def send(x, destination, tag=0, vanilla=0):
  """Wrapper for easy MPI send.
     Send x to destination with tag.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
     Also passes type and size information on as preceding message to
     simplify the receive call.
     
     The variable x can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.

  """
  import string

  control_info = get_control_info(x, vanilla)
  protocol = control_info[0]
  
  if protocol == 'array':
    send_control_info(control_info, destination)
    
    send_array(x, destination, tag)    
  elif protocol == 'string':
    send_control_info(control_info, destination)    
    
    send_string(x, destination, tag)          
  elif protocol == 'vanilla':
    from cPickle import dumps     
    s = dumps(x, 1)
    control_info[2] = str(len(s))

    send_control_info(control_info, destination)   
    
    send_string(s, destination, tag)
  else:
    raise "Unknown values for protocol: %s" %protocol    

      
def receive(source, tag=0):
  """Wrapper for easy MPI receive.
     Receive data from source with tag.
     
     Assumes preceding message containing protocol, type, size.
     Create appropriate buffer and receive data.
  """

  control_info = receive_control_info(source)
  
  protocol = control_info[0]
  typecode = control_info[1]
  size =     control_info[2]  
  
  if protocol == 'array':
    import Numeric
    x = Numeric.zeros(size,typecode)
    err, stat = receive_array(x, source, tag)    
    if not err:
      status.set_values(stat)
    else:
      raise 'receive_array failed with error code %d' %err  
  elif protocol == 'string':
    x = ' '*size
    err, stat = receive_string(x, source, tag)       
    if not err:
      status.set_values(stat)
    else:
      raise 'receive_string failed with error code %d' %err  
  elif protocol == 'vanilla':
    from cPickle import loads 
    s = ' '*size    
    err, stat = receive_string(s, source, tag)
    if not err:
      status.set_values(stat)
    else:
      raise 'receive_string failed with error code %d' %err  
    
    x = loads(s)
  else:
    raise "Unknown values for protocol: %s" %protocol
      
  return x  

def bcast(x, source, vanilla=0):
  """Wrapper for MPI bcast.
     Broadcast x from source.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
     
     The variable x can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.

  """

  protocol = get_control_info(x, vanilla)[0]
  if protocol == 'array':
    bcast_array(x, source)    
  elif protocol == 'string':
    bcast_string(x, source)          
  elif protocol == 'vanilla':
    from cPickle import loads, dumps 
    s = dumps(x, 1)
    s = s + ' '*int(0.1*len(s)) #safety 
    bcast_string(s, source)
    x = loads(s)
  else:
    raise "Unknown values for protocol: %s" %protocol  
    
  return x         

def raw_scatter(s, nums, d, source, vanilla=0):
  """Wrapper for MPI scatter.
     Scatter s in nums element to d (of the same nums size) from source.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
     
     The variable s can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.

  """

  protocol = get_control_info(s, vanilla)[0]
  if protocol == 'array':
    scatter_array(s, nums, d, source)    
  elif protocol == 'string':
    scatter_string(s, nums, d, source)          
  elif protocol == 'vanilla':
    raise "Protocol: %s unsupported for scatter" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol  
    
  return d         


def scatter(s, nums, source, vanilla=0):
  """Wrapper for easy MPI Scatter receive.
     Receive data from source with tag.
     
     Create appropriate buffer and receive data.
  """

  control_info = get_control_info(s)
  
  protocol = control_info[0]
  typecode = control_info[1]
  size =  nums   
  
  if protocol == 'array':
    import Numeric
    x = Numeric.zeros(size,typecode)
    scatter_array(s, size, x, source)    
  elif protocol == 'string':
    x = ' '*size
    scatter_string(s, size, x, source)          
  elif protocol == 'vanilla':
    raise "Protocol: %s unsupported for scatter" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol
      
  return x  

def raw_gather(s, nums, d, source, vanilla=0):
  """Wrapper for MPI gather.
     Gather s in nums element to d (of the same nums size) from source.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
     
     The variable s can be any (picklable) type, but
     Numeric variables and text strings will most efficient.
     Setting vanilla = 1 forces vanilla mode for any type.

  """

  protocol = get_control_info(s, vanilla)[0]
  if protocol == 'array':
    gather_array(s, nums, d, source)    
  elif protocol == 'string':
    gather_string(s, nums, d, source)          
  elif protocol == 'vanilla':
    raise "Protocol: %s unsupported for gather" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol  
    
  return d         


def gather(s, nums, source, vanilla=0):
  """Wrapper for easy MPI Gather receive.
     Receive data from source with tag.
     
     Create appropriate buffer and receive data.
  """

  control_info = get_control_info(s)
  
  protocol = control_info[0]
  typecode = control_info[1]
  s_size =  nums   
  
  if protocol == 'array':
    import Numeric
    x = Numeric.zeros(s_size * size(),typecode)
    gather_array(s, s_size, x, source)    
  elif protocol == 'string':
    x = ' '*s_size*size()
    gather_string(s, s_size, x, source)          
  elif protocol == 'vanilla':
    raise "Protocol: %s unsupported for gather" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol
      
  return x  


def raw_reduce(s, d, nums, op, source, vanilla=0):
  """Wrapper for MPI_Reduce.
     Reduce s in nums element to d (of the same nums size) at source
     applying operation op.
     
     Automatically determine appropriate protocol
     and call corresponding send function.
    
  """

  protocol = get_control_info(s, vanilla)[0]
  if protocol == 'array':
    reduce_array(s, d, nums, op, source)          
  elif (protocol == 'vanilla' or protocol == 'string'):
    raise "Protocol: %s unsupported for reduce" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol  
    
  return d         


def reduce(s, nums, op, source, vanilla=0):
  """Wrapper for easy MPI Gather receive.
     Receive data from source with tag.
     
     Create appropriate buffer and receive data.
  """

  control_info = get_control_info(s)
  
  protocol = control_info[0]
  typecode = control_info[1]
  s_size =  nums   
  
  if protocol == 'array':
    import Numeric
    x = Numeric.zeros(s_size * size(),typecode)
    reduce_array(s, x, s_size, op, source)    
  elif (protocol == 'vanilla' or protocol == 'string'):
    raise "Protocol: %s unsupported for reduce" %protocol
  else:
    raise "Unknown values for protocol: %s" %protocol
      
  return x  


     

#---------------------------------------------------------
# AUXILIARY FUNCTIONS
#---------------------------------------------------------
def get_control_info(x, vanilla=0):
  """Determine which protocol to use for communication:
     (Numeric) arrays, strings, or vanilla based x's type.

     There are three protocols:
     'array':   Numeric arrays of type 'i', 'l', 'f', or 'd' can be communicated
                with mpi.send_array and mpi.receive_array.
     'string':  Text strings can be communicated with mpi.send_string and
                mpi.receive_string.
     'vanilla': All other types can be communicated using the scripts send_vanilla and
                receive_vanilla provided that the objects can be serialised using
                pickle (or cPickle). The latter mode is less efficient than the
                first two but it can handle complex structures.

     Rules:
     If keyword argument vanilla == 1, vanilla is chosen regardless of 
     x's type.
     Otherwise if x is a string, the string protocol is chosen
     If x is an array, the 'array' protocol is chosen provided that x has one
     of the admissible typecodes.
  """

  protocol = 'vanilla'
  typecode = ' '
  size = '0'
  
  if not vanilla:
    #if type(x).__name__ == 'string':  # OK in Python 2.1 but not 2.2
    if type(x).__name__[0:3] == 'str': # Fixed by Gian Paolo Ciceri 10/2/2   
      protocol = 'string'
      typecode = 'c'
      size = len(x)
    elif type(x).__name__ == 'array':
      try:
        import Numeric

        typecode = x.typecode() 
        if typecode in ['i', 'l', 'f', 'd']:
          protocol = 'array'
          size = len(x)
        else:    
          print "WARNING (pypar.py): Numeric object type %s is not supported."\
                %(x.typecode())
          print "Only types 'i', 'l', 'f', 'd' are supported,",
          print "Reverting to vanilla mode."
          protocol = 'vanilla'
  
      except:
        print "WARNING (pypar.py): Numeric module could not be imported,",
        print "reverting to vanilla mode"
        protocol = 'vanilla'


  control_info = [protocol, typecode, str(size)]

  return control_info



#----------------------------------------------

def send_control_info(control_info, destination):
  """Send control info to destination
  """
  import string
  
  msg = string.join(control_info,',')
  send_string(msg, destination, control_tag)

  
def receive_control_info(source):
  """Receive control info from source
  """
  import string
  
  msg = ' '*control_data_max_size
  err, stat = receive_string(msg, source, control_tag)
  if err:
    raise Exception
  #NB: Do not set status block here  

  control_info = msg.split(',')
  assert len(control_info) == 3
  control_info[2] = int(control_info[2])

  return control_info




#
# Used only by raw communication
#
def send_vanilla(x, destination, tag=0):
  from cPickle import dumps
  from mpi import send_string as send
  
  s=dumps(x, 1)
  send(s, destination, tag)     
  return len(s)


def receive_vanilla(x, source, tag=0):
  from cPickle import loads, dumps
  from mpi import receive_string as receive


  #Create buffer of the right size
  #(assuming that x is similar to sent x).

  if options['vanilla_bufsize']:
    s = ' '*options['vanilla_bufsize']
  else:  
    s = dumps(x, 1)
    s = s + ' '*int(0.1*len(s)) #safety
                  
  receive(s, source, tag)

  return loads(s)
   

#----------------------- 
# ADM
#----------------------- 

def check_C_extension(sequential_allowed):
  """Verify existence of mpi.so.
  """   
  

  try:
    if sequential_allowed:  
      import mpi
      # This will crash on systems like the Alpha Server or the Sun
      # if program is run sequentially      
    else:  
      # A more general test suitable for the Alpha Server is
      fid = open('mpi.so', 'r')
      fid.close()
      # On the other hand, when using pypar as a package, we do 
      # not have access to mpi.so so we need link from parent directory to
      # pypar/mpi.so or pypar/mpi.c
  except:
    try:
      import compile
      compile.compile('mpi.c', 'mpicc', verbose = 0)
    except:
      raise "ERROR: Please compile C extension mpi.c - python install.py"
  


#----------------------------------------------------------------------------
# Initialise module
#----------------------------------------------------------------------------

# Determine if MPI program is allowed to run sequentially
# Attempting to check this automatically may case some systems to hang.
#
if sys.platform in ['osf1V5', 'sunos5']:  #Compaq AlphaServer or Sun
  sequential_allowed = 0
else:  
  sequential_allowed = 1

#print sequential_allowed  

# Attempt to compile mpi.so if necessary
#

check_C_extension(sequential_allowed) 


#Check if MPI module can be initialised
#
# Attempt to import and initialise mpi.so
# If this fails, define a rudimentary interface suitable for
# sequential execution.
#
#    

#


if sys.platform == 'osf1V5':
  error = os.system('python -c "import mpi" >/dev/null 2>/dev/null') 
  
  # On the sun we can't do this test - even in the parallel case...
  #
  #if sys.platform == 'osf1V5':
  #  error = os.system('python -c "import mpi" >/dev/null 2>/dev/null') 
  #else:
  #  error = 0  

  # The check is performed in a separate shell.
  # Reason: The Alpha server or the Sun cannot recover from a 
  # try: 
  #   import mpi
  # if mpi.so isn't there or if the program is run as sequential.
  # On LAM/Linux, this test may cause system to hang.
else:
  error = 0  


if error:
  print "WARNING: MPI library could not be initialised - running sequentially"
  parallel = 0

  # Define rudimentary functions to keep sequential programs happy
  
  def size(): return 1
  def rank(): return 0

  def Get_processor_name():
    import os
    try:
      hostname = os.environ['HOST']
    except:
      try:  
        hostname = os.environ['HOSTNAME']  
      except:
        hostname = 'Unknown'  

    return hostname
      

  def Abort():
    import sys
    sys.exit()

  def Finalize(): pass
  
  def Barrier(): pass  

  def Wtime():
    import time
    return time.time()
else:    
  #Import the C-extension and initialise MPI
  #
  parallel = 1
  from mpi import size, rank, Barrier, Wtime, Get_processor_name, Finalize, Abort, \
                  send_string, receive_string,\
                  send_array, receive_array, \
                  bcast_string, bcast_array, \
                  scatter_string, scatter_array, \
                  gather_string, gather_array, \
                  reduce_array, \
                  MPI_ANY_TAG as any_tag, MPI_ANY_SOURCE as any_source, \
                  mpi_MAX, mpi_MIN, mpi_SUM, mpi_PROD, mpi_LAND, mpi_BAND, \
                  mpi_LOR, mpi_BOR, mpi_LXOR, mpi_BXOR, mpi_MAXLOC, mpi_MINLOC, mpi_REPLACE

		  
  print "Proc %d: MPI initialised OK" %rank()    
		  





