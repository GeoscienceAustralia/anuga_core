#########################################################
#
#  Example: Running Python in parallel using pypar (MPI).
#
#  Author:  Ole Nielsen,  SMS, ANU, Jan. 2002.
#
#########################################################
#
# The purpose of this code is to demonstrate how Python can be
# used to communicate among processes using pypar.
#
# This demo passes messages on in a ring from processor n-1 to n starting
# and ending with processor 0.
# Each processor adds some text to the message before passing it on
# and writes a log statement to the screen.
#
# To execute:
#
#   mpirun -np 4 demo.py


import pypar    # The Python-MPI interface
import numpy
typ = 'd'  #'i', 'l', 'f', 'd'

numproc = pypar.size()
myid =    pypar.rank()
node =    pypar.Get_processor_name()

print "I am proc %d of %d on node %s" %(myid, numproc, node)

if numproc < 2:
  print "Demo must run on at least 2 processors to continue"
  pypar.abort()

if myid == 0:
  msg = numpy.array([0], typ)

  print 'Processor 0 sending message "%s" to processor %d' %(str(msg), 1)
  pypar.send(msg, 1)

  msg, status = pypar.receive(numproc-1, return_status=True)
  print 'Processor 0 received message "%s" from processor %d' %(str(msg), numproc-1)
  print 'Size of msg was %d bytes' %(status.bytes())

else:
  source = myid-1
  destination = (myid+1)%numproc

  msg, status = pypar.receive(source, return_status=True)
  print 'Processor %d received message "%s" from processor %d'\
        %(myid, str(msg), source)
  print 'Size of msg was %d bytes' %(status.bytes())

  msg = list(msg)
  msg.append(myid)
  print msg
  msg = numpy.array(msg, typ)  #Update message
  print 'Processor %d sending msg "%s" to %d' %(myid, msg, destination)
  pypar.send(msg, destination)

pypar.finalize()



