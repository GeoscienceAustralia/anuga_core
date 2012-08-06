#########################################################
#
#  Example: Distributing work in parallel using pypar (MPI).
#
#  Author:  Ole Nielsen, ANU, March 2002.
#
#########################################################
#
"""Skeleton for using pypar to distribute work across processors.
   It is assumed that the function work knows how to perform different
   tasks (or the same task on different data) based on values of
   myid and numprocs.
"""



def work(myid, numprocs, data):
  """This simple example function that slices up the data
     based on values of numproc and myid.
  """

  import numpy

  # Identify local slice and process it
  #
  interval = len(data)
  myinterval = interval/numprocs

  mylower = myid*myinterval
  if myid == numprocs-1:
    myupper = interval+1
  else:
    myupper = mylower + myinterval

  mydata = data[mylower:myupper]

  # Computation (average)
  #
  myavg = float(numpy.sum(mydata))/len(mydata)
  print "P%d: %s Local avg=%.4f" %(myid, str(mydata), myavg)

  return myavg*len(mydata)




###################################################
# Main program - communication takes place here
#
import pypar, numpy


# Get data. Here it is just generated but it could be read
# from file or given as an input parameter.
#
lower = 100
upper = 121
data = numpy.array(range(lower,upper))

#
# Get parallel parameters
#
numprocs = pypar.size()    # Number of processors
myid = pypar.rank()        # Id of this processor
node = pypar.get_processor_name()

print "I am proc %d of %d on node %s" %(myid, numprocs, node)



#
# Do work in parallel
#
x = work(myid, numprocs, data)   #Do work on all processors
print "Proc %d finished working" %myid

#
# Communication
#
if numprocs > 1:
  #
  # Processor 0 gathers all results and merge them
  #
  if myid == 0:
    for id in range(1,numprocs):
      print "P%d receving from P%d" %(0, id)
      x = x + pypar.receive(id)  #Add up (would be more complex in general)

  # All other processors send their results back to processor 0
  #
  else:
    print "P%d sending to P%d" %(myid, 0)
    pypar.send(x, 0)

print "Proc %d after communication" %myid
#
# Compute overall average and report
#

if myid == 0:
  avg = x/len(data)
  print "Global average is %.4f" %avg

pypar.finalize()


