#!/usr/bin/env python
# Inital test of MPI module 'pypar' for Python
# 
# Run as 
#   python test_init.py
# or 
#   mpirun -np 2 test_init.py
#  (perhaps try number of processors more than 2)
#
# To verify bandwidth of your architecture please run pytiming (and ctiming) 
#
# OMN April 2002


import pypar

myid =    pypar.rank()
numproc = pypar.size()
node =    pypar.Get_processor_name()

print "I am processor %d of %d on node %s" %(myid, numproc, node)
if myid == 0:
    print 'Maximal tag is %d' %pypar.max_tag
pypar.Finalize()

