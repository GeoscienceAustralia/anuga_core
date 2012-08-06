# -*- coding: iso-8859-1 -*-

#########################################################################
#########################################################################

__title__ = "bsend examples"
__author__ = "Constantinos Makassikis"
__date__ = "18/09/2007"

import numpy
import string
import time
import sys, os
try:
	import pypar
except:
	raise 'Could not find module pypar'

#
# Hereafter, can be found a set of functions implementing a broadcast.
#
# Besides being used for testing purposes, they illustrate different
# ways of using 'bsend()' function. w/o and w/ buffers.
# 
# Understanding how 'MPI_Bsend()' is used in MPI might be useful. 
#

#
# In short, when one wants to send data using 'MPI_Bsend()' he has to 
# allocate a buffer big enough to host all the data that needs to be
# send.
#
# buffer = malloc();
# MPI_Attach(buffer);
# MPI_Bsend(buffer, data1, ...);
# MPI_Bsend(buffer, data2, ...);
# MPI_Detach();
# free(buffer);
#

def simple_broadcast_v1(x, root, myid, numprocs):
	""" Broadcast implementation using bsend() w/o explicitly defining 
	receive buffers.

	Input parameters:
	 - x: data to broadcast.
	 - root: rank of the process that initiates the broadcast.
	 - myid: rank of the process calling the function.

	Return value: the broadcasted data.
	"""

	list = range(0, root) + range(root + 1, numprocs)

	# The process of rank 'root' sends 'x' to all the other processes.
	if root == myid:
		for i in list:

			# Determine automatically the size of bsend's buffer.
			pypar.push_for_alloc(x)						

			# Allocate and attach bsend's buffer.	
			pypar.alloc_and_attach()							

			# Send data to process of rank 'i'.
			pypar.bsend(x, i)				

		 	# Detach and deallocate bsend's buffer.
			pypar.detach_and_dealloc()							

		x_ = x

	# All processes with rank distinct from 'root' start receiving.
	else:
		x_ = pypar.receive(root)
		
	return x_

def simple_broadcast_v2(x, root, myid, numprocs):
	""" Optimized version of 'simple_broadcast_v1()'.

	Since every bsend needed to achieve the broadcast requires a buffer
	of the same size, it is sufficient to allocate that buffer once 
	and to use it for all bsend calls.

	This illustrates the use of functions alloc(), attach(), dealloc()
	detach() instead of mpi_alloc_and_attach() and 
	mpi_detach_and_dealloc(). As it can be seen the former alternatives
	give a finer control.

	Note: a call to 
	
	...
	alloc_and_attach()
	...

	is equivalent to:

	...
	alloc()
	attach()
	...

	Input parameters:
	 - x: data to broadcast.
	 - root: rank of the process that initiates the broadcast.
	 - myid: rank of the process calling the function.

	Return value: the broadcasted data.
	"""
	list = range(0, root) + range(root + 1, numprocs)

	# The process of rank 'root' sends 'x' to all the other processes.
	if root == myid:
		
		# Determine automatically the size of bsend's buffer.
		pypar.push_for_alloc(x)						

		# Allocate bsend's buffer.
		pypar.alloc()							

		for i in list:

			# Attach bsend's buffer.
			pypar.attach()							

			# Send data to process of rank 'i'.
			pypar.bsend(x, i)				

		 	# Detach bsend's buffer.
			pypar.detach()							

		x_ = x

		# Deallocate bsend's buffer.
		pypar.dealloc()							

	# All processes with rank distinct from 'root' start receiving.
	else:
		x_ = pypar.receive(root)
		
	return x_

def simple_broadcast_v3(x, root, myid, numprocs):
	""" Alternative implementation to 'simple_broadcast_v1()'.

	Make all bsend calls using a single buffer and a
	single call to 'mpi_alloc_and_attach()' and 'mpi_detach_and_dealloc()'.

	This implementation illustrates how to use 'push_for_alloc()' in order
	to allocate a single buffer, big enough, to handle all sends within
	one call to 'mpi_alloc_and_attach()' and 'mpi_detach_and_dealloc()'.

	Input parameters:
	 - x: data to broadcast.
	 - root: rank of the process that initiates the broadcast.
	 - myid: rank of the process calling the function.

	Return value: the broadcasted data.
	"""

	list = range(0, root) + range(root + 1, numprocs)

	# The process of rank 'root' sends 'x' to all other processes.
	if root == myid:
		# Determine the size of bsend's buffer.
		for i in list:
			pypar.push_for_alloc(x)						

		# Allocate and attach bsend's buffer.	
		pypar.alloc_and_attach()							

		# Make all sends.
		for i in list:
			pypar.bsend(x, i)				

		# Deallocate and detach bsend's buffer.
		pypar.detach_and_dealloc()							

		x_ = x

	# All processes with rank distinct from 'root' start receiving.
	else:
		x_ = pypar.receive(root)
		
	return x_

def simple_broadcast_v4(x, root, myid, numprocs, buffer):
	""" Broadcast implementation using bsend() w/ explicit definition of
	receive buffers.

	Input parameters:
	 - x: data to broadcast.
	 - root: rank of the process that initiates the broadcast.
	 - myid: rank of the process calling the function.
	 - buffer: well-dimensioned user-defined buffer for receiving 'x'.

	Return value: the broadcasted data.
	"""
	list = range(0, root) + range(root + 1, numprocs)

	# The process of rank 'root' sends 'x' to all the other processes.
	if root == myid:
		for i in list:

			# Determine automatically the size of bsend's buffer.
			pypar.push_for_alloc(x, use_buffer=True)

			# Allocate and attach bsend's buffer.	
			pypar.alloc_and_attach()							

			# Send data to process of rank 'i'.
			pypar.bsend(x, i, use_buffer=True)				

		 	# Detach and deallocate bsend's buffer.
			pypar.detach_and_dealloc()							

		buffer = x

	# All processes with rank distinct from 'root' start receiving.
	else:
		buffer = pypar.receive(root, buffer)
		
	return buffer

def simple_broadcast_v5(x, root, myid, numprocs, buffer):
	""" Same as simple_broadcast_v4 except that it uses Pypar's
	bypass mode.

	The use of bypass mode implies that the programmer has to define his 
	own buffers on the receiving side (same as when 'use_buffer' is True) 
	and that he is limited to send numpy arrays and nothing else !

	Hence, this function works only with numpy arrays.

	Input parameters:
	 - x: data to broadcast.
	 - root: rank of the process that initiates the broadcast.
	 - myid: rank of the process calling the function.
	 - buffer: well-dimensioned user-defined buffer for receiving 'x'.

	Return value: the broadcasted data.
	"""
	list = range(0, root) + range(root + 1, numprocs)

	# The process of rank 'root' sends 'x' to all the other processes.
	if root == myid:
		for i in list:

			# Determine automatically the size of bsend's buffer.
			pypar.push_for_alloc(x, bypass=True)						

			# Allocate and attach bsend's buffer.	
			pypar.alloc_and_attach()							

			# Send data to process of rank 'i'.
			pypar.bsend(x, i, bypass=True)				

		 	# Detach and deallocate bsend's buffer.
			pypar.detach_and_dealloc()							

		buffer = x

	# All processes with rank distinct from 'root' start receiving.
	else:
		buffer = pypar.receive(root, buffer, bypass=True)
		
	return buffer

#
# Main
#

if __name__ == "__main__":

	datetime = time.ctime(time.time())

	numprocs = pypar.size()
	myid = pypar.rank()
	hostname = pypar.get_processor_name()

	#print "%s - I am proc %d of %d on %s" %(datetime, myid, numprocs, hostname)

	# Testing broadcast of an array - v1.
	if myid == 0:
		snd_tab = numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
	else:
		snd_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])

	rcv_tab = simple_broadcast_v1(snd_tab, 0, myid, numprocs) 

	print "v1 - rank " + str(myid) + " received: " + str(rcv_tab)

	pypar.barrier()

	# Testing broadcast of an array - v2.
	if myid == 0:
		snd_tab = numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
	else:
		snd_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])

	rcv_tab = simple_broadcast_v2(snd_tab, 0, myid, numprocs) 

	print "v2 - rank " + str(myid) + " received: " + str(rcv_tab)

	pypar.barrier()

	# Testing broadcast of an array - v3.
	if myid == 0:
		snd_tab = numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
	else:
		snd_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])

	rcv_tab = simple_broadcast_v3(snd_tab, 0, myid, numprocs) 

	print "v3 - rank " + str(myid) + " received: " + str(rcv_tab)

	pypar.barrier()

	# Testing broadcast of an array - v4.
	if myid == 0:
		snd_tab = numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
	else:
		snd_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])

	rcv_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])
	rcv_tab = simple_broadcast_v4(snd_tab, 0, myid, numprocs, rcv_tab) 

	print "v4 - rank " + str(myid) + " received: " + str(rcv_tab)

	pypar.barrier()

	# Testing broadcast of an array - v5.
	if myid == 0:
		snd_tab = numpy.array([[1, 2, 3, 4], [5, 6, 7, 8]])
	else:
		snd_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])

	rcv_tab = numpy.array([[0, 0, 0, 0], [0, 0, 0, 0]])
	rcv_tab = simple_broadcast_v5(snd_tab, 0, myid, numprocs, rcv_tab) 

	print "v5 - rank " + str(myid) + " received: " + str(rcv_tab)

	pypar.barrier()

	# Testing broadcast of string - v1.
	if myid == 0:
		snd_str = "<< Mes sama in corpore sano. >>"
	else:
		snd_str = "Failed !!!"
		
	rcv_str = simple_broadcast_v1(snd_str, 0, myid, numprocs) 

	print "v1 - rank " + str(myid) + " received: " + str(rcv_str)

	pypar.barrier()

	# Testing broadcast of string - v2.
	if myid == 0:
		snd_str = "<< Mes sama in corpore sano. >>"
	else:
		snd_str = "Failed !!!"
		
	rcv_str = simple_broadcast_v2(snd_str, 0, myid, numprocs) 

	print "v2 - rank " + str(myid) + " received: " + str(rcv_str)

	pypar.barrier()

	# Testing broadcast of string - v3.
	if myid == 0:
		snd_str = "<< Mes sama in corpore sano. >>"
	else:
		snd_str = "Failed !!!"
		
	rcv_str = simple_broadcast_v3(snd_str, 0, myid, numprocs) 

	print "v3 - rank " + str(myid) + " received: " + str(rcv_str)

	pypar.barrier()

	# Testing broadcast of a string - v4.
	if myid == 0:
		snd_str = "<< Mes sama in corpore sano. >>"
	else:
		snd_str = "Failed !!!"
		
	rcv_str = "0000000000000000000000000000000"
	rcv_str = simple_broadcast_v4(snd_str, 0, myid, numprocs, rcv_str) 

	print "v4 - rank " + str(myid) + " received: " + str(rcv_str)

	pypar.barrier()

	"""
	# For those who want anyway to try out sending something different
	# from numpy arrays with bypass mode enabled.
	# Testing broadcast of a string - v5.
	if myid == 0:
		snd_str = "<< Mes sama in corpore sano. >>"
	else:
		snd_str = "Failed !!!"
		
	rcv_str = "0000000000000000000000000000000"
	rcv_str = simple_broadcast_v5(snd_str, 0, myid, numprocs, rcv_str) 

	print "v5 - rank " + str(myid) + " received: " + str(rcv_str)
	"""

	pypar.barrier()
	pypar.finalize()

