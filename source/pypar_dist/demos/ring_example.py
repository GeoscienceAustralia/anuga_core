import pypar                                       # Import module and initialise MPI 

proc = pypar.size()                                # Number of processes as specified by mpirun
myid = pypar.rank()                                # Id of of this process (myid in [0, proc-1]) 
node = pypar.get_processor_name()                  # Host name on which current process is running

print 'I am proc %d of %d on node %s' %(myid, proc, node)

if myid == 0:                                      # Actions for process 0:
  msg = 'P0'  
  pypar.send(msg, destination=1)                   # Send message to proces 1 (right hand neighbour)
  msg = pypar.receive(source=proc-1)               # Receive message from last process
      
  print 'Processor 0 received message "%s" from processor %d' %(msg, proc-1)

else:                                              # Actions for all other processes:

  source = myid-1                                  # Source is the process to the left
  destination = (myid+1)%proc                      # Destination is process to the right
                                                   # wrapped so that last processor will 
						   # send back to proces 0  
  
  msg = pypar.receive(source)                      # Receive message from source 
  msg = msg + 'P' + str(myid)                      # Update message     
  pypar.send(msg, destination)                     # Send message to destination   

pypar.finalize()                                   # Stop MPI 
