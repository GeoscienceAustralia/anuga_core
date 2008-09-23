#####################################
# Compilation of C extension mpi.c  #
# for use with pypar                #
#	                            #
# Ole Nielsen, Jan 2002             #
##################################### 

CC          =   mpicc #MPI C compiler


### Dependencies ###
#
default: mpi.so

mpi.so: mpi.c Makefile install.py
	python install.py

ctiming: ctiming.c Makefile
	$(CC) ctiming.c -o ctiming

clean:
	rm -f *.o so_locations *~ runpytiming.o* runpytiming.e*

veryclean: clean
	rm -f *.so *.pyc ctiming
