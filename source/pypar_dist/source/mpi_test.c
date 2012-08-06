/*
  Simple MPI communication test.
  Ole Moller Nielsen - 2011

  To compile
  mpicc mpi_test.c -lmpi -Wall

  To run on one node
  mpirun --hostfile /etc/mpihosts -host node5 -npernode 2 a.out

  To run on multiple nodes
  mpirun --hostfile /etc/mpihosts -host node5,node10 -npernode 1 a.out
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


#define M  5000     /* Data size */


int main(int argc, char **argv) {
   int repeats = 3, msgid = 0;
   int myid, procs;
   int j, k;

   double A[M];

   int  namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];

   MPI_Status stat;


   /* Initialize */
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Get_processor_name(processor_name, &namelen);

   if (myid == 0) {
     printf("Number of processes = %d\n", procs);
     printf("Test repeated %d times for reliability\n", repeats);
   }

   if (procs < 2) {
     printf("Program needs at least two processors - aborting\n");
     MPI_Abort(MPI_COMM_WORLD,999);
   }

   /* Create the data */
   for (j=0; j<M; j++) {
      A[j]=rand();
   }

   /* Synchronize */
   MPI_Barrier(MPI_COMM_WORLD);
   printf("I am process %d on node %s\n", myid, processor_name);

   /* Pass msg circularly a number of times*/
   for (k=0; k<repeats; k++) {
     if (myid == 0) {
       printf("Run %d of %d\n", k+1, repeats);
     }

     /* Communicate*/
     if (myid == 0) {
	 printf("P%i: Sending to P%i\n", myid, 1);
         MPI_Send(&A[0], M, MPI_DOUBLE, 1, msgid, MPI_COMM_WORLD);
	 printf("P%i: Waiting to receive from P%i\n", myid, procs-1);
         MPI_Recv(&A[0], M, MPI_DOUBLE, procs-1, msgid, MPI_COMM_WORLD, &stat);
	 printf("P%i: Received from to P%i\n", myid, procs-1);
       } else {
	 printf("P%i: Waiting to receive from to P%i\n", myid, myid-1);
         MPI_Recv(&A[0], M, MPI_DOUBLE, myid-1, msgid, MPI_COMM_WORLD, &stat);
	 printf("P%i: Sending to to P%i\n", myid, (myid+1)%procs);
         MPI_Send(&A[0], M, MPI_DOUBLE, (myid+1)%procs, msgid, MPI_COMM_WORLD);
     }
   }


   printf("P%i: Done\n", myid);

   MPI_Finalize();
   return 0;
}

