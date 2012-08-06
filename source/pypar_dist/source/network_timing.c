/*
  Estimate bandwidth and latency of a parallel computer using MPI.
  Ole Moller Nielsen - 1998
*/
	
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


#define MAXI  10         /* Number of blocks */
#define MAXM  500000     /* Largest block */
#define BLOCK MAXM/MAXI  /* Block size */


double linfit(double* x, double* y, int N, double* a, double* b) 
{
  /* Given vectors y and x fit a and b to the model y = ax + b */
  
  double Sx=0, Sy=0, SxoN=0, SSoN=0, t=0;
  double res, varest=0, norm=0; 
  int i;
  
  for (i=0; i<N; i++)
  {
    /*printf("x,y = %f, %f\n",x[i],y[i]);*/
    Sx  = Sx + x[i];
    Sy  = Sy + y[i];
  }

  SxoN = Sx/N;
  
  *a = 0.0; 
  for (i=0; i<N; i++)
  {
    t    = x[i] - SxoN;
    SSoN = SSoN + t*t;
    *a   = *a + t*y[i];
  }

  *a = (*a)/SSoN;          /* a = (N Sxy - SxSy)/(NSxx - Sx^2) */
  *b = (Sy - Sx*(*a))/N;
  
  /* Quality - variance estimate \sum_i r_i^2 /(m-n) */
  for (i=0; i<N; i++)
  {
    norm = norm + x[i]*x[i];
    res = y[i] - (*a)*x[i] - (*b);
    varest = varest + res*res;
  } 
  varest = varest/norm/(N-2); 
  return(varest);
}

main(int argc, char **argv) 
{
   int repeats = 10, msgid = 0;
   int myid, procs;
   int i,j,k,m;

   double t1, t2, cpuOH; 
   double Tbw, Tlat;
   double varest;
      
   int noelem[MAXI];
   double bytes[MAXI];   
   double mintime[MAXI];   
   double maxtime[MAXI];      
   double avgtime[MAXI];         
   double A[MAXM]; 

   int  namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   MPI_Status stat;
  
   
   /* Initialize */

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&procs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   MPI_Get_processor_name(processor_name,&namelen);
   
   if (myid==0)
   {
     printf("MAXM = %d, number of processors = %d\n",MAXM,procs);         
     printf("Measurements are repeated %d times for reliability\n",repeats);
   } 
   
   if (procs < 2) {
     printf("Program needs at least two processors - aborting\n");
     MPI_Abort(MPI_COMM_WORLD,999);
   }
   
   MPI_Barrier(MPI_COMM_WORLD); /* Synchronize */   
   printf("I am process %d on %s\n",myid,processor_name);   
      
   for (j=0; j<MAXM; j++) 
   {
      A[j]=rand();
   }   
   for (i=0; i<MAXI; i++) 
   {
      avgtime[i] =  0;         
      mintime[i] =  1000000;      
      maxtime[i] = -1000000;            
   }
   
   /* Determine timer overhead */
   if (myid == 0) {   
     cpuOH = 1.0;
     for (k=0; k<repeats; k++)   /* Repeat to get reliable timings */
     { 
       t1 = MPI_Wtime();
       t2 = MPI_Wtime();
       if (t2-t1 < cpuOH) cpuOH = t2-t1;
     }  
     printf("Timing overhead is %f seconds\n\n", cpuOH);              
   }

   
        
   /* Pass msg circularly */
     
   for (k=0; k<repeats; k++) {
     if (myid == 0) {  
       printf("Run %d of %d\n", k+1, repeats);
     }  

     for (i=0; i<MAXI; i++) {
       /*m=BLOCK*(i+1);*/
       m=BLOCK*i+1;       
      
       noelem[i] = m;
      
       MPI_Barrier(MPI_COMM_WORLD); /* Synchronize */
      
       if (myid == 0) {
         t1=MPI_Wtime();
         MPI_Send(&A[0],m,MPI_DOUBLE,1,msgid,MPI_COMM_WORLD);
         MPI_Recv(&A[0],m,MPI_DOUBLE,procs-1,msgid,MPI_COMM_WORLD,&stat);
         t2=MPI_Wtime() - t1 - cpuOH;
	 t2 = t2/procs;
	 avgtime[i] = avgtime[i] + t2;
         if (t2 < mintime[i]) mintime[i] = t2;
         if (t2 > maxtime[i]) maxtime[i] = t2;	 
       } else {
         MPI_Recv(&A[0],m,MPI_DOUBLE,myid-1,msgid,MPI_COMM_WORLD,&stat);
         MPI_Send(&A[0],m,MPI_DOUBLE,(myid+1)%procs,msgid,MPI_COMM_WORLD);
       }
     } 
   }

   if (myid == 0) {
     printf("Bytes transferred   time (micro seconds)\n");
     printf("                    min        avg        max \n");     
     printf("----------------------------------------------\n");      
        
     for (i=0; i<MAXI; i++) {
       avgtime[i] = avgtime[i]/repeats*1.0e6; /*Average micro seconds*/
       mintime[i] = mintime[i]*1.0e6;         /*Min micro seconds*/       
       maxtime[i] = maxtime[i]*1.0e6;         /*Min micro seconds*/              
             
       m = noelem[i];
       bytes[i] = (double) 8*noelem[i];       
	 
       /* printf("m=%d, time(min)=%lf, time(avg)=%lf, time(max)=%lf\n",
	       m,mintime[i],avgtime[i],maxtime[i]); */
       printf("%10d    %10d %10d %10d\n",
	 (int) bytes[i], (int) mintime[i], (int) avgtime[i], (int)maxtime[i]); 
     }		 
   
     varest=linfit(bytes, mintime, MAXI, &Tbw, &Tlat);
     printf("\nLinear regression on best timings (t = t_l + t_b * bytes):\n");

     printf("  t_b = %f\n  t_l = %f\n", Tbw, Tlat);
     printf("  Estimated relative variance = %.9f\n\n",varest);     

     printf("Estimated bandwith (1/t_b):  %.3f Mb/s\n", (1.0/Tbw));         
     printf("Estimated latency:           %d micro s\n", 
             (int) (mintime[0] - (float) bytes[0]* (float)Tbw));         
     
   }

   MPI_Finalize();
}

