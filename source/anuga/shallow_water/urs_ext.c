/*
gcc -fPIC -c urs_ext.c -I/usr/include/python2.4 -o urs_ext.o -Wall -O
gcc -shared urs_ext.o  -o urs_ext.so
*/
#include "Python.h"
#include "Numeric/arrayobject.h"
#include "math.h"
#include <stdio.h>
#include <float.h>
#include <time.h>
#include "structure.h"

#define MAX_FILE_NAME_LENGTH 128
#define NODATA 99.0
#define EPSILON  0.00001

#define DEBUG 0

#define POFFSET 5 //Number of site_params

void fillDataArray(int, int, int, int, int *, int *, float *, int *, int *, float *);
long getNumData(int*, int*, int);
char isdata(float);
void wrttxt(char *, float, int, float *, float, float, float, float, float, int, int, int);
float** _read_mux2(int, char **, float *, double *, int);

PyObject *read_mux2(PyObject *self, PyObject *args){
/*Read in mux 2 file
   
    Python call:
    read_mux2(numSrc,filenames,weights,file_params,write)

    NOTE:
    A Python int is equivalent to a C long
    A Python double corresponds to a C double
*/
  PyObject *filenames;
  PyArrayObject *pyweights,*file_params;
  PyArrayObject *pydata;
  PyObject *fname;

  char **muxFileNameArray;
  float *weights;
  long numSrc;
  long write;

  float **cdata;
  int dimensions[2];
  int nsta0;
  int nt;
  double dt;
  int i,j;
  int start_tstep;
  int finish_tstep;
  int it,time;
  int num_ts;
  
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "iOOOi",
			&numSrc,&filenames,&pyweights,&file_params,&write)) {			
			
    PyErr_SetString(PyExc_RuntimeError, 
		    "Input arguments to read_mux2 failed");
    return NULL;
  }

  if(!PyList_Check(filenames)) {
    PyErr_SetString(PyExc_TypeError, "get_first_elem expects a list");
    return NULL;
  }

  if(PyList_Size(filenames) == 0){
    PyErr_SetString(PyExc_ValueError, "empty lists not allowed");
    return NULL;
  }

  if (pyweights->nd != 1 || pyweights->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError,
		    "pyweights must be one-dimensional and of type double");
    return NULL; 
  }

  if(PyList_Size(filenames) != pyweights->dimensions[0]){
      PyErr_SetString(PyExc_ValueError, "Must specify one weight for each filename");
      return NULL;
  }

  muxFileNameArray = (char **) malloc((int)numSrc*sizeof(char *));
  if (muxFileNameArray == NULL) {
     printf("ERROR: Memory for muxFileNameArray could not be allocated.\n");
     exit(-1);
  }
  for (i=0;i<PyList_Size(filenames);i++){
    muxFileNameArray[i] = (char *) malloc((MAX_FILE_NAME_LENGTH+1)*sizeof(char));
    if (muxFileNameArray[i] == NULL) {
      printf("ERROR: Memory for muxFileNameArray could not be allocated.\n");
      exit(-1);
    }
    fname=PyList_GetItem(filenames, i);
    if (!PyString_Check(fname)) {
      PyErr_SetString(PyExc_ValueError, "filename not a string");
    }
    muxFileNameArray[i]=PyString_AsString(fname);
  }

  if (file_params->nd != 1 || file_params->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_ValueError,
		    "file_params must be one-dimensional and of type double");
    return NULL; 
  }

  weights = (float *) malloc((int)numSrc*sizeof(float));
  for (i=0;i<(int)numSrc;i++){
    weights[i]=(float)(*(double *)(pyweights->data+i*pyweights->strides[0]));
  }

  cdata=_read_mux2((int)numSrc,muxFileNameArray,weights,(double*)file_params->data,(int)write);


  // Allocate space for return vector
  nsta0=(int)*(double *) (file_params -> data+0*file_params->strides[0]);
  dt=*(double *) (file_params -> data+1*file_params->strides[0]);
  nt=(int)*(double *) (file_params -> data+2*file_params->strides[0]);

  start_tstep=nt+1;
  finish_tstep=-1;
  for (i=0;i<nsta0;i++){
    if ((int)cdata[i][nt+3] < start_tstep){
      start_tstep=(int)cdata[i][nt+3];
    }
    if ((int)cdata[i][nt+4] > finish_tstep){
      finish_tstep=(int)cdata[i][nt+4]; 
    }
  }

  if ((start_tstep>nt) | (finish_tstep < 0)){
    printf("ERROR: Gauge data has incorrect start and finsh times\n");
    return NULL;
  }

  if (start_tstep>=finish_tstep){
    printf("ERROR: Gauge data has non-postive_length\n");
    return NULL;
  }

  num_ts=finish_tstep-start_tstep+1;
  dimensions[0]=nsta0;
  dimensions[1]=num_ts+POFFSET;
  pydata = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  if(pydata == NULL){
    printf("ERROR: Memory for pydata array could not be allocated.\n");
    return NULL;
  }

  for (i=0;i<nsta0;i++){
    time=0;
    for (it=0;it<finish_tstep;it++){
      if((it+1>=start_tstep)&&(it+1<=finish_tstep)){
	if (it+1>(int)cdata[i][nt+4]){
	  //This gauge has stopped recording but others are still recording
	*(double *) (pydata -> data+i*pydata->strides[0]+time*pydata->strides[1])=0.0;
	}else{
	  *(double *) (pydata -> data+i*pydata->strides[0]+time*pydata->strides[1])=cdata[i][it];
	}
	time++;
      }
    }
    for (j=0; j<POFFSET; j++){
      *(double *) (pydata -> data+i*pydata->strides[0]+(num_ts+j)*pydata->strides[1])=cdata[i][nt+j];
     }
  }

  free(weights);
  free(muxFileNameArray);
  return  PyArray_Return(pydata);

}

float** _read_mux2(int numSrc, char **muxFileNameArray, float *weights, double *params, int write)
//void _read_mux2(int numSrc, char **muxFileNameArray, float *weights, float **mydata)
{
   FILE *fp;
   int nsta, nsta0, i, isrc, ista;//numSrc;
   struct tgsrwg *mytgs, *mytgs0;
   //char *muxFileNameArray;
   char *muxFileName;                                                                  
   char outFileName[MAX_FILE_NAME_LENGTH+1];
   float *wt;
   float max, amax;
   float *data, *data0;
   int istart, istop;
   int *fros, *lros;
   char susMuxFileName;
   float **muxData;
   long numData;
   time_t   start_time, stop_time;

   float **mydata;
   
   /* note starting time */
   time(&start_time);

   /* allocate space for the names and the weights and pointers to the data*/    
   wt = (float *) malloc(numSrc*sizeof(float));
   muxData = (float**) malloc(numSrc*sizeof(float*));
   
   /*read the mux file names and the associated weight from stdin*/     
   susMuxFileName=0;
   for(isrc=0;isrc<numSrc;isrc++)
     { 
       muxFileName=muxFileNameArray[isrc];
       wt=weights;
       if(!susMuxFileName && strcmp(muxFileName+strlen(muxFileName)-4,"mux2")!=0){
	 susMuxFileName=1;
       }
     }
   
   if(susMuxFileName)
   {
      printf("\n**************************************************************************\n");
      printf("   WARNING: This program operates only on multiplexed files in mux2 format\n"); 
      printf("   At least one input file name does not end with mux2\n");
      printf("   Check your results carefully!\n");
      printf("**************************************************************************\n\n");
   }   
                     
   //printf("Demuxing mux files with weights:\n\n");
   for(isrc=0;isrc<numSrc;isrc++){
     muxFileName = muxFileNameArray[isrc];
     //printf("%s %f\n", muxFileName, *(wt+isrc));
   } 
   //printf("\n");
   
   /* open the first muxfile */
   if((fp=fopen(muxFileNameArray[0],"r"))==NULL)
   {
      fprintf(stderr, "cannot open file %s\n", muxFileNameArray[0]);
      exit(-1);  
   }
 
   /* read in the header */
   /*first read the number of stations*/   
   fread(&nsta0,sizeof(int),1,fp);
   /*now allocate space for, and read in, the structures for each station*/
   mytgs0 = (struct tgsrwg *) malloc(nsta0*sizeof(struct tgsrwg));
   fread(&mytgs0[0], nsta0*sizeof(struct tgsrwg), 1, fp);

   /*make an array to hold the start and stop steps for each station for each source*/   
   fros = (int *) malloc(nsta0*numSrc*sizeof(int));
   lros = (int *) malloc(nsta0*numSrc*sizeof(int));
   
   /* read the start and stop steps for source 0 into the array */   
   fread(fros,nsta0*sizeof(int),1,fp);
   fread(lros,nsta0*sizeof(int),1,fp);

   /* compute the size of the data block for source 0 */   
   numData = getNumData(fros, lros, nsta0);

   /* Burbidge: Added a sanity check here */
   if (numData < 0) {
     fprintf(stderr,"Size of data block appears to be negative!\n");
     //fprintf(stderr,"numData=%d fros=%d lros=%d nsta0=%d\n",numData,fros,lros,nsta0);
     exit(-1);
   }
   
   /* allocate space for these data, read them and close the file */   
   *muxData = (float*) malloc(numData*sizeof(float));
   fread(*muxData, numData*sizeof(float),1,fp);
   fclose(fp); 

   if(numSrc > 1){
      /* allocate space for tgsrwg for the other sources */
      mytgs = (struct tgsrwg *)malloc( nsta0*sizeof(struct tgsrwg) );
   } else {
     /* FIXME (JJ): What should happen in case the are no source files?*/
     /* If we exit here, tests will break */
     // fprintf(stderr, "No source file specified\n");
     // exit(-1);       
   }
   
   /* loop over sources, check compatibility, and read them into *muxData */
   for(isrc=1; isrc<numSrc; isrc++){
     muxFileName = muxFileNameArray[isrc];
     
     /* open the mux file */
     if((fp=fopen(muxFileName,"r"))==NULL)
       {
	 //fprintf(stderr, "%s: cannot open file %s\n", av[0], muxFileName);
	 fprintf(stderr, "cannot open file %s\n", muxFileName);
         exit(-1);  
       }
     
     /* check that the mux files are compatible */      
     fread(&nsta,sizeof(int),1,fp);
     if(nsta != nsta0){
       fprintf(stderr,"%s has different number of stations to %s\n", muxFileName, muxFileNameArray[0]);
       fclose(fp);
       exit(-1);   
     }
     fread(&mytgs[0], nsta*sizeof(struct tgsrwg), 1, fp);
     for(ista=0; ista < nsta; ista++){
       if(mytgs[ista].dt != mytgs0[ista].dt){
	 fprintf(stderr,"%s has different sampling rate to %s\n", muxFileName, muxFileNameArray[0]);
	 fclose(fp);
	 exit(-1);                 
       }   
       if(mytgs[ista].nt != mytgs0[ista].nt){
	 fprintf(stderr,"%s has different series length to %s\n", muxFileName, muxFileNameArray[0]);
	 fclose(fp);
	 exit(-1);                 
       }
     }

      /* read the start and stop times for this source */
      fread(fros+isrc*nsta0,nsta0*sizeof(int),1,fp);
      fread(lros+isrc*nsta0,nsta0*sizeof(int),1,fp);
      
      /* compute the size of the data block for this source */
      numData = getNumData(fros+isrc*nsta0, lros+isrc*nsta0, nsta0);

      /* Burbidge: Added a sanity check here */
      if (numData < 0){
	  fprintf(stderr,"Size of data block appears to be negative!\n");
	  //fprintf(stderr,"numData=%d fros=%d lros=%d nsta0=%d\n");
	  exit(-1);
      }

      /* allocate space, read the data and close the file */
      *(muxData+isrc) = (float*) malloc(numData*sizeof(float));
      fread(*(muxData+isrc), numData*sizeof(float),1,fp);
      fclose(fp);             
   }
   params[0]=(double)nsta0;
   params[1]=(double)mytgs0[0].dt;
   params[2]=(double)mytgs0[0].nt;
   
   /* make array(s) to hold the demuxed data */
   mydata = (float **)malloc (nsta0 * sizeof(float *));
   if (mydata == NULL){
     printf("ERROR: Memory for mydata could not be allocated.\n");
     exit(-1);
   }
   for (ista=0; ista<nsta0; ista++){
     mydata[ista] = (float *)malloc( (mytgs0[0].nt + POFFSET)* sizeof(float) );
     if (mydata[ista] == NULL){
       printf("ERROR: Memory for mydata could not be allocated.\n");
       exit(-1);
     }
     mydata[ista][mytgs0[0].nt]=(float)mytgs0[ista].geolat;
     mydata[ista][mytgs0[0].nt+1]=(float)mytgs0[ista].geolon;
     mydata[ista][mytgs0[0].nt+2]=(float)mytgs0[ista].z;
     mydata[ista][mytgs0[0].nt+3]=(float)fros[ista];
     mydata[ista][mytgs0[0].nt+4]=(float)lros[ista];
   }

   //data0 = (float *)malloc( mytgs0[0].nt * sizeof(float) );
   if(numSrc > 1)
      data = (float *)malloc( mytgs0[0].nt * sizeof(float) );
         
   /* loop over stations */
   for(ista=0; ista<nsta0; ista++){              
     data0= mydata[ista];
     data=data0;
     istart = -1;
     istop = -1;
     /* fill the data0 array from the first mux file, and weight it */
     isrc=0;
     //muxFileName = muxFileNameArray + isrc*(MAX_FILE_NAME_LENGTH+1);
     muxFileName = muxFileNameArray[isrc];
     //printf("Demuxing station %d\n",ista);
     //printf("   ... source %s\n",muxFileName);
     
     fillDataArray(ista, nsta0, mytgs0[ista].nt, mytgs0[ista].ig, fros+isrc*nsta0, lros+isrc*nsta0, data0, &istart, &istop, *(muxData+isrc));
     
     /* apply the weights */
     for(i=0; i<mytgs0[ista].nt; i++)
       if(isdata(*(data0+i)))
	 *(data0+i) *= *wt;
     
     /* loop over the rest of the sources */
     for(isrc=1;isrc<numSrc;isrc++)
       {
         /* fill the data array */
         //muxFileName = muxFileNameArray + isrc*(MAX_FILE_NAME_LENGTH+1);
	 muxFileName = muxFileNameArray[isrc];
         //printf("   ... source %s\n",muxFileName);
         fillDataArray(ista, nsta0, mytgs0[ista].nt, mytgs0[ista].ig, fros+isrc*nsta0, lros+isrc*nsta0, data, &istart, &istop, *(muxData+isrc)); 
         
         /* weight appropriately and add */
         for(i=0; i<mytgs0[ista].nt; i++)
	   {
	     if(isdata(*(data0+i))&&isdata(*(data+i)))
               *(data0+i) += *(data+i)* *(wt+isrc);                                                                   
	     else
               *(data0+i) = NODATA;
	   }      
	 
       }  /* end of loop over sources */ 
     
     /* now compute the maxima for this station */      
     max = 0.0;
     amax = 0.0;
     for(i=0; i<mytgs0[ista].nt; i++){
       if(isdata(*(data0+i)))
         {
	   max = ((*(data0+i) > max) ? *(data0+i):max);
	   amax = ((fabs(*(data0+i)) > amax) ? fabs(*(data0+i)):amax);  
         }   
     }
     /* write out sac file for the current station  */
     /*thomas - instead of passing beg(=0), should pass dt*/
     /*thomas - uncomment the following if you want sac files*/
     /*sprintf(outFileName,"S%03d.sac", ista );
       printf("   ... writing file=%s\n", outFileName);
       wrtsac2(outFileName, mytgs0[ista].dt, mytgs0[ista].nt, &data0[0], mytgs0[ista].dt, 
       mytgs0[ista].geolat, mytgs0[ista].geolon,
       mytgs0[ista].center_lat, mytgs0[ista].center_lon, 
       mytgs0[ista].offset, mytgs0[ista].z );*/
     
     /* write out text file for the current station  */
     /* DB added check for non-zero max */
     if (max >0) 
       {
	 /* Burbidge: added tide gauge grid id output. no .id field in tgsrwg */
	 /* thomas: instead of passing beg(=0), pass dt */   
	 /*      wrttxt(outFileName, mytgs0[ista].dt, mytgs0[ista].nt, &data0[0], beg, mytgs0[ista].geolat, 
		 mytgs0[ista].geolon, max, mytgs0[ista].z, mytgs0[ista].ig ); */
	 if (write) {
	   sprintf(outFileName,"S%5.5d.txt", ista );
	   //printf("   ... writing file=%s\n", outFileName);
	   wrttxt(outFileName, mytgs0[ista].dt, mytgs0[ista].nt, &data0[0], mytgs0[ista].dt, mytgs0[ista].geolat, 
		  mytgs0[ista].geolon, max, mytgs0[ista].z, mytgs0[ista].ig, istart, istop );
	 }
       }
     
   }  /* end of loop over stations */
   //free(data0);
   //free(mytgs0);
   
   //if(numSrc>1)
   //{
     //free(data);
     //free(mytgs);
   //} 
   //can't free arrays because I only fill array by making pointer not copy
   //for(isrc=0; isrc<numSrc;isrc++)
   //   free(*(muxData+isrc));
   //   
   //free(muxData);
   // free(muxFileNameArray);
   
   if(susMuxFileName)
   {
      printf("\n**************************************************************************\n");
      printf("   WARNING: This program operates only on multiplexed files in mux2 format\n"); 
      printf("   At least one input filename does not end with mux2\n");
      printf("   Check your results carefully!\n");
      printf("**************************************************************************\n");
   }   
   time(&stop_time);
   //fprintf(stdout, "\nElapsed time: %10.0f seconds\n", difftime(stop_time, start_time));
   return mydata;
}   


/* thomas */
void fillDataArray(int ista, int nsta, int nt, int ig, int *nst, int *nft, float *data, int *istart_p, int *istop_p, float *muxData)
{
   int it, last_it, jsta;
   long int offset=0;


   last_it=-1;
   /* make arrays of starting and finishing time steps for the tide gauges */
   /* and fill them from the file */
      
   /* update start and stop timesteps for this gauge */
   if(nst[ista]!=-1){
     if(*istart_p==-1){
         *istart_p=nst[ista];
     }else{
         *istart_p=((nst[ista]<*istart_p)?nst[ista]:*istart_p);
     }
   }
   if(nft[ista]!=-1){
     if(*istop_p==-1){
         *istop_p=nft[ista];
     }else{
         *istop_p=((nft[ista]<*istop_p)?nft[ista]:*istop_p);
     }
   }     
   if(ig==-1 || nst[ista] == -1) /* currently ig==-1 => nst[ista]==-1 */
   {
      /* gauge never started recording, or was outside of all grids, fill array with 0 */
      for(it=0; it<nt; it++)
         data[it] = 0.0;
   }   
   else
   {
      for(it=0; it<nt; it++)
      {
         last_it = it;
         /* skip t record of data block */
         offset++;
         /* skip records from earlier tide gauges */
         for(jsta=0; jsta<ista; jsta++)
            if(it+1>=nst[jsta]&&it+1<=nft[jsta])
               offset++;
                
         /* deal with the tide gauge at hand */
         if(it+1>=nst[ista]&&it+1<=nft[ista])
         /* gauge is recording at this time */
         {
            memcpy(data+it,muxData+offset,sizeof(float));
            offset++;
         }
         else if (it+1<nst[ista])
         {
            /* gauge has not yet started recording */
            data[it] = 0.0;
         }   
         else
         /* gauge has finished recording */                                            
         {
            data[it] = NODATA;
            break;
         }
   
         /* skip records from later tide gauges */
         for(jsta=ista+1; jsta<nsta; jsta++)
            if(it+1>=nst[jsta]&&it+1<=nft[jsta])
               offset++;
      }
   
      if(last_it < nt - 1)
         /* the loop was exited early because the gauge had finished recording */
         for(it=last_it+1; it < nt; it++)
            data[it] = NODATA;
   }
} 

/* Burbidge: No "offset" is sent. Replace with max. Added grid_id */
void wrttxt( fname, dt, nt, x, beg, lat, lon, max, depth, grid_id, istart, istop )
char *fname;
float dt, *x, beg, max, lat, lon, depth;
int grid_id;
int nt;
int istart, istop;
{
   int it;
   float t;
   FILE *fp;

   fp = fopen(fname,"w");
   fprintf(fp,"> lat %g lon %g depth %g max %g start_time %g stop_time %g grid_id %d\n", lat, lon, depth, max, 
            istart*dt, istop*dt, grid_id );
   for(it=0;it<nt;it++)
   {
      t=beg+it*dt;
      fprintf(fp,"%9.3f %g\n", t, x[it]);
   }
   fclose(fp);
}
  
char isdata(float x)
{
  //char value;
   if(x < NODATA + EPSILON && NODATA < x + EPSILON)
      return 0;
   else
      return 1;  
}

long getNumData(int *fros, int *lros, int nsta)
/* calculates the number of data in the data block of a mux file */
/* based on the first and last recorded output steps for each gauge */ 
{
   int ista, last_output_step;
   long numData = 0;

   last_output_step = 0;   
   for(ista=0; ista < nsta; ista++)
      if(*(fros + ista) != -1)
      {
         numData += *(lros + ista) - *(fros + ista) + 1;
         last_output_step = (last_output_step < *(lros+ista) ? *(lros+ista):last_output_step);
      }   
   numData += last_output_step*nsta; /* these are the t records */
   return numData;
}   



//-------------------------------
// Method table for python module
//-------------------------------
static struct PyMethodDef MethodTable[] = {
  {"read_mux2", read_mux2, METH_VARARGS, "Print out"},
  {NULL, NULL}
};

// Module initialisation
void initurs_ext(void){
  Py_InitModule("urs_ext", MethodTable);

  import_array(); // Necessary for handling of NumPY structures
}
