/*
gcc -fPIC -c urs_ext.c -I/usr/include/python2.5 -o urs_ext.o -Wall -O
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
PyArrayObject *_read_mux2(int, char **, float *, double *, int);

PyObject *read_mux2(PyObject *self, PyObject *args){
/*Read in mux 2 file
   
    Python call:
    read_mux2(numSrc,filenames,weights,file_params,verbose)

    NOTE:
    A Python int is equivalent to a C long
    A Python double corresponds to a C double
*/
  PyObject *filenames;
  PyArrayObject *pyweights,*file_params;
  PyArrayObject *sts_data;
  PyObject *fname;

  char **muxFileNameArray;
  float *weights;
  long numSrc;
  long verbose;

  int nsta0;
  int nt;
  double dt;
  int i;
  
  // Convert Python arguments to C
  if (!PyArg_ParseTuple(args, "iOOOi",
			&numSrc,&filenames,&pyweights,&file_params,&verbose)) {			
			
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

  // Create array for weights which are passed to read_mux2
  weights = (float *) malloc((int)numSrc*sizeof(float));
  for (i=0;i<(int)numSrc;i++){
    weights[i]=(float)(*(double *)(pyweights->data+i*pyweights->strides[0]));
  }

  // Read in mux2 data from file
  sts_data=_read_mux2((int)numSrc,muxFileNameArray,weights,(double*)file_params->data,(int)verbose);


  // Allocate space for return vector
  nsta0=(int)*(double *) (file_params -> data+0*file_params->strides[0]);
  dt=*(double *) (file_params -> data+1*file_params->strides[0]);
  nt=(int)*(double *) (file_params -> data+2*file_params->strides[0]);

  free(weights);
  free(muxFileNameArray);
  return  PyArray_Return(sts_data);

}

PyArrayObject *_read_mux2(int numSrc, char **muxFileNameArray, float *weights, double *params, int verbose)
{
   FILE *fp;
   int nsta, nsta0, i, t, isrc, ista;
   struct tgsrwg *mytgs, *mytgs0;
   char *muxFileName;                                                                  
   int istart, istop;
   int *fros, *lros;
   char susMuxFileName;
   float *muxData;
   long numData;

   int start_tstep;
   int finish_tstep;
   int num_ts;
   int dimensions[2];
   int len_sts_data;
   float *temp_sts_data;
   PyArrayObject *sts_data;
   
   long int offset;

   /* Allocate space for the names and the weights and pointers to the data*/    
   
   /* Check that the input files have mux2 extension*/
   susMuxFileName=0;
   for(isrc=0;isrc<numSrc;isrc++)
     { 
       muxFileName=muxFileNameArray[isrc];
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
                     
   if (verbose){
     printf("Reading mux header information\n");
   }
   
   /* Open the first muxfile */
   if((fp=fopen(muxFileNameArray[0],"r"))==NULL){
      fprintf(stderr, "cannot open file %s\n", muxFileNameArray[0]);
      exit(-1);  
   }
 
   /* Read in the header */
   /* First read the number of stations*/   
   fread(&nsta0,sizeof(int),1,fp);
   
   /*now allocate space for, and read in, the structures for each station*/
   mytgs0 = (struct tgsrwg *) malloc(nsta0*sizeof(struct tgsrwg));
   fread(&mytgs0[0], nsta0*sizeof(struct tgsrwg), 1, fp);

   /*make an array to hold the start and stop steps for each station for each
     source*/   
   //FIXME: Should only store start and stop times for one source
   fros = (int *) malloc(nsta0*numSrc*sizeof(int));
   lros = (int *) malloc(nsta0*numSrc*sizeof(int));
   
   /* read the start and stop steps for source 0 into the array */   
   fread(fros,nsta0*sizeof(int),1,fp);
   fread(lros,nsta0*sizeof(int),1,fp);

   /* compute the size of the data block for source 0 */   
   numData = getNumData(fros, lros, nsta0);
   num_ts = mytgs0[0].nt;
   
   /* Burbidge: Added a sanity check here */
   if (numData < 0) {
     fprintf(stderr,"Size of data block appears to be negative!\n");
     exit(-1);
   }
   fclose(fp); 

   // Allocate header array for each remaining source file.
   // FIXME: only need to store for one source and which is compared to all subsequent
   // source headers
   if(numSrc > 1){
      /* allocate space for tgsrwg for the other sources */
      mytgs = (struct tgsrwg *)malloc( nsta0*sizeof(struct tgsrwg) );
   } else {
     /* FIXME (JJ): What should happen in case the are no source files?*/
     /* If we exit here, tests will break */
     // fprintf(stderr, "No source file specified\n");
     // exit(-1);       
   }
   
   /* loop over remaining sources, check compatibility, and read them into 
    *muxData */
   for(isrc=1; isrc<numSrc; isrc++){
     muxFileName = muxFileNameArray[isrc];
     
     /* open the mux file */
     if((fp=fopen(muxFileName,"r"))==NULL)
       {
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
       if(mytgs[ista].nt != num_ts){
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
	  exit(-1);
      }
      fclose(fp);             
   }
   params[0]=(double)nsta0;
   params[1]=(double)mytgs0[0].dt;
   params[2]=(double)mytgs0[0].nt;
  
    
   // Find min and max start times of all gauges
   start_tstep=num_ts+1;
   finish_tstep=-1;
   for (ista=0;ista<nsta0;ista++){
     if (fros[ista]< start_tstep){
       start_tstep=fros[ista];
     }
     if (lros[0+nsta0*ista] > finish_tstep){
       finish_tstep=lros[ista]; 
     }
   }
   
   if ((start_tstep>num_ts) | (finish_tstep < 0)){
     printf("ERROR: Gauge data has incorrect start and finsh times\n");
     return NULL;
   }
   
   if (start_tstep>=finish_tstep){
     printf("ERROR: Gauge data has non-postive_length\n");
     return NULL;
   }
   
   /* Make array(s) to hold the demuxed data */
   len_sts_data=num_ts+POFFSET;
   dimensions[0]=nsta0;
   dimensions[1]=len_sts_data;
   sts_data = (PyArrayObject *) PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
   if(sts_data == NULL){
    printf("ERROR: Memory for pydata array could not be allocated.\n");
    return NULL;
   }
   
   /* Initialise sts data to zero */
   for (ista=0;ista<nsta0;ista++){
     for (t=0;t<num_ts;t++){
       *(double *) (sts_data -> data+ista*sts_data->strides[0]+t*sts_data->strides[1])=0.0;
     }
     *(double *) (sts_data -> data+ista*sts_data->strides[0]+(num_ts)*sts_data->strides[1])=(float)mytgs0[ista].geolat;
     *(double *) (sts_data -> data+ista*sts_data->strides[0]+(num_ts+1)*sts_data->strides[1])=(float)mytgs0[ista].geolon;
     *(double *) (sts_data -> data+ista*sts_data->strides[0]+(num_ts+2)*sts_data->strides[1])=(float)mytgs0[ista].z;
     *(double *) (sts_data -> data+ista*sts_data->strides[0]+(num_ts+3)*sts_data->strides[1])=(float)fros[ista];
     *(double *) (sts_data -> data+ista*sts_data->strides[0]+(num_ts+4)*sts_data->strides[1])=(float)lros[ista];
   } 

   temp_sts_data = (float *)calloc(len_sts_data, sizeof(float) );
      
   /* Loop over all sources */
   //FIXME: remove istart and istop they are not used.
   istart = -1;
   istop = -1;
   for (isrc=0;isrc<numSrc;isrc++){
     /* Read in data block from mux2 file */
     muxFileName = muxFileNameArray[isrc];
     if((fp=fopen(muxFileName,"r"))==NULL)
       {
	 //fprintf(stderr, "%s: cannot open file %s\n", av[0], muxFileName);
	 fprintf(stderr, "cannot open file %s\n", muxFileName);
         exit(-1);  
       }
       
     if (verbose){
       printf("Reading mux file %s\n",muxFileName);
     }

     offset=sizeof(int)+nsta0*(sizeof(struct tgsrwg)+2*sizeof(int));
     fseek(fp,offset,0);
     
     numData = getNumData(fros+isrc*nsta0, lros+isrc*nsta0, nsta0);
     muxData = (float*) malloc(numData*sizeof(float));
     fread(muxData, numData*sizeof(float),1,fp); 
     fclose(fp); 	 

     /* loop over stations */
     for(ista=0; ista<nsta0; ista++){               
       
       /* fill the data0 array from the mux file, and weight it */
       fillDataArray(ista, nsta0, num_ts, mytgs0[ista].ig, fros+isrc*nsta0, lros+isrc*nsta0, temp_sts_data, &istart, &istop, muxData);
       
       /* weight appropriately and add */
       for(i=0; i<num_ts; i++){
	 if((isdata(*(double *) (sts_data -> data+ista*sts_data->strides[0]+i*sts_data->strides[1]))) && isdata(temp_sts_data[i])){
	   *(double *) (sts_data -> data+ista*sts_data->strides[0]+i*sts_data->strides[1]) += temp_sts_data[i] * weights[isrc];
	 }else{
	   *(double *) (sts_data -> data+ista*sts_data->strides[0]+i*sts_data->strides[1]) = 0.0;
	 }
       }
     }
   }
   free(muxData);
   free(temp_sts_data);
   return sts_data;
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
