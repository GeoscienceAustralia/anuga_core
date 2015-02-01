/* this appears to be necessary to run on Linux */

#undef _STDLIB_H_
#define _STDLIB_H_


/*

    This code interfaces pmesh directly with "triangle", a general        
    purpose triangulation code. In doing so, Python numeric data structures
     are passed to C for  use by "triangle" and the results are then          
    converted back. This was accomplished using the Python/C API.            
                                                                           
    I. Arrays of points,edges,holes, and regions must normally be 
     created for proper triangulation. However, if there are no holes,
     the  hole list need not be specified (however, a 
     measure of control over the shape of the convex hull is lost; 
     triangle decides). 
     
     points array: [(x1,y1),(x2,y2),...] ( doubles)                  
     edge array: [(Ea,Eb),(Ef,Eg),...] ( integers)                   
     hole array: [(x1,y1),...]( doubles, one inside each hole region)
     region array: [ (x1,y1,index,area),...] ( doubles)                
     pointattribute array:[ (a1,a2..),...] (doubles) 
     mode: String - the commands sent to triangle
          
     This is all that is needed to accomplish an intial triangulation.       
     A dictionary of the Triangle output is returned.
            
     I thought this function was leaking. That the returned data 
     structure didn't get garbage collected.  I decided this using  
     gc.get_objects() and seeing the dic still hanging around.
     I'm not so sure this means its leaking.  Check by looking at 
     the overall memory use.  Anyhow, I delete the lists that are 
     returned in the dic structre after they are used now, in alpha
     shape and mesh_engine.

          
     Precondition
     End list in the pointattributelist has to have the same length
     
                                                 */


#ifdef SINGLE  
#define REAL float 
#else  
#define REAL double 
#endif 

#include "triangle.h" 

#ifndef _STDLIB_H_ 
#ifndef _WIN32
extern "C" void *malloc(); 
extern "C" void free(); 
#endif
#endif  

#include "Python.h" 

#include "util_ext.h"
#include "numpy_shim.h"


//#define PY_ARRAY_UNIQUE_SYMBOL API_YEAH 

#define PY_ARRAY_UNIQUE_SYMBOL API_YEAH
//#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"
#include <sys/types.h>

 static PyObject *triang_genMesh(PyObject *self, PyObject *args){
     
  /* Mesh generation routine
     pre-conditions:
     mode must have a 'z' switch, to avoid segmentation faults
     
     shallow_water_ext.c gravity function
  */
    
  struct triangulateio in,out;

  struct triangulateio in_test;
  
  PyArrayObject *pointlist;
  PyArrayObject *seglist;
  PyArrayObject *holelist;
  PyArrayObject *regionlist;
  PyObject *mode;
  PyArrayObject *pointattributelist;
  PyArrayObject *segmarkerlist;
  PyArrayObject *test;
  
  
  PyArrayObject *gentrianglelist;
  PyArrayObject *genpointlist;
  PyArrayObject *genseglist;
  PyArrayObject *genpointmarkerlist;
  PyArrayObject *genpointattributelist;
  PyArrayObject *gentriangleattributelist;
  PyArrayObject *gensegmentlist;
  PyArrayObject *gensegmentmarkerlist;
  PyArrayObject *genneighborlist;
  PyObject *R;

  int dimensions[2];
    
  REAL Attr;
  int i, j, iatt, n, write_here,N;
  int a,b,c;
  int marker;
  int tuplesize;
    
  int index = 0;
  double x,y;
  PyObject *elems;
  PyObject *pair;
     
  char *mod;
    
  if(!PyArg_ParseTuple(args,(char *)"OOOOOOO",&pointlist,&seglist,&holelist,&regionlist,&pointattributelist,&segmarkerlist,&mode)){
    return NULL;
  }

  // check that numpy array objects arrays are C contiguous memory
  CHECK_C_CONTIG(pointlist);
  CHECK_C_CONTIG(seglist);
  CHECK_C_CONTIG(holelist);
  CHECK_C_CONTIG(regionlist);
  CHECK_C_CONTIG(pointattributelist);
  CHECK_C_CONTIG(segmarkerlist);
  
  /* Initialize  points */
  in.numberofpoints =  pointlist-> dimensions[0]; 
  in.pointlist = (double *) pointlist -> data;
  in.pointmarkerlist = (int *)NULL; 
  
  /* Initialize and fill up region list */
  in.numberofregions = regionlist-> dimensions[0]; 
  
  in.regionlist = (double *) regionlist -> data;
  
  /*Initialize segments and segment markers */
  in.numberofsegments = seglist -> dimensions[0]; 
  in.segmentlist = (int *) seglist -> data;
  in.segmentmarkerlist = (int *) segmarkerlist -> data;
  //in.segmentmarkerlist = (int *) NULL;
  
  /*Initialize triangles */
  in.numberoftriangles = 0;
  in.trianglelist = (int *)NULL;
  in.numberoftriangleattributes = 0;  
  in.triangleattributelist  = (REAL *)NULL; 
  in.trianglearealist  = (REAL *)NULL; 
  in.neighborlist = (int *)NULL;
  in.numberofcorners = 0;
     
  /*Initialize holes */
  in.numberofholes = holelist  -> dimensions[0];
  if (in.numberofholes != 0) {
    in.holelist =  (double *) holelist -> data; 
  } else {
    in.holelist = (REAL *)NULL;
  }     
  
  /* Point attribute list */
   /*printf ("in.pointattributelist -> dimensions[0] %i\n", pointattributelist -> dimensions[0]);
  printf ("in.pointattributelist -> dimensions[1] %i\n", pointattributelist -> dimensions[1]); */
  if (0 == pointattributelist -> dimensions[0]) {
    in.numberofpointattributes = 0;
    in.pointattributelist =  (double *) NULL;
  } else {
    if (0 == pointattributelist -> dimensions[1]) {
    in.numberofpointattributes = 0;
    in.pointattributelist =  (double *) NULL;
    } else {
      in.numberofpointattributes = pointattributelist -> dimensions[1];
      in.pointattributelist =  (double *) pointattributelist -> data;
    }
  }
   
  /* DEBUGGING 
  printf(" ***  hello world\n" );  
 
  printf ("numberofpoints %i\n", in.numberofpoints);
  printf ("numberofpointattributes %i\n", in.numberofpointattributes);
  for(i = 0; i < in.numberofpoints; i++) {
    printf("(%f,%f)\n",in.pointlist[i* 2+ 0],in.pointlist[i* 2+ 1]);
    for(j = 0; j < in.numberofpointattributes; j++) {
      printf("point att (%f)\n",in.pointattributelist[i* in.numberofpointattributes + j]);
    }
    
  }
      
  printf ("numberofregions %i\n", in.numberofregions);
  for(i = 0; i < in.numberofregions; i++) {
    printf("(%f,%f)  ",in.regionlist[i* 4+ 0],in.regionlist[i* 4+ 1]);
    printf("index %f Area %f\n",in.regionlist[i* 4+ 2],in.regionlist[i* 4+ 3]);
  }
   
  printf ("numberofsegments %i\n", in.numberofsegments);
  for(i = 0; i < in.numberofsegments; i++) {
      printf("(%i,%i)\n",in.segmentlist[i* 2+ 0],in.segmentlist[i* 2+ 1]);
      printf("Segment marker (%i)\n",in.segmentmarkerlist[i + 0]);
      }
 
 
  printf ("numberofholess %i\n", in.numberofholes);
  for(i = 0; i < in.numberofholes; i++) {
    printf("(%f,%f)\n",in.holelist[i* 2+ 0],in.holelist[i* 2+ 1]);
  }
     
      printf ("numberoftriangles %i\n", in.numberoftriangles);
      for(i = 0; i < in.numberoftriangles; i++) {
      printf("(%i,%i,%i)",in.trianglelist[i* 3+ 0],in.trianglelist[i* 3+ 1], in.trianglelist[i* 3+ 2]);
      } 
  printf(" ***  see ya world\n" );       
      */
  /* set up the switch arguments to the triangulation routine */
  mod = PyString_AS_STRING(mode);
     
         
  out.pointlist = (REAL *)NULL;
  out.pointmarkerlist = (int *)NULL;
  //out.numberofpointattributes = in.numberofpointattributes;
  out.pointattributelist = (REAL *)NULL;
  
  out.trianglelist = (int *)NULL;
  out.triangleattributelist = (REAL *)NULL;                   
  out.trianglearealist = (REAL *)NULL;                        
  out.neighborlist = (int *)NULL;
     
  out.segmentlist = (int *)NULL;
  out.segmentmarkerlist = (int *)NULL;
     
  out.edgelist = (int *)NULL;
  out.edgemarkerlist = (int *)NULL;
     
  out.holelist = (REAL *)NULL;
  out.regionlist = (REAL *)NULL;

  triangulate(mod, &in, &out, (struct triangulateio *)NULL );
  
  
    
  /*
    ------- Pass point numbers,coordinates and neighbors back to Python ------
    we send back a dictionary:                                               
    { index : [ coordinates, [connections], Attribute ] } 
  */  
 
  /*
  
     to return numeric arrays, check how it is done in 
     abs quantity_ext.c compute_gradients
     
  PyArray_FromDims allolws you to create a numeric array with unitialized data.
   The first argument is the size of the second argument (
   the dimensions array).
    The dimension array argument is just a 1D C array where each element of
     the array is the size of that dimension. 
     (int dimensions[2] = { 4, 3 }; defines a 4 by 3 array.) 
     The third argument is just the desired type.
  */
  
  
  /* Add triangle list */
  dimensions[0] = out.numberoftriangles;
  dimensions[1] = 3;   
  gentrianglelist = (PyArrayObject *) anuga_FromDimsAndData(2,
					    dimensions,
					    PyArray_INT, 
  					    (char*) out.trianglelist); 
    
  /* Add pointlist */
  dimensions[0] = out.numberofpoints;
  dimensions[1] = 2;   
  genpointlist = (PyArrayObject *) anuga_FromDimsAndData(2,
					 dimensions,
					 PyArray_DOUBLE, 
					 (char*) out.pointlist);
					   
 
  /* Add point marker list */
  dimensions[0] = out.numberofpoints;
  genpointmarkerlist = (PyArrayObject *) anuga_FromDimsAndData(1, 
				  	 dimensions, 
				         PyArray_INT,
				        (char*) out.pointmarkerlist);
  
  /* Add point attribute list */
  dimensions[0] = out.numberofpoints;
  dimensions[1] = out.numberofpointattributes;   
  genpointattributelist = (PyArrayObject *) anuga_FromDimsAndData(2, 
					  dimensions, 
					  PyArray_DOUBLE,
					  (char*) out.pointattributelist);
  
  
 
  /* Add triangle attribute list */
  dimensions[0] = out.numberoftriangles;
  dimensions[1] = out.numberoftriangleattributes;   
  gentriangleattributelist = (PyArrayObject *) anuga_FromDimsAndData(2, 
					   dimensions, 
					   PyArray_DOUBLE,
					  (char*)out.triangleattributelist);
  
  /* Add segment list */
  dimensions[0] = out.numberofsegments;
  dimensions[1] = 2;   
  gensegmentlist = (PyArrayObject *) anuga_FromDimsAndData(2, 
						    dimensions, 
						    PyArray_INT,
						    (char*)out.segmentlist);
  
  
  /* Add segment marker list */
  dimensions[0] = out.numberofsegments;
  gensegmentmarkerlist = (PyArrayObject *) anuga_FromDimsAndData(1, 
					    dimensions, 
					    PyArray_INT,
					   (char*)out.segmentmarkerlist);
  
  /* Add triangle neighbor list */
  if (out.neighborlist != NULL) {
    dimensions[0] = out.numberoftriangles;
    dimensions[1] = 3;   
    genneighborlist = (PyArrayObject *) anuga_FromDimsAndData(2, 
					 	 dimensions, 
					      	 PyArray_INT,
					       	 (char*)out.neighborlist);
  }else{ 
    dimensions[0] = 0;
    dimensions[1] = 0;   
    genneighborlist = (PyArrayObject *) anuga_FromDims(2, 
						       dimensions, 
						       PyArray_INT);
  }
  
  
  R = Py_BuildValue((char *)"OOOOOOOO"
		    ,PyArray_Return(gentrianglelist)
		    ,PyArray_Return(genpointlist)
		    ,PyArray_Return(genpointmarkerlist)
		    ,PyArray_Return(genpointattributelist)
		    ,PyArray_Return(gentriangleattributelist)
		    ,PyArray_Return(gensegmentlist)
		    ,PyArray_Return(gensegmentmarkerlist)
		    ,PyArray_Return(genneighborlist)
		    );
		    
  Py_DECREF(gentrianglelist);
  Py_DECREF(genpointlist);
  Py_DECREF(genpointmarkerlist);
  Py_DECREF(genpointattributelist);
  Py_DECREF(gentriangleattributelist);
  Py_DECREF(gensegmentlist);
  Py_DECREF(gensegmentmarkerlist);
  Py_DECREF(genneighborlist);
  
  
  /* These memory blocks are passed into numeric arrays 
   so don't free them  */
  /*
  if(!out.trianglelist){    
    free(out.trianglelist); out.trianglelist=NULL;
    
    
  if(!out.pointlist){
    free(out.pointlist);  out.pointlist=NULL;
  }  
  
  if(!out.pointmarkerlist){    
    free(out.pointmarkerlist); out.pointmarkerlist=NULL;
  }
  
  if(!out.pointattributelist){    
    free(out.pointattributelist); out.pointattributelist=NULL;
  }   
  
  if(!out.triangleattributelist){    
    free(out.triangleattributelist); out.triangleattributelist=NULL;
  }
  if(!out.segmentlist){
    free(out.segmentlist); out.segmentlist =NULL;
  }
  if(!out.segmentmarkerlist){
    free(out.segmentmarkerlist); out.segmentmarkerlist  =NULL;
  }
  if(!out.neighborlist){    
    free(out.neighborlist); out.neighborlist=NULL;
  }
   */

    
  
  /* Free in/out structure memory */
  
  if(!out.trianglearealist){    
    free(out.trianglearealist); out.trianglearealist=NULL;
  }
  if(!out.edgelist){
    free(out.edgelist);  out.edgelist=NULL;
  }
  if(!out.edgemarkerlist){
    free(out.edgemarkerlist);  out.edgemarkerlist=NULL;
  }
  if(!out.holelist){
    free(out.holelist); out.holelist=NULL;
  }
  if(!out.regionlist){
    free(out.regionlist); out.regionlist=NULL;
  }
  return R;
}

/* end of my function*/
 
static PyMethodDef triang_methods[] = {
  {(char *)"genMesh", triang_genMesh, 1, (char *)"Passes hull points to triangle.All arguments are lists.(<hull>,<seglist>,<holelist>,<regionlist>)-><hull>"}, 
  {NULL,NULL}
};

void initmesh_engine_c_layer(){
  Py_InitModule((char *)"mesh_engine_c_layer",triang_methods);
  
  import_array(); // Necessary for handling of NumPY structures
}    
