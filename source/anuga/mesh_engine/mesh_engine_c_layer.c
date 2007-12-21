/* this appears to be necessary to run on Linux */

#undef _STDLIB_H_
#define _STDLIB_H_


/*

    This code interfaces pmesh directly with "triangle", a general        
    purpose triangulation code. In doing so, Python Numeric data structures
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

     to return numeric arrays, check how it is done in 
     quantity_ext.c compute_gradients
          
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

//#define PY_ARRAY_UNIQUE_SYMBOL API_YEAH 

#define PY_ARRAY_UNIQUE_SYMBOL API_YEAH
//#define NO_IMPORT_ARRAY
#include "Numeric/arrayobject.h"
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
  
  
  PyArrayObject *r_test;
  PyObject *R;

  int dimensions[1];
    
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
  PyObject *holder, *holderlist, *attributelist;
  int listsize,attsize ;
  PyObject  *ii;
  
  /* used for testing numeric arrays*/
  int n_test;     
  if(!PyArg_ParseTuple(args,(char *)"OOOOOOOO",&pointlist,&seglist,&holelist,&regionlist,&pointattributelist,&segmarkerlist,&mode,&test)){
    return NULL;
  }
  
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
  /* in.segmentlist = (int *) test -> data; */
  in.segmentmarkerlist = (int *) segmarkerlist -> data;
  
  /*Initialize triangles */
  in.numberoftriangles = 0;
  in.trianglelist = (int *)NULL;
  in.numberoftriangleattributes = 0;  
  in.triangleattributelist  = (REAL *)NULL; 
     
  /*Initialize holes */
  in.numberofholes = holelist  -> dimensions[0];
  if (in.numberofholes != 0) {
    in.holelist =  (double *) holelist -> data; 
  } else {
    in.holelist = (REAL *)NULL;
  }     
  
  /* Point attribute list */
  if (0 == pointattributelist -> dimensions[0]) {
    in.numberofpointattributes = 0;
    in.pointattributelist =  NULL;
  } else {
    if (0 == pointattributelist -> dimensions[1]) {
    in.numberofpointattributes = 0;
    in.pointattributelist =  NULL;
    } else {
      in.numberofpointattributes = pointattributelist -> dimensions[1];
      in.pointattributelist =  (double *) pointattributelist -> data;
    }
  }
   
  /* DEBUGGING 
  printf(" ***  hello world\n" );  
 
  printf ("numberofpoints %i\n", in.numberofpoints);
  for(i = 0; i < in.numberofpoints; i++) {
    printf("(%f,%f)\n",in.pointlist[i* 2+ 0],in.pointlist[i* 2+ 1]);
  }
      
  printf ("numberofregions %i\n", in.numberofregions);
  for(i = 0; i < in.numberofregions; i++) {
    printf("(%f,%f)\n",in.regionlist[i* 2+ 0],in.regionlist[i* 2+ 1]);
  }
   
  printf ("numberofsegments %i\n", in.numberofsegments);
  for(i = 0; i < in.numberofsegments; i++) {
      printf("(%i,%i)",in.segmentlist[i* 2+ 0],in.segmentlist[i* 2+ 1]);
      }
 
 
  printf ("numberofholess %i\n", in.numberofholes);
  for(i = 0; i < in.numberofholes; i++) {
    printf("(%f,%f)\n",in.holelist[i* 2+ 0],in.holelist[i* 2+ 1]);
  }
     
      printf ("numberoftriangles %i\n", in.numberoftriangles);
      for(i = 0; i < in.numberoftriangles; i++) {
      printf("(%i,%i,%i)",in.trianglelist[i* 3+ 0],in.trianglelist[i* 3+ 1], in.trianglelist[i* 3+ 2]);
      } 
  printf(" ***  see ya world\n" );       */
      
  /* set up the switch arguments to the triangulation routine */
  mod = PyString_AS_STRING(mode);
     
         
  out.pointlist = (REAL *)NULL;
  out.pointmarkerlist = (int *)NULL;
  out.numberofpointattributes = in.numberofpointattributes;
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
    
  
  /*printf("\n\nTriangulate input args: %s \n\n", mod); */
  triangulate(mod, &in, &out, (struct triangulateio *)NULL );
  
  
  /*
  PyArray_FromDims allolws you to create a Numeric array with unitialized data.
   The first argument is the size of the second argument (
   the dimensions array).
    The dimension array argument is just a 1D C array where each element of
     the array is the size of that dimension. 
     (int dimensions[2] = { 4, 3 }; defines a 4 by 3 array.) 
     The third argument is just the desired type.
  */
  
  //Py_Initialize();
  // Testing passing a numeric array out
  dimensions[0] = 4;
  // Allocate space for return vectors a and b (don't DECREF)
  r_test = (PyArrayObject *) PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
  
  /* printf(" ***  back from triangulate\n" );    */
  /*
    ------- Pass point numbers,coordinates and neighbors back to Python ------
    we send back a dictionary:                                               
    { index : [ coordinates, [connections], Attribute ] } 
  */
  holder = PyDict_New();    
  
  /* Add triangle list */
  listsize = out.numberoftriangles;
  /* printf(" out.numberoftriangles %i\n", out.numberoftriangles ); */
  holderlist = PyList_New(listsize);
  for(i=0; i<listsize;i++){
    PyObject *mlist = Py_BuildValue((char *)"(i,i,i)", 
				    out.trianglelist[i*3  ], out.trianglelist [i*3+1], out.trianglelist [i*3+2]);    
    PyList_SetItem(holderlist,i, mlist);
  }    
  ii=PyString_FromString("generatedtrianglelist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);
     
  /* Add pointlist */
  listsize = out.numberofpoints;
  holderlist = PyList_New(listsize);
     
  for(i=0; i<out.numberofpoints;i++){
    PyObject *mlist = Py_BuildValue((char *)"(d,d)", 
				    out.pointlist[i*2  ], out.pointlist[i*2+1]);  
    PyList_SetItem(holderlist,i, mlist);
  }  
  ii=PyString_FromString("generatedpointlist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);
  
  /* Add point marker list */
  listsize = out.numberofpoints;
  holderlist = PyList_New(listsize);
  
  for(i=0; i<listsize;i++){
    PyObject *mlist = Py_BuildValue((char *)"d", 
				    out.pointmarkerlist[i]);     
    PyList_SetItem(holderlist,i, mlist); 
  }  
  ii=PyString_FromString("generatedpointmarkerlist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);
    
  /* Add point attribute list */
  listsize = out.numberofpoints;
  holderlist = PyList_New(listsize);
  
  for(i=0; i<listsize;i++){
    attsize = out.numberofpointattributes;
    attributelist = PyList_New(attsize);
    for(iatt=0; iatt<attsize;iatt++){
      PyObject *mlist = Py_BuildValue((char *)"d", 
				      out.pointattributelist[i*attsize + iatt]);     
      PyList_SetItem(attributelist,iatt, mlist);
    }      
    PyList_SetItem(holderlist,i, attributelist);
  }  
  ii=PyString_FromString("generatedpointattributelist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);  
 
  
  /* Add triangle attribute list */
  listsize = out.numberoftriangles;
  holderlist = PyList_New(listsize);
     
  for(i=0; i<listsize; i++){
    attsize = out.numberoftriangleattributes;
    attributelist = PyList_New(attsize);       
    for(iatt=0; iatt<attsize;iatt++){
      PyObject *mlist = Py_BuildValue((char *)"d",out.triangleattributelist[i*attsize + iatt]);  
      PyList_SetItem(attributelist,iatt, mlist);
    }      
    PyList_SetItem(holderlist,i, attributelist);
  }  
  ii=PyString_FromString("generatedtriangleattributelist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);
  
  /* Add segment list */
  listsize = out.numberofsegments;
  holderlist = PyList_New(listsize);
  for(i=0; i<listsize;i++){
    PyObject *mlist = Py_BuildValue((char *)"(i,i)", 
				    out.segmentlist[i*2  ],
				    out.segmentlist [i*2+1]);
    PyList_SetItem(holderlist,i, mlist);
  }    
  ii=PyString_FromString("generatedsegmentlist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);  
  
  /* Add segment marker list */
  listsize = out.numberofsegments;
  holderlist = PyList_New(listsize);
  for(i=0; i<listsize;i++){
    PyObject *mlist = Py_BuildValue((char *)"i", 
				    out.segmentmarkerlist[i]);
    PyList_SetItem(holderlist,i, mlist);
  }    
  ii=PyString_FromString("generatedsegmentmarkerlist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);  
  /* Add triangle neighbor list */
  if (out.neighborlist != NULL) {
    listsize = out.numberoftriangles;
    holderlist = PyList_New(listsize);
    for(i=0; i<listsize;i++){
      PyObject *mlist = Py_BuildValue((char *)"(i,i,i)", 
				      out.neighborlist[i*3  ],
				      out.neighborlist [i*3+1],
				      out.neighborlist [i*3+2]);    
      PyList_SetItem(holderlist,i, mlist);
    }    
    ii=PyString_FromString("generatedtriangleneighborlist");
    PyDict_SetItem(holder, ii, holderlist);
    Py_DECREF(ii); Py_DECREF(holderlist);
  }   
  
  
  
  /* Free in/out structure memory */
  
  /* OUT */

  if(!out.pointlist){
    free(out.pointlist);  out.pointlist=NULL;
  }
  if(!out.pointmarkerlist){    
    free(out.pointmarkerlist); out.pointmarkerlist=NULL;
  }
  if(!out.pointattributelist){    
    free(out.pointattributelist); out.pointattributelist=NULL;
  }   
  if(!out.trianglelist){    
    free(out.trianglelist); out.trianglelist=NULL;
  }
  if(!out.triangleattributelist){    
    free(out.triangleattributelist); out.triangleattributelist=NULL;
  }
  if(!out.trianglearealist){    
    free(out.trianglearealist); out.trianglearealist=NULL;
  }
  if(!out.neighborlist){    
    free(out.neighborlist); out.neighborlist=NULL;
  }
  if(!out.segmentlist){
    free(out.segmentlist); out.segmentlist =NULL;
  }
  if(!out.segmentmarkerlist){
    free(out.segmentmarkerlist); out.segmentmarkerlist  =NULL;
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
  
  /* R = Py_BuildValue((char *)"O", holder); */
  R = Py_BuildValue((char *)"OO", holder, PyArray_Return(r_test));
  Py_DECREF(holder); /** This fixed a  memory problem ticket#189 */
  Py_DECREF(r_test);
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
