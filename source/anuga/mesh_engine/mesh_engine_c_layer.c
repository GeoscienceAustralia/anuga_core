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

#include "Numeric/arrayobject.h"


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
  REAL Attr;
  int i, j, iatt, n, write_here;
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
  
  int *points_connected;
  int *lone_verts;
  
  /* used for testing numeric arrays*/
  int n_test;     
  if(!PyArg_ParseTuple(args,(char *)"OOOOOOOO",&pointlist,&seglist,&holelist,&regionlist,&pointattributelist,&segmarkerlist,&mode,&test)){
    return NULL;
  }
  /*if(!PyList_Check(hull)){
    PyErr_SetString(PyExc_TypeError, 
		    "incorrect first argument for 'genMesh': list required.");
    return NULL;
    } */
     
  /* printf("starting\n");   
  n_test=test -> dimensions[0];
  printf("test -> dimensions[0] %i\n", n_test); 
  printf("test -> data %f\n", ((double *) test -> data)[0]); /*returns 1.0*/
     
  /* Initialize  points 
  in.numberofpoints =  PyList_Size(hull); 
  in.pointlist = (REAL *)malloc(in.numberofpoints*2*sizeof(REAL));
  in.pointmarkerlist = (int *)NULL; /* malloc(in.numberofpoints*sizeof(int)); */
  
  /* Initialize  points */
  in.numberofpoints =  pointlist-> dimensions[0]; 
  in.pointlist = (double *) pointlist -> data;
  in.pointmarkerlist = (int *)NULL; 
  
  /* Initialize and fill up region list */
  in.numberofregions = regionlist-> dimensions[0]; 
  in.regionlist = (double *) regionlist -> data;
  
  /* Initialize and fill up region list 
  in.numberofregions = PyList_Size(regionlist);
  
  /*printf(" numberofregions %i\n",in.numberofregions );
  in.regionlist = (REAL *)malloc(in.numberofregions*4*sizeof(REAL));
  for(i=0; i<in.numberofregions;i++){
    elems = PyList_GetItem(regionlist,i);
    tuplesize = PySequence_Length(elems);
    /*printf("tuplesize %i\n",tuplesize);
    if (tuplesize>0){
      in.regionlist[4*i] = PyFloat_AsDouble(PyTuple_GetItem(elems,0));
    } else {
      PyErr_SetString(PyExc_TypeError,
		      "data sent to trianglemodule error: A region has no x value.");
      return NULL;
    }
    if (tuplesize>1){
      in.regionlist[4*i+1] = PyFloat_AsDouble(PyTuple_GetItem(elems,1));
    } else {
      PyErr_SetString(PyExc_TypeError,
		    "data sent to trianglemodule error: A region has no y value.");
      return NULL;
    }
    if (tuplesize>2){
      /* 2 is attribute
      in.regionlist[4*i+2] = PyFloat_AsDouble(PyTuple_GetItem(elems,2));
    } else {
      PyErr_SetString(PyExc_TypeError,
		    "data sent to trianglemodule error: A region has no attribute value.");
      return NULL;
    }
    if (tuplesize>3){
      in.regionlist[4*i+3] = PyFloat_AsDouble(PyTuple_GetItem(elems,3));
    } else {
      /* 0.0 is interpreted as no local max area by triangle  
      in.regionlist[4*i+3] = (REAL)0.0;
    } 
  }*/
  
  /*Initialize segments and segment markers */
  in.numberofsegments = seglist -> dimensions[0]; 
  in.segmentlist = (int *) seglist -> data;
  in.segmentmarkerlist = (int *) segmarkerlist -> data;
  
  /*Initialize segments 
  in.numberofsegments = PyList_Size(seglist);
  in.segmentlist = (int *)malloc(in.numberofsegments*2*sizeof(int));
  in.segmentmarkerlist = (int *)malloc(in.numberofsegments*sizeof(int)); */

    
  /*Initialize triangles */
  in.numberoftriangles = 0;
  in.trianglelist = (int *)NULL;
  in.numberoftriangleattributes = 0;  
  in.triangleattributelist  = (REAL *)NULL; 
     
  /*Initialize holes
  in.numberofholes = PyList_Size(holelist);
  if (in.numberofholes != 0) {
    in.holelist = (REAL *)malloc(in.numberofholes*2*sizeof(REAL));
  } else {
    in.holelist = (REAL *)NULL;
  }      */
     
  /*Initialize holes */
  in.numberofholes = holelist  -> dimensions[0];
  if (in.numberofholes != 0) {
    in.holelist =  (double *) holelist -> data; 
  } else {
    in.holelist = (REAL *)NULL;
  }     
  
      
  /* Fill up segment list 
  index = 0;
  if (in.numberofsegments != PyList_Size(segmarkerlist)) {
    PyErr_SetString(PyExc_TypeError,
		    "data sent to trianglemodule error:number of segments doesn't = number of segment attributes.");
    return NULL;
  }
    
  for(i = 0; i < in.numberofsegments; i++) { */
    /*
    pair =  PyList_GetItem(seglist,i);
    if (2 != PySequence_Length(pair)) {
      PyErr_SetString(PyExc_TypeError,
		      "data sent to trianglemodule error:Found a segment without 2 vertices.");
      return NULL;
    }
    a = (int)PyInt_AsLong(PyTuple_GetItem(pair,0));
    in.segmentlist[index] = a;
    index++;
    b = (int)PyInt_AsLong(PyTuple_GetItem(pair,1));
    in.segmentlist[index] = b;
    index++;
    */
    /* Fill up segment marker list 
    marker =  (int)PyInt_AsLong(PyList_GetItem(segmarkerlist,i));
    in.segmentmarkerlist[i] = marker; 
  } */
  /* Fill up hole list 
  index = 0;
  for(i = 0; i < in.numberofholes; i++) {
    pair =  PyList_GetItem(holelist,i);
    if (2 != PySequence_Length(pair)) {
      PyErr_SetString(PyExc_TypeError,
		      "data sent to trianglemodule error:Found a pair without 2 vertices.");
      return NULL;
    }
    x = PyFloat_AsDouble(PyTuple_GetItem(pair,0));
    in.holelist[index] = x;
    index++;
    y = PyFloat_AsDouble(PyTuple_GetItem(pair,1));
    in.holelist[index] = y;
    index++;
  } */
     
  /* Fill up point list
  index = 0;
  for(i = 0; i < in.numberofpoints; i++) {
    pair =  PyList_GetItem(hull,i);
    if (2 != PySequence_Length(pair)) {
      PyErr_SetString(PyExc_TypeError,
		      "data sent to trianglemodule error:Found a vertex without an x,y values.");
      return NULL;
    }
    x = PyFloat_AsDouble(PyTuple_GetItem(pair,0));
    in.pointlist[index] = x;
    index++;
    y = PyFloat_AsDouble(PyTuple_GetItem(pair,1));
    in.pointlist[index] = y;
    index++;
  }  */
  
  /* for(i = 0; i < in.numberofpoints * 2; i++) {
    printf("in.pointlist[index] %f\n", in.pointlist[i]);  
    } */
  /*
  for(i = 0; i < in.numberofpoints * 2; i++) {
    printf("in_test.pointlist[index] %f\n", in_test.pointlist[i]);  
  }
  
  printf("in_test.numberofpoints %i\n", in_test.numberofpoints); 
  printf("test -> data %f\n", ((double *) test -> data)[1,1]); /*returns 1.0*/

  
  /* Initialize  points
  in.numberofpoints = in_test.numberofpoints;
  in.pointlist = in_test.pointlist;
  in.pointmarkerlist =  in_test.pointmarkerlist;  */
  
  /* Initialize  point attributes */
/*      firstPointAttr = PyList_GetItem(pointattributelist,0); */
/*      if (firstPointAttr != NULL) { */
/*        //in.numberofpointattributes = PyList_Size(firstPointAttr); */
/*        printf("in.numberofpointattributes =%i",PyList_Size(firstPointAttr)); */
/*        //printf("in.numberofpointattributes =%i",in.numberofpointattributes); */
/*      } else { */
/*        printf("(firstPointAttr is NULL\n"); */
/*        in.pointattributelist = NULL; */
/*        in.numberofpointattributes = 0; */
/*      } */

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
   
  /* 
  n = PyList_Size(pointattributelist);
  if (n < 0) {
    in.numberofpointattributes = 0; 
    in.pointattributelist = (REAL *)NULL;
  } else {
    /* Get the # of attributes based on the size of the first point
    firstPointAttr = PyList_GetItem(pointattributelist,0);
    if (firstPointAttr == NULL) {
      in.numberofpointattributes = 0; 
      in.pointattributelist = (REAL *)NULL;
    } else {
      in.numberofpointattributes = PyList_Size(firstPointAttr);
      in.pointattributelist = (REAL *)malloc(in.numberofpoints
					     *in.numberofpointattributes
					     *sizeof(REAL)); 
     
     
      /* fill the point attribute list 
      if (in.numberofpointattributes != 0) {
	for(j = 0; j < in.numberofpoints; j++){
	  Attrlist = PyList_GetItem(pointattributelist,j);
	  if (in.numberofpointattributes != PySequence_Length(Attrlist)) {
	    PyErr_SetString(PyExc_TypeError,
			    "data sent to trianglemodule error:List size of attributes is different from vertex to vertex.");
	    return NULL;
	  }
	  for(i = 0; i < in.numberofpointattributes; i++){
	    Attr = PyFloat_AsDouble(PyList_GetItem(Attrlist,i));
	    in.pointattributelist[in.numberofpointattributes*j+i]= (REAL)Attr;
	    /* 	 printf("i =%i\n",i);  
	    /* 	 printf("j =%i\n",j);  
	    /* 	 printf("in.pointattributelist[4*j+i]
		 =%f\n",in.pointattributelist[4*j+i]); 
	  }
	}
      } 
    }
    } */
    
     
  /* this half from triangulate*/
  /* set up the switch arguments to the triangulation routine */

  /* DEBUGGING
  printf(" ***  hello world\n" );  
  
  printf ("numberofholess %i\n", in.numberofholes);
  for(i = 0; i < in.numberofholes; i++) {
    printf("(%f,%f)\n",in.holelist[i* 2+ 0],in.holelist[i* 2+ 1]);
  }
     
      printf ("numberoftriangles %i\n", in.numberoftriangles);
      for(i = 0; i < in.numberoftriangles; i++) {
      printf("(%i,%i,%i)",in.trianglelist[i* 3+ 0],in.trianglelist[i* 3+ 1], in.trianglelist[i* 3+ 2]);
      }
 
  printf(" ***  see ya world\n" );       */
 
      
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
    ------- Pass point numbers,coordinates and neighbors back to Python ------
    we send back a dictionary:                                               
    { index : [ coordinates, [connections], Attribute ] } 
  */
  holder = PyDict_New();
     
  /* list of int's, used to keep track of which verts are connected to
     triangles. */
  points_connected = (int *)malloc(out.numberofpoints*sizeof(int));
  /* lone_verts = (int *)malloc(out.numberofpoints*sizeof(int)); */ 
  
  /* Initialise lone vert list */
  for(i=0; i<out.numberofpoints;i++){
    points_connected[i] = 0;
    /* lone_verts[i] = 0; */
  }  
  
  /* Add triangle list */
  listsize = out.numberoftriangles;
  holderlist = PyList_New(listsize);
  for(i=0; i<listsize;i++){
    PyObject *mlist = Py_BuildValue((char *)"(i,i,i)", 
				    out.trianglelist[i*3  ], out.trianglelist [i*3+1], out.trianglelist [i*3+2]);    
    PyList_SetItem(holderlist,i, mlist);
    /* printf(" A vert index %i\n",out.trianglelist[i*3] );
    printf(" A vert index %i\n",out.trianglelist[i*3+1] );
    printf(" A vert index %i\n",out.trianglelist[i*3+2] ); */
    points_connected[out.trianglelist[i*3]] = 1;
    points_connected[out.trianglelist[i*3+1]] = 1;
    points_connected[out.trianglelist[i*3+2]] = 1;
    /* lone_verts[out.trianglelist[i*3]] = 1;
    lone_verts[out.trianglelist[i*3+1]] = 1;
    lone_verts[out.trianglelist[i*3+2]] = 1; */
    /*  printf(" Add triangle list \n");*/
  }    
  ii=PyString_FromString("generatedtrianglelist");
  PyDict_SetItem(holder, ii, holderlist); Py_DECREF(ii); Py_DECREF(holderlist);
  
  /* convert the points_connected vector from a true(1) false(0) vector, where
     index is the vert, to a vector of the lone verts, at the beggining
     of the vector. */
  write_here = 0;   
  for(i=0; i<out.numberofpoints;i++){
    /* lone_verts[i-write_here] = lone_verts[i]; */
    if (0 == points_connected[i]) {
      points_connected[write_here] = i;
      write_here ++;
    }
  }   
  /* printf(" ******************** \n" );
  for(i=0; i<write_here;i++){
    printf(" A vert index %i\n",points_connected[i] );
    } */
  
  listsize = write_here;
  holderlist = PyList_New(listsize);
  for(i=0; i<listsize;i++){
   PyObject *mlist = Py_BuildValue((char *)"i", points_connected[i]);    
    PyList_SetItem(holderlist,i, mlist);
  }
  ii=PyString_FromString("lonepointlist");
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
  if(!points_connected ){
    free(points_connected ); points_connected =NULL;
  }
    
  
  
  /*  IN */
  /* These are all numeric structures, so I shouldn't free them,
  but should I decrease the ref?
  Based on examples of this stuff I don't have to decrease
  the ref
    
  if(!in.pointlist){
    free(in.pointlist);  in.pointlist   =NULL;             
  }
  if(!in.pointattributelist){
    free(in.pointattributelist); in.pointattributelist =NULL;       
  }
  if(!in.pointmarkerlist){
    free(in.pointmarkerlist); in.pointmarkerlist    =NULL;        
  }
  if(!in.segmentlist){
    free(in.segmentlist);  in.segmentlist    =NULL;         
  }
  if(!in.segmentmarkerlist){
    free(in.segmentmarkerlist); in.segmentmarkerlist    =NULL;   
  }
  if(!in.regionlist){
    free(in.regionlist); in.regionlist  =NULL;            
  }
  if(!in.holelist){
    free(in.holelist); in.holelist=NULL;
    } */

  return Py_BuildValue((char *)"O", holder);
}

/* end of my function*/


   
    
static PyMethodDef triang_methods[] = {
  {(char *)"genMesh", triang_genMesh, 1, (char *)"Passes hull points to triangle.All arguments are lists.(<hull>,<seglist>,<holelist>,<regionlist>)-><hull>"}, 
  {NULL,NULL}
};

void initmesh_engine_c_layer(){
  Py_InitModule((char *)"mesh_engine_c_layer",triang_methods);
}    
