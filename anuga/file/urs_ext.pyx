#cython: wraparound=False, boundscheck=False, cdivision=True, profile=False, nonecheck=False, overflowcheck=False, cdivision_warnings=False, unraisable_tracebacks=False
import cython
from libc.stdlib cimport malloc, free
# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern from "urs.c":
  float** _read_mux2(int numSrc, char** muxFileNameArray, float* weights, double* params, int* number_of_stations, long* permutation, int verbose)

def read_mux2(int numSrc,\
              list filenames,\
              np.ndarray[double, ndim=1, mode="c"] pyweights not None,\
              np.ndarray[double, ndim=1, mode="c"] file_params not None,\
              np.ndarray[long, ndim=1, mode="c"] permutation not None,\
              int verbose):
  # Make sure filenames are cast as str().encode(). This will work in both Python2 and Python3 (Ole)	     

  cdef char** muxFileNameArray
  cdef float** cdata
  cdef float* weights
  cdef int dimensions[2]
  cdef int total_number_of_stations
  cdef int number_of_selected_stations
  cdef int nt
  cdef double dt
  cdef int i, j, start_tstep, finish_tstep, it, time, num_ts
  cdef int POFFSET = 5
  cdef np.ndarray[double, ndim=2, mode="c"] pydata
  cdef np.ndarray[long, ndim=1, mode="c"] tmp_perm

  tmp_perm = np.array([0],int)

  assert len(filenames) > 0, "empty lists not allowed"

  assert len(filenames) == pyweights.shape[0], "Must specify one weight for each filename"

  muxFileNameArray = <char** > malloc(numSrc * sizeof(char*))
  assert muxFileNameArray != NULL, "ERROR: Memory for muxFileNameArray could not be allocated."

  for i in xrange(numSrc):
    #assert isinstance(filenames[i], str), "filename not a string"  # Nor should it be ;-)
    #print(filenames[i], type(filenames[i]))
    muxFileNameArray[i] = filenames[i]
  
  weights = <float* > malloc(numSrc * sizeof(float))
  assert weights != NULL, "ERROR: Memory for weights could not be allocated."

  for i in xrange(numSrc):
    weights[i] = pyweights[i]
  
  number_of_selected_stations = permutation.shape[0]

  if number_of_selected_stations == 0:
    permutation = tmp_perm

  cdata = _read_mux2(numSrc,\
		     muxFileNameArray,\
                     weights,\
                     &file_params[0],\
                     &number_of_selected_stations,\
                     &permutation[0],\
                     verbose)
  
  assert cdata != NULL, "No STS_DATA returned"

  total_number_of_stations = int(file_params[0])
  dt = file_params[1]
  nt = int(file_params[2])

  # Find min and max start times of all gauges
  start_tstep = nt + 1
  finish_tstep = -1
  for i in xrange(number_of_selected_stations):
    if int(cdata[i][nt + 3]) < start_tstep:
      start_tstep = int(cdata[i][nt + 3])
    if int(cdata[i][nt + 4]) > finish_tstep:
      finish_tstep = int(cdata[i][nt + 4])
  
  if start_tstep > nt or finish_tstep < 0:
    print ("ERROR: Gauge data has incorrect start and finish times:")
    print ("   start_tstep = %d, max_number_of_steps = %d" % (start_tstep, nt))
    print ("   finish_tstep = %d, min_number_of_steps = %d" % (finish_tstep, 0))

    assert start_tstep <= nt and finish_tstep >= 0, "Incorrect start and finish times"

    free(weights)
    free(muxFileNameArray)

    for i in xrange(number_of_selected_stations):
      free(cdata[i])
    free(cdata)

  if start_tstep >= finish_tstep:
    assert start_tstep < finish_tstep, "ERROR: Gauge data has non-positive length"
    free(weights)
    free(muxFileNameArray)
    for i in xrange(number_of_selected_stations):
      free(cdata[i])
    free(cdata)
  
  num_ts = finish_tstep - start_tstep + 1
  dimensions[0] = number_of_selected_stations
  dimensions[1] = num_ts + POFFSET

  # Each gauge begins and ends recording at different times. When a gauge is
  # not recording but at least one other gauge is.
  # Pad the non recording gauge array with zeros.
  pydata = np.zeros((dimensions[0],dimensions[1]),dtype=float)
  for i in xrange(number_of_selected_stations):
    time = 0
    for it in xrange(finish_tstep):
      if it + 1 >= start_tstep and it + 1 <= finish_tstep:
        if it + 1 <= int(cdata[i][nt + 4]):
          pydata[i,time] = cdata[i][it]
        time += 1
    for j in xrange(POFFSET):
      pydata[i,num_ts + j] = cdata[i][nt + j]
  
  free(weights)
  free(muxFileNameArray)
  for i in xrange(number_of_selected_stations):
    free(cdata[i])
  free(cdata)

  return pydata

 
