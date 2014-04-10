struct domain* get_python_domain(struct domain *D, PyObject *domain, double timestep) {
  PyArrayObject *neighbours, *neighbour_vertices, *normals, *areas, *max_speed_array;    
    
  D->timestep = timestep;
  D->params->cfl = get_python_double(domain, "CFL");
  D->params->epsilon = get_python_double(domain, "epsilon");
  
  neighbours = get_consecutive_array(domain, "neighbours");
  D->neighbours = (long *) neighbours->data;
  neighbour_vertices = get_consecutive_array(domain, "neighbour_vertices"); 
  D->neighbour_vertices = (long *) neighbour_vertices->data;
  normals = get_consecutive_array(domain, "normals");
  D->normals = (double *) normals->data;
  areas = get_consecutive_array(domain, "areas");    
  D->areas = (double *) areas->data;
  max_speed_array = get_consecutive_array(domain, "max_speed_array");
  D->max_speed_array = (double *) max_speed_array->data;
  
  D->number_of_elements = get_python_integer(domain, "number_of_elements");

  Py_DECREF(neighbours);
  Py_DECREF(neighbour_vertices);
  Py_DECREF(normals);
  Py_DECREF(areas);
  Py_DECREF(max_speed_array);

  return D;
}
