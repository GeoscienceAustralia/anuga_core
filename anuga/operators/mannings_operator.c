// Python - C extension module for shallow_water.py
//
// To compile (Python2.3):
//  gcc -c domain_ext.c -I/usr/include/python2.3 -o domain_ext.o -Wall -O
//  gcc -shared domain_ext.o  -o domain_ext.so
//
// or use python compile.py
//
// See the module shallow_water_domain.py for more documentation on 
// how to use this module
//
//
// Ole Nielsen, GA 2004


#include "math.h"
#include "util_ext.h"
#include <stdio.h>

const double pi = 3.14159265358979;

void _manning_friction_flat(double g, double eps, int N,
               double* w, double* zv,
               double* uh, double* vh,
               double* eta, double* xmom, double* ymom) {

  int k, k3;
  double S, h, z, z0, z1, z2;

  for (k=0; k<N; k++) {
    if (eta[k] > eps) {
      k3 = 3*k;
      // Get bathymetry
      z0 = zv[k3 + 0];
      z1 = zv[k3 + 1];
      z2 = zv[k3 + 2];
      z = (z0+z1+z2)/3.0;
      h = w[k]-z;
      if (h >= eps) {
        S = -g * eta[k]*eta[k] * sqrt((uh[k]*uh[k] + vh[k]*vh[k]));
        S /= pow(h, 7.0/3);      //Expensive (on Ole's home computer)
        //S /= exp((7.0/3.0)*log(h));      //seems to save about 15% over manning_friction
        //S /= h*h*(1 + h/3.0 - h*h/9.0); //FIXME: Could use a Taylor expansion


        //Update momentum
        xmom[k] += S*uh[k];
        ymom[k] += S*vh[k];
      }
    }
  }
}


void _manning_friction_sloped(double g, double eps, int N,
               double* x, double* w, double* zv,
               double* uh, double* vh,
               double* eta, double* xmom_update, double* ymom_update) {

  int k, k3, k6;
  double S, h, z, z0, z1, z2, zs, zx, zy;
  double x0,y0,x1,y1,x2,y2;

  for (k=0; k<N; k++) {
    if (eta[k] > eps) {
      k3 = 3*k;
      // Get bathymetry
      z0 = zv[k3 + 0];
      z1 = zv[k3 + 1];
      z2 = zv[k3 + 2];

      // Compute bed slope
      k6 = 6*k;  // base index
  
      x0 = x[k6 + 0];
      y0 = x[k6 + 1];
      x1 = x[k6 + 2];
      y1 = x[k6 + 3];
      x2 = x[k6 + 4];
      y2 = x[k6 + 5];

      _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);

      zs = sqrt(1.0 + zx*zx + zy*zy);
      z = (z0+z1+z2)/3.0;
      h = w[k]-z;
      if (h >= eps) {
        S = -g * eta[k]*eta[k] * zs * sqrt((uh[k]*uh[k] + vh[k]*vh[k]));
        S /= pow(h, 7.0/3);      //Expensive (on Ole's home computer)
        //S /= exp((7.0/3.0)*log(h));      //seems to save about 15% over manning_friction
        //S /= h*h*(1 + h/3.0 - h*h/9.0); //FIXME: Could use a Taylor expansion


        //Update momentum
        xmom_update[k] += S*uh[k];
        ymom_update[k] += S*vh[k];
      }
    }
  }
}



void _chezy_friction(double g, double eps, int N,
               double* x, double* w, double* zv,
               double* uh, double* vh,
               double* chezy, double* xmom_update, double* ymom_update) {

  int k, k3, k6;
  double S, h, z, z0, z1, z2, zs, zx, zy;
  double x0,y0,x1,y1,x2,y2;

  for (k=0; k<N; k++) {
    if (chezy[k] > eps) {
      k3 = 3*k;
      // Get bathymetry
      z0 = zv[k3 + 0];
      z1 = zv[k3 + 1];
      z2 = zv[k3 + 2];

      // Compute bed slope
      k6 = 6*k;  // base index

      x0 = x[k6 + 0];
      y0 = x[k6 + 1];
      x1 = x[k6 + 2];
      y1 = x[k6 + 3];
      x2 = x[k6 + 4];
      y2 = x[k6 + 5];

      _gradient(x0, y0, x1, y1, x2, y2, z0, z1, z2, &zx, &zy);

      zs = sqrt(1.0 + zx*zx + zy*zy);
      z = (z0+z1+z2)/3.0;
      h = w[k]-z;
      if (h >= eps) {
        S = -g * chezy[k] * zs * sqrt((uh[k]*uh[k] + vh[k]*vh[k]));
        S /= pow(h, 2.0);

        //Update momentum
        xmom_update[k] += S*uh[k];
        ymom_update[k] += S*vh[k];
      }
    }
  }
}
