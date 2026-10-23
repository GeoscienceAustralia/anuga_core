#ifndef ANUGA_CONSTANTS_H
#define ANUGA_CONSTANTS_H

static const double TINY =1.0e-100;
static const double EPSILON = 1.0e-12;
static const double single_precision_epsilon = 1.0e-6;
static const double pi = 3.14159265358979;
static const double ETA_SMALL = 1.0e-15;

/* Exponents in the Manning friction terms, S = g n^2 |q| / h^(7/3), and the
   centroid average of three vertex values. Macros rather than static
   consts so they can seed firstprivate locals inside OpenMP target regions. */
#define ANUGA_ONE_THIRD    (1.0 / 3.0)
#define ANUGA_SEVEN_THIRDS (7.0 / 3.0)

/* Wet/dry limiting in the second-order extrapolation. The limiter weight
   hfactor ramps from 0 to 1 as the depth ratio (min or centroid depth over
   the max or centroid depth) rises from ANUGA_HFACTOR_B to ANUGA_HFACTOR_A:
       hfactor = clamp(c * ratio + d, 0, 1),
       c = 1 / (A - B),  d = 1 - c * A.
   Cells drier than the ramp are extrapolated first order. */
#define ANUGA_HFACTOR_A 0.3
#define ANUGA_HFACTOR_B 0.1

#endif
