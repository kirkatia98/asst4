#ifndef RUTIL_H
/*
  C library implementing math and timing functions for rat simulator.
  Provides key functions found in file rmath.py
*/

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/*
  Random number generator uses 64-bit arithmetic,
  but seed and random values guaranteed to fit in 32 bits
*/
typedef uint32_t random_t;

/* Standard parameters */
/* Default seed value */
#define DEFAULTSEED 618

/* Reinitialize seed based on list of seeds, where list has length len */
void reseed(random_t *seedp, random_t seed_list[], size_t len);

/* Generate double in range [0.0, upperlimit) */
double next_random_float(random_t *seedp, double upperlimit);

/* Compute weight function */
double mweight(double val);

#define RUTIL_H
#endif 
