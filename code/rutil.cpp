/*
  C library implementing math functions for rat simulator.
  Provides same functionality as Python code in rutil.py
*/

#include "rutil.h"

/* Standard parameters */
#define GROUPSIZE 2147483647
#define MVAL  48271
#define VVAL  16807
#define INITSEED  418


static inline random_t rnext(random_t *seedp, random_t x) {
    uint64_t s = (uint64_t) *seedp;
    uint64_t xlong = (uint64_t) x;
    random_t val = ((xlong+1) * VVAL + s * MVAL) % GROUPSIZE;
    *seedp = (random_t) val;
    return val;
}

/* Reinitialize seed based on list of seeds, where list has length len */
void reseed(random_t *seedp, random_t seed_list[], size_t len) {
    *seedp = INITSEED;
    size_t i;
    for (i = 0; i < len; i++)
	rnext(seedp, seed_list[i]);
}

/* Generate double in range [0.0, upperlimit) */
double next_random_float(random_t *seedp, double upperlimit) {
    random_t val = rnext(seedp, 0);
    return ((double) val / (double) GROUPSIZE) * upperlimit;
}

/* Parameters for computing weights that guide next-move selection */
#define COEFF 0.5
#define OPTVAL 1.5

double mweight(double val) {
    double arg = 1.0 + COEFF * (val - OPTVAL);
    double lg = log(arg) * M_LOG2E;
    double denom = 1.0 + lg * lg;
    return 1.0/denom;
}
