#ifndef CRUN_H

/* Defining variable MPI enables use of MPI primitives */
#ifndef MPI
#define MPI 0
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#if MPI
#include <mpi.h>
#endif


/* Optionally enable debugging routines */
#ifndef DEBUG
#define DEBUG 0
#endif

#include "rutil.h"
#include "cycletimer.h"

/*
  Definitions of all constant parameters.  This would be a good place
  to define any constants and options that you use to tune performance
*/

/* What is the maximum line length for reading files */
#define MAXLINE 1024

/* What is the batch size as a fraction of the number of rats */
#define BATCH_FRACTION 0.02


/* Update modes */
typedef enum { UPDATE_SYNCHRONOUS, UPDATE_BATCH, UPDATE_RAT } update_t;

/* All information needed for graphrat simulation */

/* Parameter abbreviations
   N = number of nodes
   M = number of edge
   R = number of rats
   B = batch size
   T = number of threads
 */


/* Representation of graph */
typedef struct {
    /* General parameters */
    int nnode;
    int nedge;
    int nrow;  /* == sqrt(nnode) */
    int tile_max;  /* Maximum number of consecutive rows having non-grid connections */

    /* Graph structure representation */
    // Adjacency lists.  Includes self edge. Length=M+N.  Combined into single vector
    int *neighbor;
    // Starting index for each adjacency list.  Length=N+1
    int *neighbor_start;

} graph_t;

/* Representation of simulation state */
typedef struct {
    graph_t *graph;

    int nrat;
    /* MPI processes & process id */
    int nprocess;
    int process_id;

    /* Random seed controlling simulation */
    random_t global_seed;

    /* State representation */
    // Node Id for each rat.  Length=R
    int *rat_position;
    // Next node Id for each rat.  Length=R
    int *next_rat_position;
    // Rat seeds.  Length = R
    random_t *rat_seed;

    /* Redundant encodings to speed computation */
    // Count of number of rats at each node.  Length = N.
    int *rat_count;
    /* Computed parameters */
    double load_factor;  // nrat/nnnode
    update_t update_mode; 
    int batch_size;   // Batch size for batch mode

} state_t;
    

/*** Functions in graph.c. ***/
graph_t *new_graph(int nnode, int nedge, int tile_max);

void free_graph();

graph_t *read_graph(FILE *gfile);

#if DEBUG
void show_graph(graph_t *g);
#endif


/*** Functions in simutil.c ***/
/* Print message on stderr */
void outmsg(char *fmt, ...);


/* Allocate and zero arrays of int/double */
int *int_alloc(size_t n);
double *double_alloc(size_t n);


/* Read rat file and initialize simulation state */
state_t *read_rats(graph_t *g, FILE *infile, random_t global_seed);


/* Generate done message from simulator */
void done();

/* Print state of simulation */
/* show_counts indicates whether to include counts of rats for each node */
void show(state_t *s, bool show_counts);

/*** Functions in sim.c ***/

/* Run simulation */
void simulate(state_t *s, int count, update_t update_mode, int dinterval, bool display);


#define CRUN_H
#endif /* CRUN_H */
