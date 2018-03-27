#include "crun.h"



void outmsg(const char *fmt, ...) {
    va_list ap;
    bool got_newline = fmt[strlen(fmt)-1] == '\n';
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    if (!got_newline)
	fprintf(stderr, "\n");
}



/* Allocate n int's and zero them out. */
int *int_alloc(size_t n) {
    return (int *) calloc(n, sizeof(int));
}

/* Allocate n doubles's and zero them out. */
double *double_alloc(size_t n) {
    return (double *) calloc(n, sizeof(double));
}

/* Allocate n random number seeds and zero them out.  */
static random_t *rt_alloc(size_t n) {
    return (random_t *) calloc(n, sizeof(random_t));
}


/* Allocate simulation state */
state_t *new_rats(graph_t *g, int nrat, random_t global_seed) {
    int nnode = g->nnode;

    state_t *s = (state_t*) malloc(sizeof(state_t));

    if (s == NULL) {
	outmsg("Couldn't allocate storage for state\n");
	return NULL;
    }

    s->g = g;
    s->nrat = nrat;
    s->nprocess = 1;
    s->process_id = 0;
    s->global_seed = global_seed;
    s->load_factor = (double) nrat / nnode;

    /* Compute batch size as max(BATCH_FRACTION * R, sqrt(R)) */
    int rpct = (int) (BATCH_FRACTION * nrat);
    int sroot = (int) sqrt(nrat);
    if (rpct > sroot)
	s->batch_size = rpct;
    else
	s->batch_size = sroot;
    s->update_mode = UPDATE_BATCH;

    // Allocate data structures
    bool ok = true;
    s->rat_position = int_alloc(nrat);
    ok = ok && s->rat_position != NULL;

    s->next_rat_position = int_alloc(nrat);
    ok = ok && s->next_rat_position != NULL;

    s->rat_seed = rt_alloc(nrat);
    ok = ok && s->rat_seed != NULL;

    s->rat_count = int_alloc(nnode);
    ok = ok && s->rat_count != NULL;

    s->pre_computed = double_alloc(nrat+1);
    ok = ok && s->pre_computed != NULL;

#if MPI
    s->sendcounts = int_alloc(s->nprocess);
    ok = ok && s->sendcounts != NULL;

    s->disp = int_alloc(s->nprocess);
    ok = ok && s->disp != NULL;
#endif
    if (!ok) {
	outmsg("Couldn't allocate space for %d rats", nrat);
	return NULL;
    }

    return s;
}

/* Set seed values for the rats.  Maybe you could use multiple threads ... */
static void seed_rats(state_t *s) {
    random_t global_seed = s->global_seed;
    int nrat = s->nrat;
    int r;
    for (r = 0; r < nrat; r++) {
	random_t seeds[2];
	seeds[0] = global_seed;
	seeds[1] = r;
	reseed(&s->rat_seed[r], seeds, 2);
    }
}

/* See whether line of text is a comment */
static inline bool is_comment(char *s) {
    int i;
    int n = strlen(s);
    for (i = 0; i < n; i++) {
	char c = s[i];
	if (!isspace(c))
	    return c == '#';
    }
    return false;
}

/* Read in rat file */
state_t *read_rats(graph_t *g, FILE *infile, random_t global_seed) {
    char linebuf[MAXLINE];
    int r, nnode, nid, nrat;

    // Read header information
    while (fgets(linebuf, MAXLINE, infile) != NULL) {
	if (!is_comment(linebuf))
	    break;
    }
    if (sscanf(linebuf, "%d %d", &nnode, &nrat) != 2) {
	outmsg("ERROR. Malformed rat file header (line 1)\n");
	return false;
    }
    if (nnode != g->nnode) {
	outmsg("Graph contains %d nodes, but rat file has %d\n", g->nnode, nnode);
	return NULL;
    }
    
    state_t *s = new_rats(g, nrat, global_seed);


    for (r = 0; r < nrat; r++) {
	while (fgets(linebuf, MAXLINE, infile) != NULL) {
	    if (!is_comment(linebuf))
		break;
	}
	if (sscanf(linebuf, "%d", &nid) != 1) {
	    outmsg("Error in rat file.  Line %d\n", r+2);
	    return false;
	}
	if (nid < 0 || nid >= nnode) {
	    outmsg("ERROR.  Line %d.  Invalid node number %d\n", r+2, nid);
	    return false;
	}
	s->rat_position[r] = nid;
    }

    //calculate pre-computed mweights
    int i;
    for(i = 0; i <= s->nrat; i++)
    {
        s->pre_computed[i] = mweight((double) i/s->load_factor);
    }

    memset(s->rat_count, 0, nnode * sizeof(int));
    //for each rat, look at its position and increment the correct node
    int ri;
    for (ri = 0; ri < nrat; ri++) {
        s->rat_count[s->rat_position[ri]] ++;
    }

#if MPI
    if(g->tile_size == 1)
    {
        g->tile_size = 10; //set the tiles for divisibility
    }

    g->tiles_per_side = g->nrow/g->tile_size;

    //if nrow not divisible by tile size, extend the graph to be multiple
    // of tile size
    if(g->nrow % g->tile_size != 0)
    {
        g->tiles_per_side++;
        s->side_length = g->tiles_per_side * g->tile_size;
        int* new_rat_counts = int_alloc(s->side_length * s->side_length);

        int i, j;
        for(i = 0 ; i < g->nrow; i++)
        {
            for(j = 0; j < g->nrow; j++)
            {
                new_rat_counts[i*s->side_length +j] = s->rat_count[i*g->nrow + j];
            }
        }
        free(s->rat_count);
        s->rat_count = new_rat_counts;
    }
    else //use side length from now on
        s->side_length = g->nrow;
#endif

    seed_rats(s);
    outmsg("Loaded %d rats\n", nrat);
    return s;
}

/* print state of nodes */
void show(state_t *s, bool show_counts) {
    int nid;
    graph_t *g = s->g;
    printf("STEP %d %d\n", g->nnode, s->nrat);
    if (show_counts) {
        int i, j;


        for(i = 0 ; i < g->nrow; i++)
        {
            for(j = 0; j < g->nrow; j++)
            {
                printf("%d\n", s->rat_count[i*s->side_length + j]);
            }
        }
    }
    printf("END\n");
}

/* Print final output */
void done() {
    printf("DONE\n");
}

