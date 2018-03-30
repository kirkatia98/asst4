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
    return (int *) malloc(n* sizeof(int));
}
short *short_alloc(size_t n)
{
    return (short*) malloc(n* sizeof(short));
}

/* Allocate n doubles's and zero them out. */
double *double_alloc(size_t n) {
    return (double *) malloc(n* sizeof(double));
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

#if MPI
    MPI_Comm_size(MPI_COMM_WORLD, &s->nprocess);
#else
    s->nprocess = 1;
#endif

    s->g = g;
    s->nrat = nrat;

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
    s->rat_position = short_alloc(nrat);
    ok = ok && s->rat_position != NULL;

    s->next_position = short_alloc(nrat);
    ok = ok && s->next_position != NULL;

    s->rat_seed = rt_alloc(nrat);
    ok = ok && s->rat_seed != NULL;

    s->delta = int_alloc(s->my_nodes);
    ok = ok && s->delta != NULL;

#if MPI
    s->send = int_alloc(s->nprocess+10);
    ok = ok && s->send != NULL;

    s->disp = int_alloc(s->nprocess + 11);
    ok = ok && s->disp != NULL;
#endif

    if (!ok) {
	outmsg("Couldn't allocate space for %d rats", nrat);
	exit;
    }

    return s;
}

void free_state(state_t *s)
{
    if(s->process_id == 0)
    {
        free(s->rat_count);
        free(s->pre_computed);
    }
    free(s->local_rat_count);
    free(s->rat_position);
    free(s->next_position);
    free(s->rat_seed);
    free(s->delta);
#if MPI
    free(s->send);
    free(s->disp);
#endif
    free(s);
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

    //allocate for Master process
    state_t *s = new_rats(g, nrat, global_seed);

    s->pre_computed = double_alloc(nrat+1);
    s->rat_count = int_alloc(nnode);

    if(s->rat_count == NULL || s->pre_computed == NULL)
        outmsg("Master process couldn't allocate rat count or pre-computed", nrat);


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
        g->tile_size = (10 < g->nrow ? 10 : g->nrow ); //set the tiles for
        // divisibility
    }

    g->tiles_per_side = g->nrow/g->tile_size;
    if(g->nrow%g->tile_size != 0)
    {
        g->tiles_per_side++;
    }
#endif

    seed_rats(s);
    outmsg("Loaded %d rats\n", nrat);
    return s;
}

/* print state of nodes */
void show(state_t *s, bool show_counts) {
    graph_t *g = s->g;
    printf("STEP %d %d\n", g->nnode, s->nrat);
    if (show_counts) {
        int i;
        for(i = 0 ; i < g->nnode; i++)
        {
            printf("%d\n", s->rat_count[i]);
        }
    }
    printf("END\n");
}

/* Print final output */
void done() {
    printf("DONE\n");
}

