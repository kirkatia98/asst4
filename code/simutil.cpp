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
short *short_alloc(size_t n)
{
    return (short*) calloc(n, sizeof(short));
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

#if MPI
    array_t* m = (array_t*)malloc(sizeof(array_t));
    ok = ok && m != NULL;

    m->rsend = int_alloc(s->nprocess);
    ok = ok && m->rsend != NULL;

    m->rdisp = int_alloc(s->nprocess+1);
    ok = ok && m->rdisp != NULL;

    m->nsend = int_alloc(s->nprocess);
    ok = ok && m->rsend != NULL;

    m->ndisp = int_alloc(s->nprocess+1);
    ok = ok && m->rdisp != NULL;

    m->gsend = int_alloc(s->nprocess);
    ok = ok && m->rsend != NULL;

    m->gdisp = int_alloc(s->nprocess+1);
    ok = ok && m->rdisp != NULL;
    s->mpi = m;
#endif

    if (!ok) {
	outmsg("Couldn't allocate space for %d rats", nrat);
    }

    return s;
}

//DISPLACEMENTS AND SEND COUNTS
void send_disp(state_t* s)
{
    array_t* m = s->mpi;
    graph_t* g = s->g;
    int ideal = 4;

    //divide rats, for synchronous mode only
    int p;
    int per_process = s->nrat/ideal;
    int rem =  s->nrat%s->nprocess;
    int sum = 0;

    m->rdisp[0] = 0;
    for (p = 0; p < ideal; p++) {
        m->rsend[p] = per_process;

        if (rem > 0) {
            m->rsend[p]++;
            rem--;
        }
        sum += m->rsend[p];
        m->rdisp[p+1] = sum;
    }

    //divide nodes and gsums
    per_process = (g->tiles_per_side) /ideal;
    rem = (g->tiles_per_side) % s->nprocess;
    sum = 0;

    for (p = 0; p < ideal; p++) {
        m->nsend[p] = per_process;

        if (rem > 0) {
            m->nsend[p]++;
            rem--;
        }


        m->ndisp[p] = sum;
        m->nsend[p] *= g->nrow * g->tile_size;

        sum += m->nsend[p];
    }
    //last process may get truncated chunk
    m->ndisp[ideal] = g->nnode;
    s->my_nodes = m->ndisp[s->process_id + 1] - m->ndisp[s->process_id];


    //divide gsums based on nodes
    for (p = 0; p < ideal; p++) {

        int snode = m->gdisp[p];
        int enode = m->gsend[p+1];

        m->gdisp[p] = g->neighbor_start[snode];
        m->gsend[p] = m->gdisp[p] -  g->neighbor_start[enode];
    }

    for(p = ideal; p < s->nprocess; p++)
    {
        m->gdisp[p]= g->neighbor_start[g->nnode];
        m->ndisp[p] = g->nnode;
        m->rdisp[p] = s->nrat;
    }

    s->local_rat_count = int_alloc(s->my_nodes);
    s->delta = int_alloc(s->my_nodes);
    if(s->local_rat_count == NULL || s->delta == NULL)
    {
        outmsg("Couldn't allocate storage for local or delta\n");
    }
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
    array_t* m = s->mpi;
    free(m->rsend);
    free(m->rdisp);
    free(m->gdisp);
    free(m->gsend);
    free(m->nsend);
    free(m->ndisp);
    free(m);
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

