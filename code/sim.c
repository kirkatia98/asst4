#include "crun.h"

/* Compute weight for node nid */
static inline double compute_weight(state_t *s, int nid) {
    int count = s->rat_count[nid];
    return mweight((double) count/s->load_factor);
}

/* Compute sum of weights in region of nid */
static inline double compute_sum_weight(state_t *s, int nid) {
    graph_t *g = s->graph;
    double sum = 0.0;
    int eid;
    int eid_start = g->neighbor_start[nid];
    int eid_end  = g->neighbor_start[nid+1];
    int *neighbor = g->neighbor;
    for (eid = eid_start; eid < eid_end; eid++) {
	sum += compute_weight(s, neighbor[eid]);

    }
    return sum;
}

/** DEBUGGING CODE **/
#if DEBUG
/* Handy place to set a breakpoint */
static void panic(char *msg) {
    outmsg("PANIC: %s\n", msg);
}

static void show_weights(state_t *s) {
    int nid, eid;
    graph_t *g = s->graph;
    int nnode = g->nnode;
    int *neighbor = g->neighbor;
    outmsg("Weights\n");
    for (nid = 0; nid < nnode; nid++) {
	int eid_start = g->neighbor_start[nid];
	int eid_end  = g->neighbor_start[nid+1];
	track = true;
	outmsg("In show_weights:\n");
	compute_sum_weight(s, nid);
	track = false;
	outmsg("%d: [sum = %.3f]", nid, compute_sum_weight(s, nid));
	for (eid = eid_start; eid < eid_end; eid++) {
	    outmsg(" %.3f", compute_weight(s, neighbor[eid]));
	}
	outmsg("\n");
    }
}
#endif


/* Recompute all node counts according to rat population */
static inline void take_census(state_t *s) {
    int *rat_position = s->rat_position;
    int *rat_count = s->rat_count;
    int nrat = s->nrat;
    int ri;
    graph_t *g = s->graph;
    int nnode = g->nnode;
    memset(rat_count, 0, nnode * sizeof(int));
    for (ri = 0; ri < nrat; ri++) {
	rat_count[rat_position[ri]] ++;
    }
}



/*
  Given list of integer counts, generate real-valued weights
  and use these to flip random coin returning value between 0 and len-1
*/
static inline int next_random_move(state_t *s, int r) {
    int nid = s->rat_position[r];
    int nnid = -1;
    random_t *seedp = &s->rat_seed[r];
    double tsum = compute_sum_weight(s, nid);
    graph_t *g = s->graph;
    int eid;
    
    double val = next_random_float(seedp, tsum);

    double psum = 0.0;
    for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++) {
	psum += compute_weight(s, g->neighbor[eid]);
	if (val < psum) {
	    nnid = g->neighbor[eid];
	    break;
	}
    }
    if (nnid == -1) {
	/* Shouldn't get here */
	int degree = g->neighbor_start[nid+1] - g->neighbor_start[nid];
	outmsg("Internal error.  next_random_move.  Didn't find valid move.  Node %d. Degree = %d, Target = %.2f/%.2f.  Limit = %.2f\n",
	       nid, degree, val, tsum, psum);
	nnid = 0;
    }
    return nnid;
}







static void process_batch(state_t *s, int bstart, int bcount) {
    int rid;
    for (rid = bstart; rid < bstart+bcount; rid++)
	// Determine where rat will go
	s->next_rat_position[rid] = next_random_move(s, rid);
    for (rid = bstart; rid < bstart+bcount; rid++) {
	int onid = s->rat_position[rid];
	int nnid = s->next_rat_position[rid];
	s->rat_count[onid]--;
	s->rat_count[nnid]++;
	s->rat_position[rid] = nnid;
    }
}

static void run_step(state_t *s, int batch_size) {
    int b, bcount;
    for (b = 0; b < s->nrat; b += batch_size) {
	int rest = s->nrat - b;
	bcount = rest < batch_size ? rest : batch_size;
	process_batch(s, b, bcount);
    }
}

void simulate(state_t *s, int count, update_t update_mode, int dinterval, bool display) {
    s->update_mode = update_mode;
    int i;
    /* Compute and show initial state */
    bool show_counts = display;
    int batch_size;
    switch(update_mode) {
    case UPDATE_SYNCHRONOUS:
	batch_size = s->nrat;
	break;
    case UPDATE_RAT:
	batch_size = 1;
	break;
    case UPDATE_BATCH:
	batch_size = s->batch_size;
	break;
    default:
	outmsg("WARNING: Unknown update mode.  Using batch mode\n");
	batch_size = s->batch_size;
	break;
    }
    take_census(s);
    if (display) {
	show(s, show_counts);
    }
    for (i = 0; i < count; i++) {
	run_step(s, batch_size);
	if (display) {
	    show_counts = (((i+1) % dinterval) == 0) || (i == count-1);
	    show(s, show_counts);
	}
    }
    if (display && s->process_id == 0)
	done();
}

