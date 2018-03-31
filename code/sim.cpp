#include "crun.h"


//Fetch pre computed weight for that count
static inline double compute_weight(state_t *s, int nid) {
    int count = s->rat_count[nid];
    return s->pre_computed[count];
}

/* Compute sum of sumweights in region of nid */
static inline double compute_sum_weight(state_t *s, int nid) {
    graph_t *g = s->g;
    int eid_end = g->neighbor_start[nid+1];
    return g->gsums[eid_end - 1];
}

/** DEBUGGING CODE **/
#if DEBUG
/* Handy place to set a breakpoint */
static void panic(char *msg) {
    outmsg("PANIC: %s\n", msg);
}

static void show_weights(state_t *s) {
    int nid, eid;
    graph_t *g = s->g;
    int nnode = g->nnode;
    int *neighbor = g->neighbor;
    outmsg("Weights\n");
    for (nid = 0; nid < nnode; nid++) {
	int eid_start = g->neighbor_start[nid];
	int eid_end  = g->neighbor_start[nid+1];
	outmsg("In show_weights:\n");
	compute_sum_weight(s, nid);
	outmsg("%d: [sum = %.3f]", nid, compute_sum_weight(s, nid));
	for (eid = eid_start; eid < eid_end; eid++) {
	    outmsg(" %.3f", compute_weight(s, neighbor[eid]));
	}
	outmsg("\n");
    }
}
#endif


/* Recompute all node counts according to rat population */
void take_census(state_t *s) {
    graph_t *g = s->g;
    int nnode = g->nnode;



    //for each node, fill in its weight in the self edge index
    int nid, eid;
    for (nid = 0; nid < nnode; nid++)
    {
        eid = g->neighbor_start[nid];
        g->gsums[eid] = compute_weight(s, nid);
    }

    //for each node fill in the accumulation of the weights of its neighbors
    for (nid = 0; nid < nnode; nid++)
    {
        double sum = 0;
        for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++)
        {
            //find neighbor's weight in gsum
            int neighboredge = g->neighbor_start[g->neighbor[eid]];
            double neighborweights = g->gsums[neighboredge];

            sum += neighborweights;
            g->gsums[eid] = sum;
        }
    }
}


#define NEIGHBORS 16
/*
  Given list of integer counts, generate real-valued weights
  and use these to flip random coin returning value between 0 and len-1
*/
static inline int next_random_move(state_t *s, int r) {
    int nid = s->rat_position[r];
    int nnid = -1;
    random_t *seedp = &s->rat_seed[r];

    //bounds of search
    graph_t *g = s->g;
    int lo = g->neighbor_start[nid];
    int hi = g->neighbor_start[nid+1];

    double tsum = g->gsums[hi - 1];
    int eid;

    double val = next_random_float(seedp, tsum);

    //half linear search
    if(hi - lo <= NEIGHBORS)
    {
        //from end
        if (val > (tsum / 2.0f)) {
            eid = hi - 1;
            while (g->gsums[eid - 1] > val)
            {
                eid--;
                if (eid == lo) break;
            }
            return g->neighbor[eid];
        } //from beginning
        else {
            eid = lo;
            while (g->gsums[eid] <= val)
            {
                eid++;
            }
            return g->neighbor[eid];
        }
    }
    else
    { //binary search
        int mid;
        int beg = lo; //beginning of array
        while(lo < hi)
        {
            mid = lo + (hi - lo) / 2; // no overflow

            if (val < g->gsums[mid])
            {
                //either it's the first element or its first one strictly
                // bigger than val
                if (mid == beg || val >= g->gsums[mid - 1])
                    return g->neighbor[mid];

                else
                    hi = mid;
            }
            else
                lo = mid+1;
        }
    }
#if DEBUG
    if (nnid == -1) {
        /* Shouldn't get here */
        int degree = g->neighbor_start[nid+1] - g->neighbor_start[nid];
        outmsg("Internal error.  next_random_move.  Didn't find valid move.  Node %d. Degree = %d, Target = %.2f/%.2f",
               nid, degree, val, tsum);
        nnid = 0;
    }
#endif
    return nnid;
}

//TODO: do based on domain
static void process_batch(state_t *s, int bstart, int bcount) {
    graph_t *g = s->g;
    int rid, nid;

#if MPI
    array_t* m = s->mpi;
    int snode = m->ndisp[s->process_id];
    int enode = m->ndisp[s->process_id + 1];
#endif

    //for each rat compute the next positions (every processes)
    for (rid = bstart; rid < bstart+bcount; rid++)
    {
        s->next_position[rid] = next_random_move(s, rid);
    }

    for (rid = bstart; rid < bstart+bcount; rid++)
    {
        int onid = s->rat_position[rid];
        int nnid = s->next_position[rid];

#if MPI

        if(snode <= onid && onid < enode)
                s->delta[onid - snode]--;
        if(snode <= nnid && nnid < enode)
                s->delta[nnid - snode]++;
#else
        s->delta[onid]--;
        s->delta[nnid]++;
#endif
        s->rat_position[rid] = nnid;
    }

    for(nid = 0; nid < s->my_nodes; nid++)
    {
#if MPI
        if(g->tiles_per_side > 1)
            s->local_rat_count[nid]+= s->delta[nid];
        else
            s->rat_count[nid]+= s->delta[nid];
#else
        s->rat_count[nid]+= s->delta[nid];
#endif
        s->delta[nid] = 0;
    }
}


//TODO: do based on domain
static void process_rats(state_t *s, int bstart, int bcount) {
    graph_t *g = s->g;
    int rid, nid;

#if MPI
    array_t* m = s->mpi;
    int snode = m->ndisp[s->process_id];
    int enode = m->ndisp[s->process_id + 1];

    int srat = m->rdisp[s->process_id];
    int erat = m->rdisp[s->process_id + 1];
#else
    int srat = bstart;
    int erat = bstart + bcount;
#endif


    //for each rat compute the next positions
    for (rid = srat; rid < erat; rid++)
    {
        s->next_position[rid] = next_random_move(s, rid); 
    }
#if MPI

    //all gather for the batch of rats
    if(g->titles_per_side > 1)
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   s->next_position, m->rsend, m->rdisp, MPI_INT,
                   MPI_COMM_WORLD);
#endif
    for (rid = bstart; rid < bstart + bcount; rid++)
    {
        int onid = s->rat_position[rid];
        int nnid = s->next_position[rid];

#if MPI

        if(g->tiles_per_side > 1)
            if(snode <= onid && onid < enode)
                s->local_rat_count[onid - snode]--;
            if(snode <= nnid && nnid < enode)
                s->local_rat_count[nnid - snode]++;
        else
            s->rat_count[nid]+= s->delta[nid];
#else
        s->rat_count[onid]--;
        s->rat_count[nnid]++;
#endif
        s->rat_position[rid] = nnid;
    }
}

static void run_step(state_t *s, int batch_size) {
    int b, bcount;
    graph_t* g = s->g;
    array_t* m = s->mpi;

    for (b = 0; b < s->nrat; b += batch_size) {
        int rest = s->nrat - b;
        bcount = rest < batch_size ? rest : batch_size;

        //update gsums and distribute
        if(s->process_id == 0)
            take_census(s);
#if MPI
        if(g->tiles_per_side > 1)
            MPI_Bcast(g->gsums, g->nnode + g->nedge, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
        if(batch_size == s->nrat)
            process_rats(s, b, bcount);
        else
            process_batch(s, b, bcount);

#if MPI
        //gather all the rat counts
        if(g->tiles_per_side > 1)
            MPI_Gatherv(s->local_rat_count, s->my_nodes, MPI_INT,
                    s->rat_count, m->nsend, m->ndisp, MPI_INT,
                    0, MPI_COMM_WORLD);
#endif

    }
}

void simulate(state_t *s, int count, update_t update_mode, int dinterval, bool display) {
    bool mpi_master = (s->process_id == 0);

    s->update_mode = update_mode;
    graph_t* g = s->g;
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


    if (display && mpi_master) {
	    show(s, show_counts);
    }

    for (i = 0; i < count; i++) {

        run_step(s, batch_size);

        if (display && mpi_master) {
            show_counts = (((i+1) % dinterval) == 0) || (i == count-1);
            show(s, show_counts);
        }
#if MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    if (display && mpi_master)
	    done();
}

