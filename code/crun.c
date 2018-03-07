/* C implementation of graphrats simulator */

#include <string.h>
#include <getopt.h>

#include "crun.h"

static void usage(char *name) {
    char *use_string = "-g GFILE -r RFILE [-n STEPS] [-s SEED] [-u (r|b|s)] [-q] [-i INT]";
    outmsg("Usage: %s %s\n", name, use_string);
    outmsg("   -h        Print this message\n");
    outmsg("   -g GFILE  Graph file\n");
    outmsg("   -r RFILE  Initial rat position file\n");
    outmsg("   -n STEPS  Number of simulation steps\n");
    outmsg("   -s SEED   Initial RNG seed\n");
    outmsg("   -u UPDT   Update mode:\n");
    outmsg("             s: Synchronous.  Compute all new states and then update all\n");
    outmsg("             r: Rat order.    Compute update each rat state in sequence\n");
    outmsg("             b: Batched.      Repeatedly compute states for small batches of rats and then update\n");
    outmsg("   -q        Operate in quiet mode.  Do not generate simulation results\n");
    outmsg("   -i INT    Display update interval\n");
    done();
    exit(0);
}

int main(int argc, char *argv[]) {
    FILE *gfile = NULL;
    FILE *rfile = NULL;
    int steps = 1;
    int dinterval = 1;
    random_t global_seed = DEFAULTSEED;
    update_t update_mode = UPDATE_BATCH;
    int process_count = 1;
    int process_id = 0;
    int c;
    graph_t *g = NULL;
    state_t *s = NULL;
    bool display = true;
#if MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
#endif
    bool mpi_master = process_id == 0;
    char *optstring = "hg:r:R:n:s:u:i:q";
    while ((c = getopt(argc, argv, optstring)) != -1) {
	switch(c) {
	case 'h':
	    if (!mpi_master) break;
	    usage(argv[0]);
	    break;
	case 'g':
	    if (!mpi_master) break;
	    gfile = fopen(optarg, "r");
	    if (gfile == NULL) {
		outmsg("Couldn't open graph file %s\n", optarg);
		done();
		exit(1);
	    }
	    break;
	case 'r':
	    if (!mpi_master) break;
	    rfile = fopen(optarg, "r");
	    if (rfile == NULL) {
		outmsg("Couldn't open rat position file %s\n", optarg);
		done();
		exit(1);
	    }
	    break;
	case 'n':
	    steps = atoi(optarg);
	    break;
	case 's':
	    global_seed = strtoul(optarg, NULL, 0);
	    break;
	case 'u':
	    if (optarg[0] == 'r')
		update_mode = UPDATE_RAT;
	    else if (optarg[0] == 'b')
		update_mode = UPDATE_BATCH;
	    else if (optarg[0] == 's')
		update_mode = UPDATE_SYNCHRONOUS;
	    else {
		if (!mpi_master) exit(1);
		outmsg("Invalid update mode '%c'\n", optarg[0]);
		usage(argv[0]);
		done();
		exit(1);
	    }
	    break;
	case 'q':
	    display = false;
	    break;
	case 'i':
	    dinterval = atoi(optarg);
	    break;
	default:
	    if (!mpi_master) break;
	    outmsg("Unknown option '%c'\n", c);
	    usage(argv[0]);
	    done();
	    exit(1);
	}
    }
    if (mpi_master) {
      	if (gfile == NULL) {
	    outmsg("Need graph file\n");
	    usage(argv[0]);
	}
	if (rfile == NULL) {
	    outmsg("Need initial rat position file\n");
	    usage(argv[0]);
	}


	g = read_graph(gfile);
	if (g == NULL) {
	    done();
	    exit(1);
	}
	s = read_rats(g, rfile, global_seed);
	if (s == NULL) {
	    done();
	    exit(1);
	}
	s->nprocess = process_count;
	s->process_id = process_id;
	/* The master should distribute the graph & the rats to the other nodes */
    } else {
	/* The other nodes should get the graph & the rats from the master */
    }

    double start = currentSeconds();
    if (mpi_master)
	// Run sequential simulator on process #0
	simulate(s, steps, update_mode, dinterval, display);
    double delta = currentSeconds() - start;
    if (mpi_master) {
	outmsg("%d steps, %d rats, %.3f seconds\n", steps, s->nrat, delta);
    }
#if MPI
    MPI_Finalize();
#endif    
    return 0;
}

