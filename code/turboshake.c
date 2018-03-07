/* Code to stress the processor, getting it to find a stable level for its clock */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>


#include "cycletimer.h"

#define ARRAY_SIZE (30 * 1000)
#define DEFAULT_TURBO_SECONDS 2
#define DEFAULT_SLEEP_SECONDS 5


void usage(char *pname) {
    printf("Usage: %s [-h] [-v] [-s sleepSeconds] [-t turboSeconds] \n", pname);
    exit(1);
}

volatile double total_sum = 0;

int verbose = 0;

void turbo(int seconds) {
    double tstart = currentSeconds();
    
    double *data = calloc(ARRAY_SIZE, sizeof(double));

    int iters = 0;

    while (currentSeconds() - tstart < seconds) {
	int i;
	double sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+:sum)
	for (i = 0; i < ARRAY_SIZE; i++)
	{
	    int j = 0;
	    for (j = 0; j < ARRAY_SIZE; j++) {
		data[i] += (double) (i + j);
	    }
	    sum += data[i];
	}
	total_sum += sum;
	iters++;
    }
    if (verbose)
	printf("Seconds = %.2f.  Iterations = %d.  Total sums = %ld\n",
	       currentSeconds() - tstart, iters, (long) iters * ARRAY_SIZE * (ARRAY_SIZE + 1));
}


int main(int argc, char *argv[]) {
    int turboSeconds = DEFAULT_TURBO_SECONDS;
    int sleepSeconds = DEFAULT_SLEEP_SECONDS;
    verbose = 0;
    int c;
    while ((c = getopt(argc, argv, "hvs:t:")) != -1) {
	switch (c) {
	case 'h':
	    usage(argv[0]);
	    break;
	case 'v':
	    verbose = 1;
	    break;
	case 's':
	    sleepSeconds = atoi(optarg);
	    break;
	case 't':
	    turboSeconds = atoi(optarg);
	    break;
	default:
	    printf("Unkown option '%c'\n", c);
	    usage(argv[0]);
	    break;
	}
    }
    if (sleepSeconds > 0)
	sleep(sleepSeconds);
    turbo(turboSeconds);
    return 0;
}
