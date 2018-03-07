#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "crun.h"

graph_t *new_graph(int nnode, int nedge, int tile_max) {
    bool ok = true;
    graph_t *g = malloc(sizeof(graph_t));
    if (g == NULL)
	return NULL;
    g->nnode = nnode;
    g->nrow = (int) sqrt(g->nnode);
    g->tile_max = tile_max > 0 ? tile_max : g->nrow;
    g->nedge = nedge;
    g->neighbor = calloc(nnode + nedge, sizeof(int));
    ok = ok && g->neighbor != NULL;
    g->neighbor_start = calloc(nnode + 1, sizeof(int));
    ok = ok && g->neighbor_start != NULL;
    if (!ok) {
	outmsg("Couldn't allocate graph data structures");
	return NULL;
    }
    return g;
}

void free_graph(graph_t *g) {
    free(g->neighbor);
    free(g->neighbor_start);
    free(g);
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

/* Read in graph file and build graph data structure */
graph_t *read_graph(FILE *infile) {
    char linebuf[MAXLINE];
    int nnode, nedge;
    int tile_max = 0;
    int i, hid, tid;
    int nid, eid;

    // Read header information
    while (fgets(linebuf, MAXLINE, infile) != NULL) {
	if (!is_comment(linebuf))
	    break;
    }
    if (sscanf(linebuf, "%d %d %d", &nnode, &nedge, &tile_max) < 2) {
	outmsg("ERROR. Malformed graph file header (line 1)\n");
	return NULL;
    }
    graph_t *g = new_graph(nnode, nedge, tile_max);
    if (g == NULL)
	return g;

    nid = -1;
    // We're going to add self edges, so eid will keep track of all edges.
    eid = 0;  
    for (i = 0; i < nedge; i++) {
	while (fgets(linebuf, MAXLINE, infile) != NULL) {
	    if (!is_comment(linebuf))
		break;
	}
	if (sscanf(linebuf, "%d %d", &hid, &tid) != 2) {
	    outmsg("Line #%u of graph file malformed\n", i+2);
	    return false;
	}
	if (hid < 0 || hid >= nnode) {
	    outmsg("Invalid head index %d on line %d\n", hid, i+2);
	    return false;
	}
	if (tid < 0 || tid >= nnode) {
	    outmsg("Invalid tail index %d on line %d\n", tid, i+2);
	    return false;
	}
	if (hid < nid) {
	    outmsg("Head index %d on line %d out of order\n", hid, i+2);
	    return false;
	    
	}
	// Starting edges for new node(s)
	while (nid < hid) {
	    nid++;
	    g->neighbor_start[nid] = eid;
	    // Self edge
	    g->neighbor[eid++] = nid;
	}
	g->neighbor[eid++] = tid;
    }
    while (nid < nnode-1) {
	// Fill out any isolated nodes
	nid++;
	g->neighbor[eid++] = nid;
    }
    g->neighbor_start[nnode] = eid;
    outmsg("Loaded graph with %d nodes and %d edges\n", nnode, nedge);
#if DEBUG
    show_graph(g);
#endif
    return g;
}

#if DEBUG
void show_graph(graph_t *g) {
    int nid, eid;
    outmsg("Graph\n");
    for (nid = 0; nid < g->nnode; nid++) {
	outmsg("%d:", nid);
	for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++) {
	    outmsg(" %d", g->neighbor[eid]);
	}
	outmsg("\n");
    }
    
}
#endif
