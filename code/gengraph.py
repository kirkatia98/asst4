#!/usr/bin/python

# Code for generating and reading graphrat mazes
# Parameters:
# k: Grid size will be k * k nodes
# t T: Tile size (superimpose TxT tiles, each containing one node connected to all other nodes in tile)
# fractal: Optionally superimpose fractal hierarchy


import getopt
import sys
import math
import string
import datetime

import rutil

def usage(name):
    print "Usage: %s [-h] [-k K] [-f] [-t T] [-o OUT] [-s SEED]"
    print "\t-h     Print this message"
    print "\t-k K   Base graph as k x k grid"
    print "\t-f     Create fractal graph"
    print "\t-t T   Add tiles, each with TxT nodes"
    print "\t-o OUT Specify output file"

class RatMode:
    # Different options for specifying initial rat state
    (uniform, diagonal, upleft, lowright) =  range(4)
    modeNames = ["uniform", "diagonal", "upper-left", "lower-right"]

class Graph:
    k = 0
    nodeCount = 0
    tileMax = 1
    edges = {}  # Maps edges to True.  Include both directions
    commentList = []  # Documentation about how generated

    def __init__(self, k = 0, fractal = False, tile = 0):
        self.generate(k, fractal, tile)

    def generate(self, k = 10, fractal = False, tile = 0):
        self.commentList = []
        tgen = datetime.datetime.now()
        self.commentList.append("# Generated %s" % tgen.ctime())
        self.commentList.append("# Parameters: k = %d, %s" % (k, "fractal" if fractal else "tiled" if tile > 1 else "uniform"))
        self.k = k
        self.nodeCount = k * k
        nodes = [i for i in range(self.nodeCount)]
        self.edges = {}
        # Generate grid edges
        for r in range(k):
            for c in range(k):
                own = self.id(r, c)
                north = self.id(r-1, c)
                if north >= 0:
                    self.addEdge(own, north)
                south = self.id(r+1, c)
                if south >= 0:
                    self.addEdge(own, south)
                west = self.id(r, c-1)
                if west >= 0:
                    self.addEdge(own, west)
                east = self.id(r, c+1)
                if east >= 0:
                    self.addEdge(own, east)
        if fractal:
            # Generate fractal graph
            (x, y, w) = (0, 0, k)
            self.fracture(x, y, w)
            self.tileMax = self.k/2
        if tile > 0:
            self.tile(tile)
            self.tileMax = tile

    def fracture(self, x, y, w):
        if (w % 2) != 0:
            self.makeHubs(x, y, w, w)
            return
        nw = w / 2
        self.makeHubs(x, y, w, nw, 4, 1)
        self.makeHubs(x, y+nw, nw, nw, 1, 1)
        self.fracture(x+nw, y+nw, nw)

    def tile(self, tile):
        for x in range(0, self.k, tile):
            w = min(tile, self.k - x)
            for y in range(0, self.k-tile+1, tile):
                h = min(tile, self.k - y)
                self.makeHubs(x, y, w, h)

    def makeHubs(self, x, y, w, h, xcount=1, ycount=1):
        if w <= 2*xcount:
            wsep = w/xcount
        else:
            wsep = w/(xcount + 1)

        hsep = h/(ycount + 1)
        
        if w <= xcount:
            cxList = [x + wsep * i for i in range(xcount)]
        elif w <= 2*xcount:
            cxList = [1 + x + wsep * i for i in range(xcount)]
        else:
            cxList = [x + wsep * (i + 1) for i in range(xcount)]
        cyList = [y + hsep * (i + 1) for i in range(ycount)]
        for cx in cxList:
            for cy in cyList:
                cid = self.id(cy, cx)
                for j in range(w):
                    for i in range(h):
                        id = self.id(y+i, x+j)
                        self.addEdge(cid, id)
        
    # Check whether string is a comment
    def isComment(self, s):
        # Strip off leading whitespace
        while len(s) > 0 and s[0] in string.whitespace:
            s = s[1:]
        return len(s) == 0 or s[0] == '#'

    # Load graph from file
    def load(self, fname = ""):
        self.k = 0
        self.nodeCount = 0
        self.edges = {}
        self.tileMax = 1
        if fname == "":
            f = sys.stdin
        else:
            try:
                f = open(fname, "r")
            except:
                sys.stderr.write("Could not open file '%s'\n" % fname)
                return False
        expectedEgeCount = 0
        realEdgeCount = 0
        for line in f:
            if self.isComment(line):
                continue
            args = line.split()
            if self.k == 0:
                self.nodeCount = int(args[0])
                self.k = int(math.sqrt(self.nodeCount))
                expectedEdgeCount = int(args[1])
                if len(args) >= 3:
                    self.tileMax = int(args[2])
                else:
                    self.tileMax = self.k
            else:
                i = int(args[0])
                j = int(args[1])
                if self.addEdge(i,j):
                    # Since addEdge puts both (i,j) and (j,i) into set, only half of the
                    # edges will return True from addEdge
                    realEdgeCount += 2 
        if fname != "":
            f.close()
        if realEdgeCount != expectedEdgeCount:
            sys.stderr.write("Error reading graph file '%s'.  Expected %d edges.  Found %d\n" % (fname, expectedEdgeCount, realEdgeCount))
            return False
        else:
            sys.stderr.write("Read graph with %d nodes and %d edges\n" % (self.nodeCount, realEdgeCount))
            return True
 
    def id(self, r, c):
        if r < 0 or r >= self.k:
            return -1
        if c < 0 or c >= self.k:
            return -1
        return r * self.k + c

    def addEdge(self, i, j):
        if i < 0 or i >= self.nodeCount:
            sys.stderr.write("Error: Invalid from node id %d\n" % i)
        if j < 0 or j >= self.nodeCount:
            sys.stderr.write("Error: Invalid to node id %d\n" % j)
        if i != j and (i,j) not in self.edges:
            self.edges[(i,j)] = True
            self.edges[(j,i)] = True
            return True
        return False
            
    def edgeList(self):
        elist = [e for e in self.edges]
        elist.sort()
        return elist

    # Generate list with entry with each node, giving its degree (including self)
    def degreeList(self):
        result = [1] * self.nodeCount
        for e in self.edges:
            idx = e[0]
            result[idx] += 1
        return result

    # Store graph
    def store(self, fname = ""):
        if fname == "":
            f = sys.stdout
        else:
            try:
                f = open(fname, "w")
            except:
                sys.stderr.write("Error.  Couldn't open file '%s' for writing\n" % (fname))
                return False
        elist = self.edgeList()
        f.write("%d %d %d\n" % (self.nodeCount, len(self.edges), self.tileMax))
        for c in self.commentList:
            f.write(c + '\n')
        for e in elist:
            f.write("%d %d\n" % e)
        if fname != "":
            f.close()
        return True

    # Generate rats for graph and write to file
    def makeRats(self, fname = "", mode = RatMode.uniform, load = 1, seed = rutil.DEFAULTSEED):
        clist = []
        tgen = datetime.datetime.now()
        clist.append("# Generated %s" % tgen.ctime())
        clist.append("# Parameters: load = %d, mode = %s, seed = %d" % (load, RatMode.modeNames[mode], seed))
        rng = rutil.RNG([seed])
        if fname == "":
            f = sys.stdout
        else:
            try:
                f = open(fname, "w")
            except:
                "Couldn't open output file '%s'"
                return False
        rlist = []
        if mode == RatMode.uniform:
            rlist = range(self.nodeCount)
        elif mode == RatMode.diagonal:
            rlist = [(self.k+1)*i for i in range(self.k)]
        elif mode == RatMode.upleft:
            rlist = [0]
        elif mode == RatMode.lowright:
            rlist = [self.nodeCount-1]
        else:
            sys.stderr.write("ERROR: Invalid rat mode\n")
            return False
        factor = self.nodeCount * load / len(rlist)
        fullRlist = rlist * factor
        if len(rlist) > 0:
            fullRlist = rng.permute(fullRlist)
        # Print it out
        f.write("%d %d\n" % (self.nodeCount, len(fullRlist)))
        for c in clist:
            f.write(c + '\n')
        for id in fullRlist:
            f.write("%d\n" % id)
        if fname != "":
            f.close()
        return True

# Runtime code for graph generation
def run(name, args):
    k = 10
    fractal = False
    fname = ""
    optlist, args = getopt.getopt(args, "hk:fo:")
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
            sys.exit(0)
        if opt == '-k':
            k = int(val)
        if opt == '-f':
            fractal = True
        if opt == '-o':
            fname = val
    g = Graph(k = k, fractal = fractal)
    g.store(fname = fname)

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])


        
    
