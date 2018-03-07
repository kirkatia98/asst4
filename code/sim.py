#!/usr/bin/python

# Implementation of GraphRat simulation
# This file contains the core simulation modeling code

import sys
import datetime
import math
import string

import rutil
import gengraph


# Enumerated type for update mode:
# synchronous:  First compute all next states for all rats, and then move them
# ratOrder:     For each rat: compute its next state and move it immediately
class UpdateMode:
    ratOrder, batch, synchronous = range(3)


# Representation of single rat
class Rat:
    id = 0
    rng = None     # Random number generator
    node = None    # Graph node where rat is located
    newNode = None # Next node where rat will move

    def __init__(self, id, node, seed=rutil.DEFAULTSEED):
        self.id = id
        self.rng = rutil.RNG()
        self.node = node
        self.newNode = None
        node.addRat(self)
        self.reset(seed)

    # Reset RNG
    def reset(self, seed = rutil.DEFAULTSEED):
        self.rng.reseed([seed, self.id])

    # Next state computation
    def next(self, loadFactor = 1.0):
        loads = [float(nid.ratCount)/loadFactor for nid in self.node.region]
        idx = rutil.chooseMove(self.rng, loads)
        self.newNode = self.node.region[idx]

    # Update state
    def move(self):
        self.node.removeRat(self)
        self.node = self.newNode
        self.node.addRat(self)
        self.newNode = None

# Representation of single node
class Node:
    id = 0
    region = []  # self + neighbors
    ratCount = 0

    def __init__(self, id):
        self.id = id
        # Region is own node + adjacency list
        self.region = [self]
        self.reset()

    # Clear all rats
    def reset(self):
        self.ratCount = 0

    # Add node to adjacency list
    def addNeighbor(self, nd):
        self.region.append(nd)

    # Move rat to node
    def addRat(self, r):
        self.ratCount += 1

    # Move rat away from node
    def removeRat(self, r):
        self.ratCount -= 1


# Overall simulation.  This one only operates in "drive" or "benchmark" mode
class Simulator:
    nodes = []
    rats = []
    time = 0          # Number of steps simulated
    loadFactor = 0.0  # Ratio of rats to nodes
    batchSize = 0

    def __init__(self, graph):
        self.nodes = [Node(id) for id in xrange(graph.nodeCount)]
        for (hidx,tidx) in graph.edgeList():
            head = self.nodes[hidx]
            tail = self.nodes[tidx]
            head.addNeighbor(tail)
        self.time = 0

    # Check whether string is a comment
    def isComment(self, s):
        # Strip off leading whitespace
        while len(s) > 0 and s[0] in string.whitespace:
            s = s[1:]
        return len(s) == 0 or s[0] == '#'

    # Read rat position file
    def loadRats(self, fname = "", seed = rutil.DEFAULTSEED):
        ratPositions = []
        if fname == "":
            f = sys.stdin
        else:
            try:
                f = open(fname, "r")
            except:
                self.errorMsg("Couldn't open file '%s'" % fname)
                self.finish()
                return False
        first = True
        for line in f:
            if self.isComment(line):
                continue
            if first:
                ncount, rcount = map(int, line.split())
                if ncount != len(self.nodes):
                    self.errorMsg("Mismatch.  Graph has %d nodes.  Rat file has %d nodes.  No rats addded." % (len(self.nodes), ncount))
                    return
                first = False
            else:
               nid = int(line.split()[0])
               ratPositions.append(nid)
        f.close()
        self.restart(ratPositions, seed)
        sys.stderr.write("Loaded %d rats\n" % rcount)
        self.loadFactor = float(rcount) / len(self.nodes)
        self.batchSize = max(int(math.sqrt(rcount)), int(0.02 * rcount))
        return True

    # Restart simulation.  Use rat position array read from file
    def restart(self, ratPositions = [], seed = rutil.DEFAULTSEED):
        self.rats = []
        for n in self.nodes:
            n.reset()
        self.time = 0
        for rid in xrange(len(ratPositions)):
            nid = ratPositions[rid]
            if nid < 0 or nid >= len(self.nodes):
                self.errorMsg("Invalid rat position: %d.  Ignoring" % nid)
                return
            node = self.nodes[nid]
            rat = Rat(rid, node, seed)
            self.rats.append(rat)

    def ratCount(self):
        return len(self.rats)

    def finish(self):
        self.driveDone()

    # Write rat position file based on current state
    def storeRats(self, fname = ""):
        if fname == "":
            f = sys.stdout
        else:
            try:
                f = open(fname, "w")
            except:
                self.errorMsg("Couldn't open file '%s'" % fname)
                self.finish()
                return False
        f.write("%d %d\n" % (len(self.nodes), self.ratCount()))
        for r in self.rats:
            f.write("%d\n" % r.position.id)
        f.close()
        return True

    # Return list with count of rats for each node
    def populationList(self):
        return [nd.ratCount for nd in self.nodes]

    # Generate output suitable for reading into another copy of program running in driven or benchmark mode
    # First line is header "STEP" as identifier
    # Second line of form "N R", where N is number of nodes, and R is number of rats
    # Each successive line then lists the number of rats at each node
    # Terminate with line "END"
    def driveOut(self, f = sys.stdout, display = True):
        f.write("STEP %d %d\n" % (len(self.nodes), self.ratCount()))
        if display:
            for nd in self.nodes:
                f.write("%d\n" % nd.ratCount)
        f.write("END\n")
                
    # Final line of driver output, to indicate simulation has completed
    # It's a good idea to put this at the end of any output to signal the visualizer
    # that the program is terminating
    def driveDone(self, f = sys.stdout):
        f.write("DONE\n")

    # Should print any information messages on stderr, since stdout is being piped into another program
    def errorMsg(self, text):
        if text[-1] != '\n':
            text += '\n'
        sys.stderr.write(text)

    # Basic simulation step
    def simulate(self, stepCount = 1, update = UpdateMode.synchronous, displayInterval = 1):
        # Determine batch size
        bsize = len(self.rats)
        if update == UpdateMode.batch:
            bsize = self.batchSize
        elif update == UpdateMode.ratOrder:
            bsize = 1
        display = True
        # Emit initial state
        self.driveOut(display = display)
        for step in xrange(stepCount):
            ridx = 0
            while ridx < len(self.rats):
                bcount = min(bsize, len(self.rats) - ridx)
                for i in xrange(bcount):
                    r = self.rats[i+ridx]
                    r.next(loadFactor = self.loadFactor)
                for i in xrange(bcount):
                    r = self.rats[i+ridx]
                    r.move()
                ridx += bcount
            self.time += 1
            # Emit new state
            display = step == stepCount-1 or ((step+1) % displayInterval) == 0
            self.driveOut(display = display)
        self.driveDone()
                
                
                
                
