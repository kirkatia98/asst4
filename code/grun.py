#!/usr/bin/python

# Implementation of GraphRat simulation
import sys
import getopt
import datetime
import math

import rutil
import gengraph
import sim
import viz

def usage(name):
    print "Usage: %s [-h] [-d] [-g GFILE] [-r RFILE] [-n STEPS] [-s SEED] [-u (s|r|b)] [-i INT] [-m (q|s|d)] [-p PERIOD] [-v (a|h|b)] [-c CFILE]"
    print "\t-h        Print this message"
    print "\t-d        Operate in driven mode, serving as visualizer for another simulator"
    print "\t          In driven mode, only additional options -m, -p, -v, and -c are useful"
    print "\t-g GFILE  Graph file"
    print "\t-r RFILE  Initial rat position file"
    print "\t-n STEPS  Number of simulation steps"
    print "\t-s SEED   Initial RNG seed"
    print "\t-u UPDT   Update mode:"
    print "\t          s: Synchronous.   Compute all new states and then update all."
    print "\t          r: Rat order:     Compute and update each rat state in sequence"
    print "\t          b: Batched.       Repeatedly compute states for small batches of rats and then update"
    print "\t-i INT    Generate image only once every INT steps"
    print "\t-m MODE   Output mode:"
    print "\t          q: Quiet.  Only statistics"
    print "\t          s: Step.   Show result of each step (Default)"
    print "\t          d: Drive.  Generate data to drive another program operating as visualizer"
    print "\t-p PERIOD Target refresh period (seconds)"
    print "\t-v VIS    Visualization Mode:"
    print "\t          b: Both    Show both ways (default)"
    print "\t          a: ASCII.  Print as numbers on grid"
    print "\t          h: Heatmap Show as graphical heatmap"
    print "\t-c CFILE  Capture final state as image (extensions .jpg and .png supported)"
    sys.exit(0)

# Enumerated type for output mode
# Extends limited form of basic simulator
class OutputMode:
    drive, quiet, step, error = range(4)

    def parse(self, name):
        if len(name) != 1:
            return self.error
        elif name == 'd':
            return self.drive
        elif name == 'q':
            return self.quiet
        elif name == 's':
            return self.step
        else:
            return self.error


# Generalizer of simulator to support multiple output modes
class VizSimulator(sim.Simulator):

    vizMode = viz.VizMode.heatmap
    formatter = None
    displayInterval = 1

    def __init__(self, graph, verb = OutputMode.step, vizMode = viz.VizMode.heatmap):
        sim.Simulator.__init__(self, graph)
        self.formatter = None
        self.verb = verb
        self.vizMode = vizMode

    def finish(self, fname = ""):
        if self.formatter is not None:
            self.formatter.finish(fname)

    # Display graph
    def show(self, period = 0.0, last = False):
        if self.formatter is None:
            k = int(math.sqrt(len(self.nodes)))
            self.formatter = viz.Formatter(k, self.ratCount(), viz = self.vizMode)
        else:
            self.formatter.reset()
        self.formatter.printLine("t = %d." % self.time)
        self.formatter.show(self.populationList(), period)


    def errorMsg(self, text):
        if self.formatter:
            self.formatter.printLine(text)
        else:
            if text[-1] != '\n':
                text += '\n'
            sys.stderr.write(text)


    # Return display to normal mode
    def finishDynamic(self, wait = False):
        if self.formatter is not None:
            self.formatter.finishDynamic(wait = wait)

    def simulate(self, stepCount = 1, update = sim.UpdateMode.synchronous, period = 0.0, displayInterval = 1):
        tstart = datetime.datetime.now()
        # Determine batch size
        rcount = len(self.rats)
        bsize = rcount
        if update == sim.UpdateMode.batch:
            bsize = self.batchSize
        elif update == sim.UpdateMode.ratOrder:
            bsize = 1
        if self.verb == OutputMode.step:
            self.show(period = period)
        elif self.verb == OutputMode.drive:
            self.driveOut()
        for step in xrange(stepCount):
            ridx = 0
            while ridx < rcount:
                bcount = min(bsize, rcount - ridx)
                for i in xrange(bcount):
                    r = self.rats[i+ridx]
                    r.next(loadFactor = self.loadFactor)
                for i in xrange(bcount):
                    r = self.rats[i+ridx]
                    r.move()
                ridx += bcount
            self.time += 1
            display = step == stepCount-1 or ((step+1) % displayInterval) == 0
            if display and self.verb == OutputMode.step:
                self.show(period = period)
            elif self.verb == OutputMode.drive:
                self.driveOut()
        self.finishDynamic()
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        rops = self.ratCount() * stepCount
        mrps = 1e-6 * float(rops)/secs
        if self.verb in [OutputMode.step]:
            self.show(period = 0.0, last = True)
        if self.verb in [OutputMode.quiet]:
            print "Elapsed time = %.2f seconds.  Ran at %.2f Mega Rats Per Second" % (secs, mrps)
        if self.verb == OutputMode.drive:
            self.driveDone()

# Special class to implement simulator in "driven mode"
# This mode enables another simulator to generate data
# and this simulator only to serve as a visualization tool
class DrivenSimulator(VizSimulator):
    # In this mode, don't model or track rats
    # Must keep track of rat count, since don't maintain list of rats
    nrats = 0

    def __init__(self, verb = OutputMode.quiet, vizMode = viz.VizMode.heatmap):
        self.nrats = 0
        self.nodes = []
        self.rats = []
        self.time = 0
        self.verb = verb
        self.formatter = None
        self.vizMode = vizMode

    def restart(self):
        self.time = 0
        
    def ratCount(self):
        return self.nrats

    def loadCounts(self):
        id = -1
        for line in sys.stdin:
            if line[-1] == '\n':
                line = line[:-1]
            tokens = line.split()
            if id == -1:
                if len(tokens) >= 1 and tokens[0] == "DONE":
                    return tokens[0]
                if len(tokens) < 1 or tokens[0] != "STEP":
                    self.errorMsg("Invalid driver input.  First line contents '%s'" % line)
                    return "ERROR"
                try:
                    ncount, self.nrats = map(int, tokens[1:])
                except Exception as e:
                    self.errorMsg("Failed to receive parameter line from driver: %s.  Line contents '%s'" % (e, line))
                    return "ERROR"
                if self.nodes == []:
                    self.nodes = [sim.Node(nid) for nid in xrange(ncount)]
                else:
                    for nid in xrange(ncount):
                        self.nodes[nid].reset()
            elif len(tokens) == 1 and tokens[0] == "END":
                break
            elif len(tokens) == 1:
                try:
                    count = int(tokens[0])
                except Exception as e:
                    self.errorMsg("Failed to receive input for node %d from driver: %s.  Line contents '%s'" % (id, e, line))
                    return "ERROR"
                self.nodes[id].ratCount = count
            else:
                self.errorMsg("Failed to receive input for node %d from driver: %s.  Line contents '%s'" % (id, e, line))
            id += 1
        return "EMPTY" if id <= 1 else "OK"
                
    def finishSim(self, tstart, count):
        self.finishDynamic()
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        rops = self.ratCount() * count
        mrps = 1e-6 * float(rops)/secs
        if self.verb in [OutputMode.step]:
            self.show(period = 0.0, last = True)
        if self.verb in [OutputMode.quiet]:
            print "Elapsed time = %.2f seconds.  Ran at %.2f Mega Rats Per Second" % (secs, mrps)

    def simulate(self, stepCount = 1, update = sim.UpdateMode.synchronous, period = 0.0, displayInterval = 1):
        tstart = datetime.datetime.now()
        realStepCount = 0
        code = self.loadCounts()
        if code == "DONE":
            self.finishSim(tstart, realStepCount)
            return
        elif code == "ERROR":
            self.errorMsg("Initial load failed.  Exiting")
            self.finishDynamic()
            return
        if self.verb == OutputMode.step:
            self.show(period = 0.0)
            # Force delay after showing initial state
            if period > 0:
                self.show(period = period)
        while True:
            code = self.loadCounts()
            if code == "DONE":
                self.finishSim(tstart, realStepCount)
                return
            elif code not in  ["OK", "EMPTY"]:
                self.finishDynamic()
                self.errorMsg("Driver failed.  Return code = '%s'.  Exiting" % code)
                return
            self.time += 1
            realStepCount += 1
            if code == "OK" and self.verb == OutputMode.step:
                self.show(period = period)
                

def run(name, args):
    gfname = ""
    irfname = ""
    steps = 1
    seed = rutil.DEFAULTSEED
    period = 0.1
    drivenMode = False
    vm = OutputMode()
    verb = vm.step
    displayInterval = 1
    updateMode = sim.UpdateMode.batch
    vizm = viz.VizMode()
    vizMode = vizm.heatmap
    captureFile = ""
    optlist, args = getopt.getopt(args, "hdg:r:R:n:s:u:m:p:i:v:c:")
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
            sys.exit(0)
        if opt == '-g':
            gfname = val
        if opt == '-r':
            irfname = val
        if opt == '-n':
            steps = int(val)
        if opt == '-s':
            seed = int(val)
        if opt == '-u':
            if len(val) != 1 or val not in "bsr":
                print "Error.  Unrecognized update mode '%s'" % val
                usage(name)
                return
            if val == 's':
                updateMode = sim.UpdateMode.synchronous
            elif val == 'r':
                updateMode = sim.UpdateMode.ratOrder
            else:
                updateMode = sim.UpdateMode.batch
        if opt == '-m':
            verb = vm.parse(val)
            if verb == vm.error:
                print "Error.  Unrecognized output mode '%s'" % val
                usage(name)
                return
        if opt == '-p':
            period = float(val)
        if opt == '-d':
            drivenMode = True
        if opt == '-i':
            displayInterval = int(val)
        if opt == '-v':
            vizMode = vizm.parse(val)
            if vizMode == vizm.error:
                print "Error.  Invalid visualization mode '%s'" % val
                usage(name)
                return
        if opt == '-c':
            captureFile = val
    if drivenMode:
        s = DrivenSimulator(verb = verb, vizMode = vizMode)
    else:
        if gfname == "":
            print "Error.  Need graph file"
            usage(name)
            return
        if irfname == "":
            print "Error.  Need file of initial rat positions"
            usage(name)
            return
        g = gengraph.Graph()
        if not g.load(gfname):
            return
        s = sim.Simulator(g) if verb == vm.drive else VizSimulator(g, verb = verb, vizMode = vizMode)
        if not s.loadRats(irfname, seed):
            return
    try:
        if verb == vm.drive:
            s.simulate(steps, update = updateMode, displayInterval = displayInterval)
        else:
            s.simulate(steps, update = updateMode, period = period, displayInterval = displayInterval)
    except Exception as E:
        s.errorMsg("Error: %s" % E)
        s.finish()
        return
    if verb != vm.drive:
        s.finish(captureFile)
    
if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
