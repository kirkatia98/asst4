#!/usr/bin/python

import subprocess
import sys
import os
import getopt
import math
import datetime
import time

import grade

def usage(fname):
    ustring = "Usage: %s [-h] [-s SCALE] [-u UPDATELIST] [-f OUTFILE] [-S SSECS] [-T TSECS]" % fname
    ustring += " [-p PROCESSLIMIT]"
    print ustring 
    print "(All lists given as colon-separated text.)"
    print "    -h              Print this message"
    print "    -s SCALE        Reduce number of steps in each benchmark by specified factor"
    print "    -u UPDATELIST   Specify update modes(s):"
    print "       r: rat order"
    print "       s: synchronous"
    print "       b: batch"
    print "    -f OUTFILE      Create output file recording measurements"
    print "    -S SSECS        Let processor cool down by sleeping for SSECS seconds"
    print "    -T TSECS        Stimulate Turboboost by running job for TSECS seconds"
    print "    -p PROCESSLIMIT Specify upper limit on number of MPI processes"
    print "       If > 1, will run crun-mpi.  Else will run crun"
    sys.exit(0)

# Enumerated type for update mode:
# synchronous:  First compute all next states for all rats, and then move them
# ratOrder:     For each rat: compute its next state and move it immediately
class UpdateMode:
    ratOrder, batch, synchronous = range(3)
    flags = ['r', 'b', 's']

# General information
simProg = "./crun"
mpiSimProg = "./crun-mpi"

dataDir = "./data/"
outFile = None

# Program to stimulate Turboboost
turboProg = "./turboshake"
# How long to sleep before Turboboost
sleepSeconds = 10
# How long to run program to kick into Turboboost
turboSeconds = 2

# Dictionary of geometric means, indexed by (mode, threads)
gmeanDict = {}

def outmsg(s, noreturn = False):
    if len(s) > 0 and s[-1] != '\n' and not noreturn:
        s += "\n"
    sys.stdout.write(s)
    sys.stdout.flush()
    if outFile is not None:
        outFile.write(s)


runFlags = ["-q"]

# For computing geometric mean
logSum = 0.0
bcount = 0

# For keeping track of speedups
# Mapping from runtime parameters to MRPS
resultCache = {}

# Marker that allows filtering via grep
marker = "+++\t"
nomarker = "\t"

def reset():
    global logSum, bcount
    logSum = 0.0
    bcount = 0

# Graph/rat combinations: (graphSize, graphType, ratType, loadFactor)
benchmarkList = [
    (32400, 'u', 'u', 32),
    (32400, 't', 'u', 32),
    (32400, 'u', 'd', 32),
    (32400, 't', 'd', 32),
    ]


# Tests
synchRunList = [(1, 1000), (12, 1000)]
otherRunList = [(1, 500), (12, 500)]

def cmd(graphSize, graphType, ratType, loadFactor, stepCount, updateType, processCount, otherArgs = []):
    global bcount, logSum
    global cacheKey
    updateFlag = UpdateMode.flags[updateType]
    params = ["%5d" % graphSize, graphType, "%4d" % loadFactor, ratType, str(stepCount), updateFlag]
    cacheKey = ":".join(params)
    results = params + [str(processCount)]
    sizeName = str(graphSize)
    graphFileName = dataDir + "g-" + graphType + sizeName + ".gph"
    ratFileName = dataDir + "r-" + sizeName + '-' + ratType + str(loadFactor) + ".rats"
    clist = runFlags + ["-g", graphFileName, "-r", ratFileName, "-u", updateFlag, "-n", str(stepCount), "-i", str(stepCount)] + otherArgs
    if processCount > 1:
        gcmd = ["mpirun", "-np", str(processCount), mpiSimProg] + clist
    else:
        gcmd = [simProg] + clist
    gcmdLine = " ".join(gcmd)
    retcode = 1
    tstart = datetime.datetime.now()
    try:
        # File number of standard output
        stdoutFileNumber = 1
        simProcess = subprocess.Popen(gcmd, stderr = stdoutFileNumber)
        simProcess.wait()
        retcode = simProcess.returncode
    except Exception as e:
        print "Execution of command '%s' failed. %s" % (gcmdLine, e)
        return False
    if retcode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        rops = int(graphSize * loadFactor) * stepCount
        ssecs = "%.2f" % secs 
        results.append(ssecs)
        mrps = 1e-6 * float(rops)/secs
        if mrps > 0:
            logSum += math.log(mrps)
            bcount += 1
        smrps = "%7.2f" % mrps
        results.append(smrps)
        if cacheKey in resultCache:
            speedup = mrps / resultCache[cacheKey] 
            sspeedup = "(%5.2fX)" % speedup
            results.append(sspeedup)
        if processCount == 1:
            resultCache[cacheKey] = mrps
        pstring = marker + "\t".join(results)
        outmsg(pstring)
    else:
        print "Execution of command '%s' gave return code %d" % (gcmdLine, retcode)
        return False
    return True

def turbo():
    if turboSeconds > 0:
        if not os.path.exists(turboProg):
            print "Warning.  Cannot find program %s" % turboProg
            return
        print "Running %s to sleep for %d seconds and then go active for %d seconds" % (turboProg, sleepSeconds, turboSeconds)
        tcmd = [turboProg, '-s', str(sleepSeconds), '-t', str(turboSeconds)]
        simProcess = subprocess.Popen(tcmd)
        simProcess.wait()

def sweep(updateType, processLimit, scale, otherArgs):
    runList = synchRunList if updateType == UpdateMode.synchronous else otherRunList
    for rparams in runList:
        reset()
        (processCount, stepCount) = rparams
        if processCount > 1 and processLimit == 1:
            continue
        if processCount > processLimit:
            processCount = processLimit
        outmsg("\tNodes\tgtype\tlf\trtype\tsteps\tupdate\tprocs\tsecs\tMRPS")
        stepCount = stepCount / scale
        outmsg(nomarker + "---------" * 8)
        if processCount == 1:
            turbo()
        for bparams in benchmarkList:
            if processCount > 1:
                turbo()
            (graphSize, graphType, ratType, loadFactor) = bparams
            cmd(graphSize, graphType, ratType, loadFactor, stepCount, updateType, processCount, otherArgs)
        if bcount > 0:
            gmean = math.exp(logSum/bcount)
            updateFlag = UpdateMode.flags[updateType]
            gmeanDict[(updateFlag, processCount)] = gmean
            outmsg(marker + "Gmean\t\t\t\t\t%s\t%d\t\t%7.2f" % (updateFlag, processCount, gmean))
            outmsg(marker + "---------" * 8)

    
def run(name, args):
    global outFile, sleepSeconds, turboSeconds
    scale = 1
    updateList = [UpdateMode.batch, UpdateMode.synchronous]
    optString = "hms:u:p:f:S:T:"
    processLimit = 100
    optlist, args = getopt.getopt(args, optString)
    otherArgs = []
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
        elif opt == '-s':
            scale = float(val)
        elif opt == '-f':
            try:
                outFile = open(val, "w")
            except Exception as e:
                outFile = None
                outmsg("Couldn't open file '%s'" % val)
        elif opt == '-u':
            ulist = val.split(":")
            updateList = []
            for c in ulist:
                if c == 's':
                    updateList.append(UpdateMode.synchronous)
                elif c == 'b':
                    updateList.append(UpdateMode.batch)
                elif c == 'r':
                    updateList.append(UpdateMode.ratOrder)
                else:
                    print "Invalid update mode '%s'" % c
                    usage(name)
        elif opt == '-p':
            processLimit = int(val)
        elif opt == '-S':
            sleepSeconds = float(val)
        elif opt == '-T':
            turboSeconds = float(val)
        else:
            print "Uknown option '%s'" % opt
            usage(name)

    tstart = datetime.datetime.now()

    for u in updateList:
        sweep(u, processLimit, scale, otherArgs)
    
    delta = datetime.datetime.now() - tstart
    secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
    print "Total test time = %.2f secs." % secs

    grade.grade(gmeanDict, sys.stdout)

    if outFile:
        grade.grade(gmeanDict, outFile)
        outFile.close()

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
