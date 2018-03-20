#!/usr/bin/python

import subprocess
import sys
import os
import getopt
import math
import datetime
import socket
import random

import grade

def usage(fname):
    ustring = "Usage: %s [-h] [-a (x|o|n)] [-s SCALE] [-u UPDATELIST] [-f OUTFILE]" % fname
    ustring += " [-p PROCESSLIMIT]"
    print ustring 
    print "(All lists given as colon-separated text.)"
    print "    -h              Print this message"
    print "    -a (x|o|n)      Set processor affinity: x=None, o=Old flags, n=New flags"
    print "    -s SCALE        Reduce number of steps in each benchmark by specified factor"
    print "    -u UPDATELIST   Specify update modes(s):"
    print "       r: rat order"
    print "       s: synchronous"
    print "       b: batch"
    print "    -f OUTFILE      Create output file recording measurements"
    print "         If file name contains field of form XX..X, will replace with ID having that many digits"
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

# Additional flags to control MPI 
mpiCmd = ["mpirun"]

newMpiFlags = ["-map-by", "core", "-bind-to", "core"]
oldMpiFlags = ["-bycore", "-bind-to-core"]

# Dictionary of geometric means, indexed by (mode, threads)
gmeanDict = {}

# Host IP addresses on Latedays cluster, and their measured performance
hostDict = {
    '10.22.1.240' : 'fast',
    '10.22.1.241' : 'slow',
    '10.22.1.242' : 'slow',
    '10.22.1.243' : 'normal',
    '10.22.1.244' : 'normal',
    '10.22.1.245' : 'slow',
    '10.22.1.246' : 'fast',
    '10.22.1.247' : 'normal',
    '10.22.1.248' : 'normal',
    '10.22.1.249' : 'normal',
    '10.22.1.250' : 'slow',
    '10.22.1.251' : 'normal',
    '10.22.1.252' : 'normal',
    '10.22.1.253' : 'normal',
    '128.2.204.183' : 'normal'  # Latedays head node
    }

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

def cmd(graphSize, graphType, ratType, loadFactor, stepCount, updateType, processCount, mpiFlags, otherArgs = []):
    global bcount, logSum
    global cacheKey
    # File number of standard output
    stdoutFileNumber = 1
    updateFlag = UpdateMode.flags[updateType]
    params = ["%5d" % graphSize, graphType, "%4d" % loadFactor, ratType, str(stepCount), updateFlag]
    cacheKey = ":".join(params)
    results = params + [str(processCount)]
    sizeName = str(graphSize)
    graphFileName = dataDir + "g-" + graphType + sizeName + ".gph"
    ratFileName = dataDir + "r-" + sizeName + '-' + ratType + str(loadFactor) + ".rats"
    clist = runFlags + ["-g", graphFileName, "-r", ratFileName, "-u", updateFlag, "-n", str(stepCount), "-i", str(stepCount)] + otherArgs
    if processCount > 1:
        gcmd = mpiCmd + mpiFlags + ["-np", str(processCount), mpiSimProg] + clist
    else:
        gcmd = [simProg] + clist
    gcmdLine = " ".join(gcmd)


    # Legacy code to allow multiple trials
    reps = 1
    secs = 1e9
    for r in range(reps):
        tstart = datetime.datetime.now()
        try:
            simProcess = subprocess.Popen(gcmd, stderr = stdoutFileNumber)
            simProcess.wait()
            retcode = simProcess.returncode
        except Exception as e:
            outmsg("Execution of command '%s' failed. %s" % (gcmdLine, e))
            return False
        if retcode == 0:
            delta = datetime.datetime.now() - tstart
            tsecs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
            secs = min(secs, tsecs)
        else:
            outmsg("Execution of command '%s' gave return code %d" % (gcmdLine, retcode))
            return False

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
    return True

def sweep(updateType, processLimit, scale, mpiFlags, otherArgs):
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
        for bparams in benchmarkList:
            (graphSize, graphType, ratType, loadFactor) = bparams
            cmd(graphSize, graphType, ratType, loadFactor, stepCount, updateType, processCount, mpiFlags, otherArgs)
        if bcount > 0:
            gmean = math.exp(logSum/bcount)
            updateFlag = UpdateMode.flags[updateType]
            gmeanDict[(updateFlag, processCount)] = gmean
            outmsg(marker + "Gmean\t\t\t\t\t%s\t%d\t\t%7.2f" % (updateFlag, processCount, gmean))
            outmsg(marker + "---------" * 8)

def classifyProcessor():
    hostName = socket.gethostname()
    hostIP = socket.gethostbyname(hostName)
    status = hostDict[hostIP] if hostIP in hostDict else 'unknown'
    if status == 'unknown':
        outmsg("WARNING.  Running on unknown host %s (%s)" % (hostName, hostIP))
    elif status != 'normal':
        outmsg("WARNING.  Running on host %s (%s), which is known to be %s" % (hostName, hostIP, status))
    else:
        outmsg("INFO.  Running on normal host %s (%s)" % (hostName, hostIP))
    return "Host %s(%s)" % (hostIP, status)
               

def generateFileName(template):
    n = len(template)
    ls = []
    for i in range(n):
        c = template[i]
        if c == 'X':
            c = chr(random.randint(ord('0'), ord('9')))
        ls.append(c)
    return "".join(ls)

def run(name, args):
    global outFile
    scale = 1
    updateList = [UpdateMode.batch, UpdateMode.synchronous]
    optString = "ha:s:u:p:f:"
    processLimit = 100
    otherArgs = []
    mpiFlags = []
    optlist, args = getopt.getopt(args, optString)
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
        elif opt == '-s':
            scale = float(val)
        elif opt == '-a':
            if val == 'x':
                mpiFlags = []
            elif val == 'o':
                mpiFlags = oldMpiFlags
            elif val == 'n':
                mpiFlags = newMpiFlags
            else:
                outmsg("Invalid MPI flag specifier '%s'" % val)
                usage(name)
        elif opt == '-f':
            fname = generateFileName(val)
            try:
                outFile = open(fname, "w")
                outmsg("Writing to file '%s'" % fname)
            except Exception as e:
                outFile = None
                outmsg("Couldn't open file '%s'" % fname)
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
                    outmsg("Invalid update mode '%s'" % c)
                    usage(name)
        elif opt == '-p':
            processLimit = int(val)
        else:
            outmsg("Uknown option '%s'" % opt)
            usage(name)

    hostInfo = classifyProcessor()
    tstart = datetime.datetime.now()

    for u in updateList:
        sweep(u, processLimit, scale, mpiFlags, otherArgs)
    
    delta = datetime.datetime.now() - tstart
    secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
    outmsg("Total test time = %.2f secs." % secs)

    grade.grade(gmeanDict, sys.stdout, hostInfo)

    if outFile:
        grade.grade(gmeanDict, outFile, hostInfo)
        outFile.close()

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
