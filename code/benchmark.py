#!/usr/bin/python

import subprocess
import sys
import os
import os.path
import getopt
import math
import datetime
import socket
import random

import grade

def usage(fname):
    ustring = "Usage: %s [-h] [-a (x|o|n)] [-s SCALE] [-u UPDATELIST] [-f OUTFILE]" % fname
    ustring += "[-c] [-p PROCESSLIMIT]"
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
    print "    -c            Compare simulator output to recorded result"
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

dataDirectory = "./data/"
outFile = None
captureDirectory = "./capture"
doCheck = False

# How many mismatched lines warrant detailed report
mismatchLimit = 5

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
captureRunFlags = []

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


def captureFileName(graphSize, graphType, ratType, loadFactor, stepCount, updateFlag):
    params = (graphSize, graphType, ratType, loadFactor, stepCount, updateFlag)
    return captureDirectory + "/cap" + "-%.3d-%s-%s-%.3d-%.3d-%s.txt" % params

def openCaptureFile(graphSize, graphType, ratType, loadFactor, stepCount, updateFlag):
    if not doCheck:
        return None
    name = captureFileName(graphSize, graphType, ratType, loadFactor, stepCount, updateFlag)
    try:
        cfile = open(name, "r")
    except Exception as e:
        outmsg("Couldn't open captured result file '%s': %s" % (name, e))
        return None
    return cfile


def checkOutputs(captureFile, outputFile):
    if captureFile == None or outputFile == None:
        return True
    badLines = 0
    lineNumber = 0
    while True:
        rline = captureFile.readline()
        tline = outputFile.readline()
        lineNumber +=1
        if rline == "":
            if tline == "":
                break
            else:
                badLines += 1
                outmsg("Mismatch at line %d.  Reference file ended prematurely" % (lineNumber))
                break
        elif tline == "":
            badLines += 1
            outmsg("Mismatch at line %d.  Simulation output ended prematurely\n" % (lineNumber))
            break
        if rline[-1] == '\n':
            rline = rline[:-1]
        if tline[-1] == '\n':
            tline = tline[:-1]
        if rline != tline:
            badLines += 1
            if badLines <= mismatchLimit:
                outmsg("Mismatch at line %d.  Expected result:'%s'.  Simulation result:'%s'\n" % (lineNumber, rline, tline))
    captureFile.close()
    if badLines > 0:
        outmsg("%d total mismatches.\n" % (badLines))
    if badLines == 0:
        outmsg("Simulator output matches recorded results!")
    return badLines == 0


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
    graphFileName = dataDirectory + "g-" + graphType + sizeName + ".gph"
    ratFileName = dataDirectory + "r-" + sizeName + '-' + ratType + str(loadFactor) + ".rats"
    checkFile = openCaptureFile(graphSize, graphType, ratType, loadFactor, stepCount, updateFlag)
    recordOutput = checkFile is not None
    ok = True
    if recordOutput:
        clist = captureRunFlags + ["-g", graphFileName, "-r", ratFileName, "-u", updateFlag, "-n", str(stepCount), "-i", str(stepCount)] + otherArgs
    else:
        clist = runFlags + ["-g", graphFileName, "-r", ratFileName, "-u", updateFlag, "-n", str(stepCount), "-i", str(stepCount)] + otherArgs
    if processCount > 1:
        gcmd = mpiCmd + mpiFlags + ["-np", str(processCount), mpiSimProg] + clist
    else:
        gcmd = [simProg] + clist
    gcmdLine = " ".join(gcmd)
    retcode = 1
    tstart = datetime.datetime.now()
    try:
        if recordOutput:
            simProcess = subprocess.Popen(gcmd, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
            ok = ok and checkOutputs(checkFile, simProcess.stdout)
            # Echo any results printed by simulator on stderr onto stdout
            for line in simProcess.stderr:
                sys.stdout.write(line)
        else:
            simProcess = subprocess.Popen(gcmd, stderr = stdoutFileNumber)
        simProcess.wait()
        retcode = simProcess.returncode
    except Exception as e:
        outmsg("Execution of command '%s' failed. %s" % (gcmdLine, e))
        return False
    if retcode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
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
    return ok

def sweep(updateType, processLimit, scale, mpiFlags, otherArgs):
    runList = synchRunList if updateType == UpdateMode.synchronous else otherRunList
    ok = True
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
            ok = ok and cmd(graphSize, graphType, ratType, loadFactor, stepCount, updateType, processCount, mpiFlags, otherArgs)
        if bcount > 0:
            gmean = math.exp(logSum/bcount)
            updateFlag = UpdateMode.flags[updateType]
            gmeanDict[(updateFlag, processCount)] = gmean
            outmsg(marker + "Gmean\t\t\t\t\t%s\t%d\t\t%7.2f" % (updateFlag, processCount, gmean))
            outmsg(marker + "---------" * 8)
    return ok

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
    global outFile, doCheck
    scale = 1
    updateList = [UpdateMode.batch, UpdateMode.synchronous]
    optString = "ha:s:u:p:f:c"
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
        elif opt == '-c':
            doCheck = True
        elif opt == '-p':
            processLimit = int(val)
        else:
            outmsg("Unknown option '%s'" % opt)
            usage(name)

    hostInfo = classifyProcessor()
    tstart = datetime.datetime.now()
    ok = True

    for u in updateList:
        ok = ok and sweep(u, processLimit, scale, mpiFlags, otherArgs)
    
    delta = datetime.datetime.now() - tstart
    secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
    outmsg("Total test time = %.2f secs." % secs)

    grade.grade(ok, gmeanDict, sys.stdout, hostInfo)

    if outFile:
        grade.grade(ok, gmeanDict, outFile, hostInfo)
        outFile.close()

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
