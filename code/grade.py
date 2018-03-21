import math

# Grading targets for GraphRats benchmark

# Targets for speed (in MRPS)
gmeanTarget = {
    ('b', 1) :      25.0,
    ('s', 1) :      35.0,
    ('b', 12) :     87.5,
    ('s', 12) :    112.0,
    }

# Targets for speedup
speedupTarget = {
    ('b', 1, 12) : 3.5,
    ('s', 1, 12) : 3.2,
}

modes = ['b', 's']
parallel = [1, 12]

flagNames = {'b' : "batch", 's' : "synch"}


fullCreditThreshold = 1.0
partialCreditThreshold = 2.0

perfWeight = 17.0
speedupWeight = 6.0

def score(s, target, wt):
    mins = target/partialCreditThreshold
    maxs = target/fullCreditThreshold
    if s >= maxs:
        return wt
    if s < mins:
        return 0
    return wt * (s-mins)/(maxs-mins)

# Accept dictionary of MRPS measurements, each having (mode, parallel) as key
def grade(ok, gmeanDict, outf, info = ""):
    total = 0.0
    maxtotal = 0.0
    outf.write("---------" * 9 + "\n")
    outf.write("MRPS Scores\n")
    for t in parallel:
        for m in modes:
            wt = perfWeight
            if (m,t) not in gmeanDict:
                continue
            s = gmeanDict[(m,t)]
            maxtotal += wt
            name = flagNames[m]
            target = gmeanTarget[(m,t)]
            val = score(s, target, wt)
            outString = "   "
            if info != "":
                outString += info + ", "
            outString += "Parallel = %d, Mode = %s, Achieved = %.2f, Target = %.2f, Score = %.2f/%.2f\n" % (t, name, s, target, val, wt)
            outf.write(outString)
            total += val
    outf.write("Speedup Scores\n")
    params = speedupTarget.keys()
    params.sort()
    for (m, t1, t2) in params:
        wt = speedupWeight
        if (m, t1) not in gmeanDict or (m, t2) not in gmeanDict:
            continue
        g1 = min(gmeanDict[(m,t1)], gmeanTarget[(m,t1)])
        s = float(gmeanDict[(m,t2)]) / g1
        maxtotal += wt
        name = flagNames[m]
        target = speedupTarget[(m, t1, t2)]
        val = score(s, target, wt)
        outf.write("  Ratio = %d:%d, Mode = %s, Achieved = %.2f, Target = %.2f, Score = %.2f/%.2f\n" % (t2, t1, name, s, target, val, wt))
        total += val
    itotal = math.ceil(total)
    if not ok:
        outf.write("ERROR: One or more tests failed.  No credit given\n")
        itotal = 0
    outf.write("TOTAL = %d/%.0f\n" % (itotal, maxtotal))
    outf.write("---------" * 9 + "\n")

    
            
    
        
    
    
