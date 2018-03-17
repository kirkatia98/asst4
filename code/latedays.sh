#!/bin/bash
# This script lets you submit jobs for execution on the latedays cluster
# You should submit it using qsub:
#   'qsub latedays.sh'

# Upon completion, the output generated on stdout will show up in the
# file latedays.sh.oNNNNN where NNNNN is the job number.  The output
# generated on stderr will show up in the file latedays.sh.eNNNNN.

# Limit execution time to 20 minutes
#PBS -lwalltime=0:20:00
# Allocate all available CPUs on a single node
#PBS -l nodes=1:ppn=24

# Go to the directory from which you submitted your job
cd $PBS_O_WORKDIR

# Execute the performance evaluation program and store summary in benchmark-XXXX.out
# (XXXX will get replaced by randomly generated 4-digit sequence)
# The "-a o" option instructs it to use processor affinity directives that are compatible
# with the version of MPI running on the Latedays nodes
./benchmark.py -a o -f benchmark-XXXX.out


