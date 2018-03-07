This directory contains core code for GraphRat simulators written in Python and C.

Quickstart:

You can run a series of demos by executing:

    linux> make demoX

for X from 1 to 10.

Executable Files:
	grun.py	      Simulator.  Can also operate as visualizer for another simulator
	regress.py    Regression test C version of simulator against Python version.
	benchmark.py  Benchmark C programs and report grades

Python support Files:
	gengraph.py   Used by grun.py to load graphs
	grade.py      Implements grading logic
	rutil.py      Support for random number generation and value function calculation.
	sim.py        Core simulator implementation
	viz.py        Support for visualization of graphs using ASCII formatting and/or a heat-map representation
	
C Files:
	crun.{h,c}    Top-level control for simulator
	sim.c         Core simulation code
	simutil.c     Routines for supporting simulation
	rutil.{h,c}   Support for random number generation and value function calculation.
	cycletimer.{h,c} Implements low-overhead, fine-grained time measurements

Other Files:
        latedays.sh   Used to submit benchmarking jobs when using the Latedays cluster

FILE/OUTPUT FORMATS

All files are line-oriented text files, using decimal representations
of numbers.  Any line having the first non-whitespace character equal
to '#' is ignored.

GRAPH FILES

First line of form "N M T" where N is number of nodes, M is number of
edges, and T is the maximum tile size.  (this means that only grid
connections can occur between every T'th pair of rows.)

Remaining lines of form "I J", indicating an edge from node I to node
J.  The graphs are all undirected, and so their will be entries for
both (I, J) and (J, I).

Nodes must be numbered between 0 and N-1.  Edges must be in sorted
order, with the sorting key first in I and then in J.

RAT POSITION FILES

First line of form "N R" where N is number of nodes, and R is number of rats

Remaining lines of form "I", indicating node number of each successive rat.
I must be between 0 and N-1.

SIMULATION DRIVER

When operating in driving mode the simulator should produce the following on each step:

First line: "STEP N R", where N is the number of nodes, and R is the number of rats
Optionally:
	R lines of the form "I", giving position of each successive rat
	This information can be omitted to reduce bandwidth requirements.  Display will remain at previous state.
Last line for step: "END"

At the very end, the final line of the stream should be "DONE"

Note: Don't try to print error messages or debugging information for
the simulator on stdout, since this will be piped to grun.py.
Instead, use stderr.  If you need to perform error exit, emit "DONE"
on stdout to terminate visualization.

