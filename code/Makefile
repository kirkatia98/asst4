CC=gcc

MPI=-DMPI
MPICC = mpicc

DEBUG=0
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
LDFLAGS= -lm
DDIR = ./data

CFILES = crun.c graph.c simutil.c sim.c rutil.c cycletimer.c
HFILES = crun.h rutil.h cycletimer.h

GFILES = gengraph.py grun.py rutil.py sim.py viz.py  regress.py benchmark.py grade.py


DFILES = $(DDIR)/g-t3600.gph $(DDIR)/g-t32400.gph $(DDIR)/g-t4.gph $(DDIR)/g-t400.gph \
	 $(DDIR)/g-u3600.gph $(DDIR)/g-u32400.gph $(DDIR)/g-u4.gph $(DDIR)/g-u400.gph 

RFILES  = $(DDIR)/r-3600-d10.rats  $(DDIR)/r-3600-u10.rats \
	$(DDIR)/r-32400-d32.rats  $(DDIR)/r-32400-u32.rats \
	$(DDIR)/r-4-d1.rats  $(DDIR)/r-4-u1.rats \
	$(DDIR)/r-400-d10.rats $(DDIR)/r-400-u10.rats 

all: crun crun-mpi

crun: crun-seq
	cp -p crun-seq crun

crun-seq: $(CFILES) $(HFILES) 
	$(CC) $(CFLAGS) -o crun-seq $(CFILES) $(LDFLAGS)

crun-mpi: $(CFILES) $(XCFILES) $(HFILES) $(XHFILES)
	$(MPICC) $(CFLAGS) $(MPI) -o crun-mpi $(CFILES) $(XCFILES) $(LDFLAGS)


demo1: grun.py
	@echo "Running Python simulator with text visualization.  Synchronous mode."
	./grun.py -g data/g-u400.gph -r data/r-400-d10.rats -n 20 -u s -v a -p 0.3

demo2: grun.py
	@echo "Running Python simulator with heat-map visualization.  Synchronous mode"
	@echo "Note the oscillatory behavior"
	./grun.py -g data/g-u400.gph -r data/r-400-d10.rats -u s -n 20 -v h -p 0.3

demo3: grun.py
	@echo "Running Python simulator with heat-map visualization.  Rat-order mode"
	@echo "Note the smooth transitions"
	./grun.py -g data/g-u400.gph -r data/r-400-d10.rats -u r -n 20 -v h -p 0.3

demo4: grun.py
	@echo "Running Python simulator with heat-map visualization.  Batch mode"
	@echo "Note the smooth transitions"
	./grun.py -g data/g-u400.gph -r data/r-400-d10.rats -u b -n 20 -v h -p 0.3

demo5: crun grun.py
	@echo "Running C simulator using Python simulator as visualizer"
	./crun -g data/g-u400.gph -r data/r-400-d10.rats -u b -n 20 | ./grun.py -d -v h -p 0.3

demo6: crun grun.py
	@echo "Running on a 180x180 grid graph"
	@echo "This will run a lot faster once you speed up your code"
	./crun -g data/g-u32400.gph -r data/r-32400-d32.rats -u b -n 100 | ./grun.py -d -v h -p 0.02

demo7: crun grun.py
	@echo "Running on a 180x180 tiled graph"
	./crun -g data/g-t32400.gph -r data/r-32400-d32.rats -u b -n 100 | ./grun.py -d -v h -p 0.02

demo8: crun grun.py
	@echo "Running on a 180x180 tiled graph with rats initially distributed uniformly"
	./crun -g data/g-t32400.gph -r data/r-32400-u32.rats -u b -n 100 | ./grun.py -d -v h -p 0.02

demo9: crun grun.py
	@echo "Running on a 180x180 grid graph with rats initially distributed uniformly"
	./crun -g data/g-u32400.gph -r data/r-32400-u32.rats -u b -n 100 | ./grun.py -d -v h -p 0.02

clean:
	rm -f *~ *.pyc
	rm -rf *.dSYM
	rm -f *.tgz
	rm -f crun crun-seq crun-mpi
