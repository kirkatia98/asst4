FILES= code/

XFILE = code/data

handin.tar: clean
	tar cvf handin.tar --exclude=$(XFILE) $(FILES) 

clean:
	(cd code; make clean)
	rm -f *~ handin.tar
