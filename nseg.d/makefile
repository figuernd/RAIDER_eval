
all : nseg nmerge

nseg : nseg.c lnfac.h genwin.h genwin.o
	gcc -o nseg nseg.c genwin.o -lm

nmerge : nmerge.c genwin.h genwin.o
	gcc -o nmerge nmerge.c genwin.o -lm

genwin.o : genwin.c genwin.h
	gcc -c genwin.c

clean:
	rm -f nseg nseg.o genwin.o
