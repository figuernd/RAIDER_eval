# -*- makefile -*-
# just make headings for each thing and put the headings
# on the line for ALL
PROJECTDIR = .
BINDIR = $(PROJECTDIR)/bin

CC = g++
CFLAGS =  -g 
LDFLAGS = -lm


CLEANABLE_FILES = *.o *~

ALL = ABclass_search BBclass_search GBBclass_search

.SUFFIXES: .cc

.cc.o:
	$(CC) $(CFLAGS) -c $*.cc -o $*.o 

all:    $(ALL)

ABclass_search :    ABclass_search.o
	$(CC) -o $@ ABclass_search.o $(LDFLAGS)

BBclass_search :    BBclass_search.o
	$(CC) -o $@ BBclass_search.o $(LDFLAGS)

GBBclass_search :    GBBclass_search.o
	$(CC) -o $@ GBBclass_search.o $(LDFLAGS)

clean:
	/bin/rm -f $(CLEANABLE_FILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.

