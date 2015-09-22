CC = g++ -std=c++11 
OPTS = -O3
CXXFLAGS+=-Wall -DNDEBUG $(OPTS)


../CompositeDiscover: CompositeDiscover
	cp CompositeDiscover ../

CompositeDiscover: CompositesDiscover.cpp CompositeFind.cpp Composite.cpp Family.cpp libseq.a libsequence.a Composite.h CompositeFind.h FamilyPath.h 
	${CC} ${CXXFLAGS} CompositesDiscover.cpp CompositeFind.cpp Composite.cpp Family.cpp -o CompositeDiscover -L . -lseq -lsequence

libseq.a: seq.c nib.c util.c charvec.c
	gcc -O3 -c seq.c nib.c util.c charvec.c ${OPTS}
	ar rc libseq.a seq.o nib.o util.o charvec.o 
	ranlib libseq.a
	rm seq.o nib.o util.o charvec.o 

libsequence.a: sequence.cpp
	g++ -O3 -c sequence.cpp ${OPTS}
	ar rc libsequence.a sequence.o
	ranlib libsequence.a
	rm sequence.o

clean:
	rm libseq.a libsequence.a CompositeDiscover
