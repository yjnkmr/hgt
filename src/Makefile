### Makfile for HGT analysis

CC = gcc
CFLAGS = -lstdc++ -lm
PRGS = predict markov
INSTALLDIR = ../bin
FLAG = -D SKIP

.PHONY: all
all: $(PRGS)

# -----  Main Program -----
#
predict: main.o hmm.o viterbi.o readfile.o SEQ.o
	$(CC) -O -o predict main.o hmm.o viterbi.o readfile.o SEQ.o $(CFLAGS)
# ----- Model Construction Programs -----
markov: markov.o readfile.o SEQ.o
	$(CC) -O -o markov markov.o readfile.o SEQ.o $(CFLAGS)

main.o: main.cpp hmm.h viterbi.h readfile.h SEQ.h
	$(CC) $(DEBUG) -c main.cpp
markov.o: markov.cpp readfile.h SEQ.h
	$(CC) -c markov.cpp

hmm.o: hmm.cpp hmm.h readfile.h
	$(CC) $(DEBUG) -c hmm.cpp
viterbi.o: viterbi.cpp viterbi.h hmm.h
	$(CC) -c viterbi.cpp $(FLAG)
readfile.o: readfile.cpp readfile.h
	$(CC) -c readfile.cpp
SEQ.o: SEQ.cpp SEQ.h
	$(CC) -c SEQ.cpp

install: $(PRGS)
	install -s predict $(INSTALLDIR)
	install -s markov $(INSTALLDIR)

clean:
	rm -f core *.o *~ $(PRGS)
