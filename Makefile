CFLAGS= -O3 -Wall -Wextra
LDFLAGS= -lm
CC=gcc
SRC=viterbi.c hmmutils.c sequence.c forward.c backward.c baumwelch.c bfgs.c \
	nrutil.c hmmarc.c lbfgsb.c linesearch.c subalgorithms.c print.c linpack.c \
	miniCBLAS.c timer.c
OBJ=$(SRC:.c=.o)

all : hmmarc

hmmarc: hmmarc.o $(OBJ)
	 $(CC) -o hmmarc $(OBJ) $(CFLAGS) $(LDFLAGS)
clean:
	rm *.o hmmarc
