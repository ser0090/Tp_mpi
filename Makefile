CC=mpicc
CODE=main.c
RUNNER=mpirun
PROC_COUNT=4
CFLAGS=-fopenmp

all: $(CODE)
	$(CC) $(CFLAGS) -o main $(CODE)

run:
	$(RUNNER) -np $(PROC_COUNT) $(RUNNER_FLAGS) ./main

