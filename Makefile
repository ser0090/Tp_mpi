CC=mpicc
CODE=main.c
RUNNER=mpirun
PROC_COUNT=4

all: $(CODE)
	$(CC) -o main $(CODE)

run:
	$(RUNNER) -np $(PROC_COUNT) ./main

