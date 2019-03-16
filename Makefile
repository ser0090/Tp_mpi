CC=mpicc
RUNNER=mpirun
PROC_COUNT=4
CFLAGS=-fopenmp
EXE_SINC=main_sinc
EXE_ASINC=main_asinc
CODE_SINC=main.c
CODE_ASINC=main_asinc.c


all: $(CODE)
	$(CC) $(CFLAGS) -o $(EXE_SINC) $(OTHER_CFLAGS) $(CODE_SINC)
	$(CC) $(CFLAGS) -o $(EXE_ASINC) $(OTHER_CFLAGS) $(CODE_ASINC)

run_sinc:
	$(RUNNER) -np $(PROC_COUNT) $(RUNNER_FLAGS) $(EXE_SINC)

run_asinc:
	$(RUNNER) -np $(PROC_COUNT) $(RUNNER_FLAGS) $(EXE_ASINC)


