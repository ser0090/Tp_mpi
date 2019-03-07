CC=mpicc
MAIN_BLOQ=main.c
EXE_BLOQ=main
MAIN_NBLOQ=main_asinc.c
EXE_NBLOQ=main_asinc
RUNNER=mpirun
PROC_COUNT=4
CFLAGS=-fopenmp

all: $(CODE)
	$(CC) $(CFLAGS) $(OTHER_CFLAGS) -o $(EXE_BLOQ)  $(MAIN_BLOQ)
	$(CC) $(CFLAGS) $(OTHER_CFLAGS) -o $(EXE_NBLOQ)  $(MAIN_NBLOQ)
run_sinc:
	$(RUNNER) -np $(PROC_COUNT) $(RUNNER_FLAGS) $(EXE_BLOQ)
run_asinc:
	$(RUNNER) -np $(PROC_COUNT) $(RUNNER_FLAGS) $(EXE_NBLOQ)
