# ------------------------------------------------
# Please set your own MPICH2 Directory (MPICH_DIR)
# example:
# MPICH_DIR = /home/nick/mpich2

MPICH_DIR = 

# ------------------------------------------------
PROG = $MPICH_DIR/bin/mpif90
PROGRUN = $MPICH_DIR/bin/mpiexec
FFLAGS = -I$MPICH_DIR/include -L$MPICH_DIR/lib/ -lmpi


OBJECT1 = generate_models.o 
OBJECT2 = generate_initial.o
NUM_PROCS=18

# run mpirun
#mpirun -np $NUMBER_OF_PROCESSORS  $PROGRAM

test1: generate_model
	$(PROGRUN) -n $(NUM_PROCS) ./generate_model

generate_model: $(OBJECT1)
	$(PROG) $(FFLAGS) $< -o $@

generate_initial: $(OBJECT2)
	$(PROG) $(FFLAGS) $< -o $@


%.o : %.f90
	$(PROG) $(FFLAGS) -c  $< 
