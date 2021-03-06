# ------------------------------------------------
# Please set your own MPICH2 Directory (MPICH_DIR)
# example:
# MPICH_DIR = /home/nick/mpich2

MPICH_DIR = 

# ------------------------------------------------
PROG = $MPICH_DIR/bin/mpif90
PROGRUN = $MPICH_DIR/bin/mpiexec
FFLAGS = -I$MPICH_DIR/include -L$MPICH_DIR/lib/ -lmpi
NUM_PROCS=18
OBJECTS = fd3s_modules.o fd3s_comm_tap.o fd3s_evolution.o fd3s_init.o fd3s_input.o fd3s_main.o  fd3s_oper.o fd3s_output.o

test: main.exe
	$(PROGRUN) -n $(NUM_PROCS) ./main.exe

main.exe: $(OBJECTS)
	$(PROG) $(FFLAGS) $(OBJECTS) -o $@

%.o : %.f90
	$(PROG) $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod
