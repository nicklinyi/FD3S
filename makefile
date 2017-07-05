all:
	mpif90.mpich -mismatch -w fd3s_modules.f90 fd3s_comm_tap.f90 fd3s_evolution.f90 fd3s_init.f90 fd3s_input.f90 fd3s_main.f90 fd3s_oper.f90 fd3s_output.f90 -o main.exe 
