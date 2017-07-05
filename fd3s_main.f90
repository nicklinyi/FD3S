!*****************************************************************************
!***************** Main Program for spherical FD code ************************
!*****************************************************************************

program fd3s_main
use parameters
use variables
implicit none
include 'mpif.h'

    	integer, dimension(8) :: val
        character(len=8) :: date1,date2
	character(len=10) :: time1,time2  
	character(len=20) :: dummy, message
	real :: vr_min, vr_max
    
	!======================================================================
	! initialise mpi
	!======================================================================
	
	call mpi_init(ierr)
	call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
	call mpi_comm_size(mpi_comm_world, p, ierr)
    
	!======================================================================
	! open and write logfile
	!======================================================================
	
	call int2str(my_rank, dummy)
    
	open(unit=99, file='DATA/LOGFILES/logfile'//dummy(1:len_trim(dummy)), action='WRITE')
    
	write(99,*) '---------------------------------------'
	write(99,*) 'logfile process ', my_rank
	write(99,*) '---------------------------------------'
    
	if (my_rank==0) then
	
		call date_and_time(date1, time1)
		
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'starting date: ', date1(7:8), '. ', date1(5:6), '. ', date1(5:5), date1(4:4)
		write(*,*) 'starting time: ', time1(1:2), ':', time1(3:4), ',', time1(5:6)
		write(*,*) '------------------------------------------------------------------'
		
		call tic(val)
		
	endif
    
	!======================================================================
	! read parameters from files 'Par' and 'recfile' and 'boxfile
	!======================================================================

	call fd3s_input

	!======================================================================
	! various initializations (model space, source, receivers, boundaries)
	!======================================================================
	
	call fd3s_init
	
	call comm_tap_v
	
	!======================================================================
	! start time evolution
	!======================================================================

    	do it=1,nt
		
		!==============================================================
		! calculate FD solutions for wavelet source
		!==============================================================
		
		write(99,*) "iteration", it
		
		call fd3s_evolution
		call record_seismograms
		
		call mpi_reduce(maxval(w11),vr_max,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
		call mpi_reduce(minval(w11),vr_min,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
		
		if (my_rank==0) then
			write(*,*) 'iteration ',it,':   ', vr_min, '< vr <', vr_max
		endif
		message="elapsed time"
		call write_cpu_time_log(message)
		write(99,*) 'iteration ',it,':   ', vr_min, '< vr <', vr_max
		
		if (mod(it,ssamp)==0) then
			
			call fd3s_output(it)
		
		endif
		
	enddo
   
	call fd3s_output(nt)
	
	!======================================================================
	! date and time output
	!======================================================================

	if (my_rank==0) then
	
		call date_and_time(date2, time2)
		
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'starting date: ', date1(7:8), '. ', date1(5:6), '. ', date1(5:5), date1(4:4)
		write(*,*) 'starting time: ', time1(1:2), ':', time1(3:4), ',', time1(5:6)
		write(*,*) '------------------------------------------------------------------'
		write(*,*) 'finishing date: ', date2(7:8), '. ', date2(5:6), '. ', date2(5:5), date2(4:4)
		write(*,*) 'finishing time: ', time2(1:2), ':', time2(3:4), ',', time2(5:6)
		write(*,*) '------------------------------------------------------------------'
	  
		call toc(val)
		
	endif
	
!==============================================================================
! clean up
!==============================================================================

	close(99)
	call mpi_finalize(ierr)

END PROGRAM fd3s_main


!########### Additional subroutines ###########################################

subroutine tic(value)
implicit none

	integer, dimension(8), intent(out) :: value
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: time_zone
	
	call date_and_time(date, time, time_zone, value)

end subroutine tic

subroutine toc(value1)
implicit none

	integer, dimension(8), intent(in) :: value1
	integer, dimension(8) :: value2
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: time_zone
	real, dimension(1) :: seconds1, seconds2
	
	call date_and_time(date, time, time_zone, value2)
	
	seconds1(1:1)=real(value1(7:7))+real(value1(8:8))/1000.0+60.0*real(value1(6:6))+ &
		      60.0*60.0*real(value1(5:5))
	
	seconds2(1:1)=real(value2(7:7))+real(value2(8:8))/1000.0+60.0*real(value2(6:6))+ &
		      60.0*60.0*real(value2(5:5))

	write(*,*) 'elapsed time: ', seconds2(1:1)-seconds1(1:1), 'seconds'
	!write(*,*) 'elapsed time: ', seconds1(1:1), 'seconds'

end subroutine toc	


subroutine write_cpu_time(proc)
use parameters
use variables
implicit none
include 'mpif.h'

	real :: tmd
	integer, intent(in) :: proc

	if (my_rank==proc) then
		call cpu_time(tmd)
		write(*,*) tmd
	endif
	
end subroutine write_cpu_time

subroutine write_cpu_time_log(message)
use parameters
use variables
implicit none
include 'mpif.h'

	real :: tmd
	character(len=20), intent(in) :: message

	call cpu_time(tmd)
	write(99,*) message
	write(99,*) tmd
	
end subroutine write_cpu_time_log

!==============================================================================
! subroutine to add variable integers to a character string
!==============================================================================

SUBROUTINE add_num_to_char(file,i)
    
        IMPLICIT NONE
        CHARACTER(LEN=50) file
        CHARACTER(LEN=4) ichar

        INTEGER i,in

        IF(i.LT.10)THEN
          WRITE(ichar,'(i1)') i
          in=1
        ELSEIF(i.LT.100)THEN
          WRITE(ichar,'(i2)') i
          in=2
        ELSEIF(i.LT.1000)THEN
          WRITE(ichar,'(i3)') i
          in=3
        ELSEIF(i.LT.10000)THEN
          WRITE(ichar,'(i4)') i
          in=4
        ELSE
          WRITE(*,*)' Warning : number too large !!! '
        ENDIF
 
        file = file(1:LEN_TRIM(file))//ichar(1:in)
 
        RETURN
END SUBROUTINE add_num_to_char
