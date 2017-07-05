!=============================================================================
! Andreas Fichtner, 26. Mai 2006
!=============================================================================

PROGRAM generate_initial
IMPLICIT NONE
include 'mpif.h'

!=============================================================================
! variables 
!=============================================================================

	character(len=20) :: junk
	integer :: i, j, k, mi
	integer :: model_type
	integer :: nx, ny, nz, px, py, pz
	integer :: bnx, bny, bnz
	real :: dummy
	real :: z_max, z_min, y_max, y_min, x_max, x_min, dz, dy, dx
	real, allocatable, dimension(:,:,:) :: field
	integer, allocatable, dimension(:,:) :: min_index, max_index		! minimum and maximum index boundaries
	integer :: status, ierr, p, my_rank, is_aniso, is_formatted


!=============================================================================
! initialise mpi
!=============================================================================

	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
	
!=============================================================================
! process 0 reads input and broadcasts it to the other processes
!=============================================================================

	if (my_rank==0) then

		open(UNIT=10,FILE='model_parameters',ACTION='READ')
	
		read(10,*) junk
		read(10,*) nx
		!WRITE(*,*) nx
		read(10,*) ny
		!WRITE(*,*) ny
		read(10,*) nz
		!WRITE(*,*) nz
		read(10,*) x_min
		!write(*,*) x_min
		read(10,*) x_max
		!write(*,*) x_max
		read(10,*) y_min
		!write(*,*) y_min
		read(10,*) y_max
		!write(*,*) y_max
		read(10,*) z_min
		!WRITE(*,*) z_min
		read(10,*) z_max
		!WRITE(*,*) z_max
	
		read(10,*) junk
		read(10,*) model_type
		!WRITE(*,*) model_type
		read(10,*) dummy
		!WRITE(*,*) mu_homogen
		read(10,*) dummy
		!WRITE(*,*) lambda_homogen
		read(10,*) dummy
		!WRITE(*,*) rho_homogen
		read(10,*) dummy
		!WRITE(*,*) A_homogen
		read(10,*) dummy
		!WRITE(*,*) B_homogen
		read(10,*) dummy
		!WRITE(*,*) C_homogen
		read(10,*) is_aniso
		read(10,*) is_formatted
	
		read(10,*) junk
		read(10,*) px
		!write(*,*) px
		read(10,*) py
		!write(*,*) py
		read(10,*) pz
		!write(*,*) pz
	
		CLOSE(UNIT=10)
		
	endif
	
	call mpi_bcast(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
	call MPI_BCAST(x_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(x_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(y_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(y_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(z_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(z_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	
	call MPI_BCAST(px,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(py,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(is_formatted,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
	!=====================================================================
	! set up box geometry 
	!=====================================================================

	allocate(min_index(1:px*py*pz,1:3),STAT=status)
	allocate(max_index(1:px*py*pz,1:3),STAT=status)
	
	dz=(z_max-z_min)/nz
	dy=(y_max-y_min)/ny
	dx=(x_max-x_min)/nx
	
	bnx=floor(real((nx+1)/px))
	bny=floor(real((ny+1)/py))
	bnz=floor(real((nz+1)/pz))
	
	do k=1,pz,1
		do j=1,py,1
			do i=1,px,1
				
				! indices (multiple and single) of the processor boxes
				mi=(k-1)*py*px+(j-1)*px+i	! long index of a processor box
				!write(*,*) i,j,k,mi
				
				! index boundaries of the processor boxes
				min_index(mi,1)=(i-1)*bnx+1-int(real(1/i))
				min_index(mi,2)=(j-1)*bny+1-int(real(1/j))
				min_index(mi,3)=(k-1)*bnz+1-int(real(1/k))
				
				if (i/=px) then
					max_index(mi,1)=i*bnx
				else
					max_index(mi,1)=nx
				!write(*,*) 'hello from process ', my_rank, ', max_index(x)=', max_index(mi,1)
				endif
				
				if (j/=py) then
					max_index(mi,2)=j*bny
				else
					max_index(mi,2)=ny
				endif
				
				if (k/=pz) then
					max_index(mi,3)=k*bnz
				else
					max_index(mi,3)=nz
				endif
				
			enddo
		enddo
	enddo

	!==============================================================================
	! save to files, formatted or unformatted
	!==============================================================================

	call int2str(my_rank,junk)

	if (is_formatted==1) then
	
		OPEN(UNIT=11,FILE='v0_phi'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=12,FILE='v0_r'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=13,FILE='v0_theta'//junk(1:len_trim(junk)),ACTION='WRITE')
		
		OPEN(UNIT=21,FILE='s0_phi_phi'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=22,FILE='s0_phi_r'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=23,FILE='s0_theta_phi'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=24,FILE='s0_theta_r'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=25,FILE='s0_theta_theta'//junk(1:len_trim(junk)),ACTION='WRITE')
		OPEN(UNIT=26,FILE='s0_r_r'//junk(1:len_trim(junk)),ACTION='WRITE')
		
	else
		
		OPEN(UNIT=11,FILE='v0_phi'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=12,FILE='v0_r'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=13,FILE='v0_theta'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		
		OPEN(UNIT=21,FILE='s0_phi_phi'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=22,FILE='s0_phi_r'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=23,FILE='s0_theta_phi'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=24,FILE='s0_theta_r'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=25,FILE='s0_theta_theta'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		OPEN(UNIT=26,FILE='s0_r_r'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
		
	endif

	if (is_formatted==1) then
	
		DO k=min_index(my_rank+1,3), max_index(my_rank+1,3)
			DO j=min_index(my_rank+1,2), max_index(my_rank+1,2)
				DO i=min_index(my_rank+1,1), max_index(my_rank+1,1)
				
					WRITE(11,*) 0.00
					WRITE(12,*) 0.00
					WRITE(13,*) 0.00
					
					WRITE(21,*) 0.00
					WRITE(22,*) 0.00
					WRITE(23,*) 0.00
					WRITE(24,*) 0.00
					WRITE(25,*) 0.00
					WRITE(26,*) 0.00
				
				END DO
			END DO
		END DO
		
	else
		
		allocate(field(min_index(my_rank+1,1):max_index(my_rank+1,1),min_index(my_rank+1,2):max_index(my_rank+1,2), &
			(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status)
			
		field(:,:,:)=0.0
		
		write(11) field(:,:,:)
		write(12) field(:,:,:)
		write(13) field(:,:,:)
		
		write(21) field(:,:,:)
		write(22) field(:,:,:)
		write(23) field(:,:,:)
		write(24) field(:,:,:)
		write(25) field(:,:,:)
		write(26) field(:,:,:)
		
		deallocate(field)
		
	endif
		
	CLOSE(UNIT=11)
	CLOSE(UNIT=12)
	CLOSE(UNIT=13)
		
	CLOSE(UNIT=21)
	CLOSE(UNIT=22)
	CLOSE(UNIT=23)
	close(UNIT=24)
	close(UNIT=25)
	close(UNIT=26)

!==============================================================================
! clean up
!==============================================================================
	
	deallocate(min_index)
	deallocate(max_index)
	
	call MPI_Finalize(ierr)

END PROGRAM generate_initial


!==============================================================================
! subroutine for integer to string conversion
!==============================================================================

subroutine int2str(value, string)
implicit none

	integer, intent(in) :: value
	character(len=*), intent(inout) :: string
	
	character(len=10) :: c
	integer :: k, n, new_value, is
	real :: e
	
	
	e=1e9
	is=0
	
	if (value==0) then
		string(1:1)='0'
		string(2:10)=' '
	else
	
		new_value=value
	
		do k=1,10
			c(k:k)=char(floor(new_value/e)+48)
		
			if ((floor(new_value/e)==0) .and. (is==0)) then
				n=k
			else
				is=1
			endif

			new_value=new_value-e*floor(new_value/e)
			e=e/10
			string(k:k)=' '
		
		enddo
		
		string(1:10-n)=c(n+1:10)
		
	endif

	if (len(string)>10) then
		string(11:len(string))=' '
	endif

return
end subroutine int2str
