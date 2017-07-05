!==============================================================================
! communicate and taper boundaries of velocity field
!==============================================================================

subroutine comm_tap_v
use parameters
use variables
implicit none
include 'mpif.h'

	integer :: mi_rec		! index of the receiver, rank_receiver=mi_receiver-1
	integer :: send_tag, rec_tag, num
        integer :: status(MPI_STATUS_SIZE)
	
!------------------------------------------------------------------------------

	!======================================================================
	! ix_multi=1 -> taper left theta boundary, communicate right boundary
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(ny+1)*(nz+1)
	
	if (ix_multi==1) then
		
		! tapering of left boundary
	
		w11(0:nabs,:,:)=w11(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w12(0:nabs,:,:)=w12(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w13(0:nabs,:,:)=w13(0:nabs,:,:)*window_x_left(0:nabs,:,:)
	
		w11(-bndry:-1,:,:)=0.0
		w12(-bndry:-1,:,:)=0.0
		w13(-bndry:-1,:,:)=0.0
	
		!write(*,*) 'theta left boundary tapered'
		
		! send right boundary to process with ix_multi+1, iy_multi, iz_multi
		! and receive left boundary from ix_multi+1, iy_multi, iz_multi
		! only if a right neighbor exists
	
		if (px > 1) then
			
			mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			
			call mpi_send(w11((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
			call mpi_send(w13((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
			call mpi_recv(w11((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
		endif

	endif
	
	!======================================================================
	! ix_multi=px -> taper right theta boundary, communicate left boundary
	!======================================================================
	
	!call write_cpu_time
	
	if (ix_multi==px) then
		
		! tapering right boundary
	
		w11((nx-nabs):nx,:,:)=w11((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w12((nx-nabs):nx,:,:)=w12((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w13((nx-nabs):nx,:,:)=w13((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		
		w11((nx+1):(nx+bndry),:,:)=0.0
		w12((nx+1):(nx+bndry),:,:)=0.0
		w13((nx+1):(nx+bndry),:,:)=0.0
		
		!write(99,*) 'theta right boundary tapered'
		
		! send the left boundary to process with ix_multi-1, iy_multi, iz_multi
		! and receive right boundary from the same process
		! only if a left neighbour exists
		
		if (px > 1) then
			
			mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w11(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
			call mpi_send(w13(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
		
			rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
			call mpi_recv(w11((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'theta boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! ix_multi /= px AND /= 1 -> communicate both boundaries
	!======================================================================
	
	!call write_cpu_time
	
	if ((ix_multi < px) .and. (ix_multi > 1)) then
		
		! send right boundary to right neighbor, receive left boundary from right neighbor 
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
	
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w11((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
		call mpi_send(w13((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
		call mpi_recv(w11((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'theta boundary exchange with process', mi_rec-1
		
		! send left boundary to left neighbor, receive right boundary from left neighbor
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w11(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
		call mpi_send(w13(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		call mpi_recv(w11((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'theta boundary exchange with process', mi_rec-1
		
	endif
	
!------------------------------------------------------------------------------
	
	!======================================================================
	! iy_multi=1 -> taper left phi boundary, communicate right boundary
	if (my_rank==0) then
		write(*,*) '.'
	endif
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(nx+1)*(nz+1)
	
	if (iy_multi==1) then
	
		! tapering of left boundary
	
		w11(:,0:nabs,:)=w11(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w12(:,0:nabs,:)=w12(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w13(:,0:nabs,:)=w13(:,0:nabs,:)*window_y_left(:,0:nabs,:)
	
		w11(:,-bndry:-1,:)=0.0
		w12(:,-bndry:-1,:)=0.0
		w13(:,-bndry:-1,:)=0.0
		
		!write(99,*) 'phi left boundary tapered'
	
		! send right boundary to process with ix_multi, iy_multi+1, iz_multi
		! and receive left boundary from ix_multi, iy_multi+1, iz_multi
		! only if a right neighbor exists
	
		if (py > 1) then
		
			mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
	
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w11(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
			call mpi_send(w13(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
			call mpi_recv(w11(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'phi boundary exchange with process', mi_rec-1
			
		endif
		
		! AN DIESER STELLE MUESSTE MAN PRUEFEN, OB DIE FELDER RICHTIG KOPIERT WERDEN
		
	endif
	
	!======================================================================
	! iy_multi=py -> taper right phi boundary, communicate left boundary
	!======================================================================
	
	if (iy_multi==py) then
		
		! tapering right boundary
	
		w11(:,(ny-nabs):ny,:)=w11(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w12(:,(ny-nabs):ny,:)=w12(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w13(:,(ny-nabs):ny,:)=w13(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		
		w11(:,(ny+1):(ny+bndry),:)=0.0
		w12(:,(ny+1):(ny+bndry),:)=0.0
		w13(:,(ny+1):(ny+bndry),:)=0.0
		
		!write(99,*) 'phi right boundary tapered'
		
		! send the left boundary to process with ix_multi, iy_multi-1, iz_multi
		! and receive right boundary from the same process
		! only if a left neighbour exists
		
		if (py > 1) then
			
			mi_rec=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w11(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
			call mpi_send(w13(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
			call mpi_recv(w11(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'phi boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! iy_multi /= py AND /= 1 -> communicate both boundaries
	!======================================================================
	
	if ((iy_multi < py) .and. (iy_multi > 1)) then
		
		! send right boundary to right neighbor, receive left boundary from right neighbor 
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
	
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w11(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
		call mpi_send(w13(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
		call mpi_recv(w11(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
	
		!write(99,*) 'phi boundary exchange with process', mi_rec-1
		
		! send left boundary to left neighbor, receive right boundary from left neighbor
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w11(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
		call mpi_send(w13(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		call mpi_recv(w11(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'phi boundary exchange with process', mi_rec-1
		
	endif
		
!------------------------------------------------------------------------------

	!======================================================================
	! iz_multi=1 -> taper lower z boundary, communicate upper boundary
	if (my_rank==0) then
		write(*,*) '.'
	endif
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(nx+1)*(ny+1)
	
	if (iz_multi==1) then
	
		! tapering of lower boundary (large z-indeces)
	
		w11(:,:,(nz-nabs):nz)=w11(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w12(:,:,(nz-nabs):nz)=w12(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w13(:,:,(nz-nabs):nz)=w13(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
	
		w11(:,:,(nz+1):(nz+bndry))=0.0
		w12(:,:,(nz+1):(nz+bndry))=0.0
		w13(:,:,(nz+1):(nz+bndry))=0.0
		
		!write(99,*) 'z lower boundary tapered'
	
		! send upper boundary to process with ix_multi, iy_multi, iz_multi+1
		! and receive lower boundary from ix_multi, iy_multi, iz_multi+1
		! only if a upper neighbor exists
	
		if (pz > 1) then
		
			!write(99,*) 'starting boundary exchange iz_multi=1'
			
			mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
			
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi !=my_rank+1
		        call mpi_send(w11(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+p, mpi_comm_world,ierr)
			call mpi_send(w13(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+p+p, mpi_comm_world,ierr)
			
			rec_tag=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
		        call mpi_recv(w11(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
			write(99,*) 'z boundary exchange with process', mi_rec-1
			
		endif
		
		! AN DIESER STELLE MUESSTE MAN PRUEFEN, OB DIE FELDER RICHTIG KOPIERT WERDEN
		
	endif
	
	!======================================================================
	! iz_multi=pz -> leave upper boundary unchanged, communicate lower boundary
	!======================================================================
	
	if (iz_multi==pz) then
		
		! send the lower boundary to process with ix_multi, iy_multi, iz_multi-1
		! and receive upper boundary from the same process
		! only if a lower neighbour exists
		
		if (pz > 1) then
			
			!write(99,*) 'starting boundary exchange iz_multi=pz'
			
			mi_rec=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi 
			call mpi_send(w11(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
			call mpi_send(w12(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
			call mpi_send(w13(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_recv(w11(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
			call mpi_recv(w12(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
			call mpi_recv(w13(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
			
			write(99,*) 'z boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! iz_multi /= pz AND /= 1 -> communicate both boundaries
	!======================================================================
	
	if ((iz_multi < pz) .and. (iz_multi > 1)) then
		
		!write(99,*) 'starting boundary exchange pz>iz_multi>1'
		
		! send upper boundary to upper neighbor, receive lower boundary from upper neighbor 
		
		mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
			
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi !=my_rank+1
		call mpi_send(w11(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+p, mpi_comm_world,ierr)
		call mpi_send(w13(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+p+p, mpi_comm_world,ierr)
			
		rec_tag=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_recv(w11(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)
		
		write(99,*) 'z boundary exchange with process', mi_rec-1
		
		! send lower boundary to lower neighbor, receive upper boundary from lower neighbor
		
		mi_rec=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi 
		call mpi_send(w11(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag,mpi_comm_world,ierr)
		call mpi_send(w12(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+p,mpi_comm_world,ierr)
		call mpi_send(w13(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+p+p,mpi_comm_world,ierr)
			
		rec_tag=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_recv(w11(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag,mpi_comm_world,status,ierr)
		call mpi_recv(w12(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+p,mpi_comm_world,status,ierr)
		call mpi_recv(w13(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+p+p,mpi_comm_world,status,ierr)	
		
		write(99,*) 'z boundary exchange with process', mi_rec-1
		
	endif
		
end subroutine comm_tap_v


!==============================================================================
! communicate and taper boundaries of stress field
!==============================================================================

subroutine comm_tap_s
use parameters
use variables
implicit none
include 'mpif.h'

	integer :: mi_rec		! index of the receiver, rank_receiver=mi_receiver-1
	integer :: send_tag, rec_tag, num
	integer :: status(MPI_STATUS_SIZE)
	!real :: D=0.4		! damping parameter
	
	!integer :: time_dummy1, time_dummy2
	
	
!------------------------------------------------------------------------------

	!======================================================================
	! ix_multi=1 -> taper left theta boundary, communicate right boundary
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(ny+1)*(nz+1)
	
	if (ix_multi==1) then
	
		! tapering of left boundary
	
		w14(0:nabs,:,:)=w14(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w15(0:nabs,:,:)=w15(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w16(0:nabs,:,:)=w16(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w17(0:nabs,:,:)=w17(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w18(0:nabs,:,:)=w18(0:nabs,:,:)*window_x_left(0:nabs,:,:)
		w19(0:nabs,:,:)=w19(0:nabs,:,:)*window_x_left(0:nabs,:,:)
	
		w14(-bndry:-1,:,:)=0.0
		w15(-bndry:-1,:,:)=0.0
		w16(-bndry:-1,:,:)=0.0
		w17(-bndry:-1,:,:)=0.0
		w18(-bndry:-1,:,:)=0.0
		w19(-bndry:-1,:,:)=0.0
	
		!write(99,*) 'theta left boundary tapered'
		
		! send right boundary to process with ix_multi+1, iy_multi, iz_multi
		! and receive left boundary from ix_multi+1, iy_multi, iz_multi
		! only if a right neighbor exists
	
		if (px > 1) then
		
			mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
	
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w14((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
			call mpi_recv(w14((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'theta boundary exchange with process', mi_rec-1
			
		endif
		
		! AN DIESER STELLE MUESSTE MAN PRUEFEN, OB DIE FELDER RICHTIG KOPIERT WERDEN
		
	endif
	
	!======================================================================
	! ix_multi=px -> taper right theta boundary, communicate left boundary
	!======================================================================
	
	if (ix_multi==px) then
		
		! tapering right boundary
	
		w14((nx-nabs):nx,:,:)=w14((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w15((nx-nabs):nx,:,:)=w15((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w16((nx-nabs):nx,:,:)=w16((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w17((nx-nabs):nx,:,:)=w17((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w18((nx-nabs):nx,:,:)=w18((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		w19((nx-nabs):nx,:,:)=w19((nx-nabs):nx,:,:)*window_x_right(0:nabs,:,:)
		
		w14((nx+1):(nx+bndry),:,:)=0.0
		w15((nx+1):(nx+bndry),:,:)=0.0
		w16((nx+1):(nx+bndry),:,:)=0.0
		w17((nx+1):(nx+bndry),:,:)=0.0
		w18((nx+1):(nx+bndry),:,:)=0.0
		w19((nx+1):(nx+bndry),:,:)=0.0
		
		!write(99,*) 'theta right boundary tapered'
		
		! send the left boundary to process with ix_multi-1, iy_multi, iz_multi
		! and receive right boundary from the same process
		! only if a left neighbour exists
		
		if (px > 1) then
			
			mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w14(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
		
			rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
			call mpi_recv(w14((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'theta boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! ix_multi /= px AND /= 1 -> communicate both boundaries
	!======================================================================
	
	if ((ix_multi < px) .and. (ix_multi > 1)) then
		
		! send right boundary to right neighbor, receive left boundary from right neighbor 
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
	
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w14((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19((nx-bndry+1):nx,0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi+1
		call mpi_recv(w14((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19((nx+1):(nx+bndry),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'theta boundary exchange with process', mi_rec-1
		
		! send left boundary to left neighbor, receive right boundary from left neighbor
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w14(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19(0:(bndry-1),0:ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi-1
		call mpi_recv(w14((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19((-bndry):(-1),0:ny,0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'theta boundary exchange with process', mi_rec-1
		
	endif
	
!------------------------------------------------------------------------------
	
	!======================================================================
	! iy_multi=1 -> taper left phi boundary, communicate right boundary
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(nx+1)*(nz+1)
	
	if (iy_multi==1) then
	
		! tapering of left boundary
	
		w14(:,0:nabs,:)=w14(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w15(:,0:nabs,:)=w15(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w16(:,0:nabs,:)=w16(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w17(:,0:nabs,:)=w17(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w18(:,0:nabs,:)=w18(:,0:nabs,:)*window_y_left(:,0:nabs,:)
		w19(:,0:nabs,:)=w19(:,0:nabs,:)*window_y_left(:,0:nabs,:)
	
		w14(:,-bndry:-1,:)=0.0
		w15(:,-bndry:-1,:)=0.0
		w16(:,-bndry:-1,:)=0.0
		w17(:,-bndry:-1,:)=0.0
		w18(:,-bndry:-1,:)=0.0
		w19(:,-bndry:-1,:)=0.0
		
		!write(99,*) 'phi left boundary tapered'
	
		! send right boundary to process with ix_multi, iy_multi+1, iz_multi
		! and receive left boundary from ix_multi, iy_multi+1, iz_multi
		! only if a right neighbor exists
	
		if (py > 1) then
		
			mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
	
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w14(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
			call mpi_recv(w14(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'phi boundary exchange with process', mi_rec-1
			
		endif
		
		! AN DIESER STELLE MUESSTE MAN PRUEFEN, OB DIE FELDER RICHTIG KOPIERT WERDEN
		
	endif
	
	!======================================================================
	! iy_multi=py -> taper right phi boundary, communicate left boundary
	!======================================================================
	
	if (iy_multi==py) then
		
		! tapering right boundary
	
		w14(:,(ny-nabs):ny,:)=w14(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w15(:,(ny-nabs):ny,:)=w15(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w16(:,(ny-nabs):ny,:)=w16(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w17(:,(ny-nabs):ny,:)=w17(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w18(:,(ny-nabs):ny,:)=w18(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
		w19(:,(ny-nabs):ny,:)=w19(:,(ny-nabs):ny,:)*window_y_right(:,0:nabs,:)
	
	
		!do i=0,nabs
		!
		!	w14(:,ny-nabs+i,:)=w14(:,ny-nabs+i,:)*tap1(nabs-i)
		!	w15(:,ny-nabs+i,:)=w15(:,ny-nabs+i,:)*tap1(nabs-i)
		!	w16(:,ny-nabs+i,:)=w16(:,ny-nabs+i,:)*tap1(nabs-i)
		!	w17(:,ny-nabs+i,:)=w17(:,ny-nabs+i,:)*tap1(nabs-i)
		!	w18(:,ny-nabs+i,:)=w18(:,ny-nabs+i,:)*tap1(nabs-i)
		!	w19(:,ny-nabs+i,:)=w19(:,ny-nabs+i,:)*tap1(nabs-i)
		!
		!enddo
		
		w14(:,(ny+1):(ny+bndry),:)=0.0
		w15(:,(ny+1):(ny+bndry),:)=0.0
		w16(:,(ny+1):(ny+bndry),:)=0.0
		w17(:,(ny+1):(ny+bndry),:)=0.0
		w18(:,(ny+1):(ny+bndry),:)=0.0
		w19(:,(ny+1):(ny+bndry),:)=0.0
		
		!write(99,*) 'phi right boundary tapered'
		
		! send the left boundary to process with ix_multi, iy_multi-1, iz_multi
		! and receive right boundary from the same process
		! only if a left neighbour exists
		
		if (py > 1) then
			
			mi_rec=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_send(w14(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
			call mpi_recv(w14(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			!write(99,*) 'phi boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! iy_multi /= py AND /= 1 -> communicate both boundaries
	!======================================================================
	
	if ((iy_multi < py) .and. (iy_multi > 1)) then
		
		! send right boundary to right neighbor, receive left boundary from right neighbor 
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
	
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w14(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19(0:nx,(ny-bndry+1):ny,0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi)*px+ix_multi
		call mpi_recv(w14(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19(0:nx,(ny+1):(ny+bndry),0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
	
		!write(99,*) 'phi boundary exchange with process', mi_rec-1
		
		! send left boundary to left neighbor, receive right boundary from left neighbor
		
		mi_rec=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_send(w14(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19(0:nx,0:(bndry-1),0:nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
		
		rec_tag=(iz_multi-1)*py*px+(iy_multi-2)*px+ix_multi
		call mpi_recv(w14(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19(0:nx,(-bndry):(-1),0:nz),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
		
		!write(99,*) 'phi boundary exchange with process', mi_rec-1
		
	endif
		
!------------------------------------------------------------------------------

	!======================================================================
	! iz_multi=1 -> taper lower z boundary, communicate upper boundary
	!======================================================================
	
	!call write_cpu_time
	
	num=bndry*(nx+1)*(ny+1)
	
	if (iz_multi==1) then
	
		! tapering of lower boundary (large z-indeces)
	
		w14(:,:,(nz-nabs):nz)=w14(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w15(:,:,(nz-nabs):nz)=w15(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w16(:,:,(nz-nabs):nz)=w16(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w17(:,:,(nz-nabs):nz)=w17(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w18(:,:,(nz-nabs):nz)=w18(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
		w19(:,:,(nz-nabs):nz)=w19(:,:,(nz-nabs):nz)*window_z(:,:,0:nabs)
	
		w14(:,:,(nz+1):(nz+bndry))=0.0
		w15(:,:,(nz+1):(nz+bndry))=0.0
		w16(:,:,(nz+1):(nz+bndry))=0.0
		w17(:,:,(nz+1):(nz+bndry))=0.0
		w18(:,:,(nz+1):(nz+bndry))=0.0
		w19(:,:,(nz+1):(nz+bndry))=0.0
		
		!write(99,*) 'z lower boundary tapered'
	
		! send upper boundary to process with ix_multi, iy_multi, iz_multi+1
		! and receive lower boundary from ix_multi, iy_multi, iz_multi+1
		! only if a upper neighbor exists
	
		if (pz > 1) then
		
			!write(99,*) 'starting boundary exchange iz_multi=1'
			
			mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
			
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi !=my_rank+1
		        call mpi_send(w14(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
		        call mpi_recv(w14(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			write(99,*) 'z boundary exchange with process', mi_rec-1
			
		endif
		
		! AN DIESER STELLE MUESSTE MAN PRUEFEN, OB DIE FELDER RICHTIG KOPIERT WERDEN
		
	endif
	
	!======================================================================
	! iz_multi=pz -> leave upper boundary unchanged, communicate lower boundary
	!======================================================================
	
	if (iz_multi==pz) then
		
		! send the lower boundary to process with ix_multi, iy_multi, iz_multi-1
		! and receive upper boundary from the same process
		! only if a lower neighbour exists
		
		if (pz > 1) then
			
			!write(99,*) 'starting boundary exchange iz_multi=pz'
			
			mi_rec=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			
			send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi 
			call mpi_send(w14(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
			call mpi_send(w15(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
			call mpi_send(w16(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
			call mpi_send(w17(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
			call mpi_send(w18(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
			call mpi_send(w19(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
			rec_tag=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			call mpi_recv(w14(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
			call mpi_recv(w15(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
			call mpi_recv(w16(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
			call mpi_recv(w17(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
			call mpi_recv(w18(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
			call mpi_recv(w19(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
			
			write(99,*) 'z boundary exchange with process', mi_rec-1
			
		endif
		
	endif
	
	!======================================================================
	! iz_multi /= pz AND /= 1 -> communicate both boundaries
	!======================================================================
	
	if ((iz_multi < pz) .and. (iz_multi > 1)) then
		
		!write(99,*) 'starting boundary exchange pz>iz_multi>1'
		
		! send upper boundary to upper neighbor, receive lower boundary from upper neighbor 
		
		mi_rec=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
			
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi !=my_rank+1
		call mpi_send(w14(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19(0:nx,0:ny,0:(bndry-1)),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
		rec_tag=(iz_multi)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_recv(w14(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19(0:nx,0:ny,(-bndry):(-1)),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
		
		write(99,*) 'z boundary exchange with process', mi_rec-1
		
		! send lower boundary to lower neighbor, receive upper boundary from lower neighbor
		
		mi_rec=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
			
		send_tag=(iz_multi-1)*py*px+(iy_multi-1)*px+ix_multi 
		call mpi_send(w14(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+4*p,mpi_comm_world,ierr)
		call mpi_send(w15(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+5*p,mpi_comm_world,ierr)
		call mpi_send(w16(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+6*p,mpi_comm_world,ierr)
		call mpi_send(w17(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+7*p,mpi_comm_world,ierr)
		call mpi_send(w18(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+8*p,mpi_comm_world,ierr)
		call mpi_send(w19(0:nx,0:ny,(nz-bndry+1):nz),num,mpi_real,mi_rec-1,send_tag+9*p,mpi_comm_world,ierr)
			
		rec_tag=(iz_multi-2)*py*px+(iy_multi-1)*px+ix_multi
		call mpi_recv(w14(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+4*p,mpi_comm_world,status,ierr)
		call mpi_recv(w15(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+5*p,mpi_comm_world,status,ierr)
		call mpi_recv(w16(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+6*p,mpi_comm_world,status,ierr)
		call mpi_recv(w17(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+7*p,mpi_comm_world,status,ierr)
		call mpi_recv(w18(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+8*p,mpi_comm_world,status,ierr)
		call mpi_recv(w19(0:nx,0:ny,(nz+1):(nz+bndry)),num,mpi_real,mi_rec-1,rec_tag+9*p,mpi_comm_world,status,ierr)
		
		write(99,*) 'z boundary exchange with process', mi_rec-1
		
	endif
		
end subroutine comm_tap_s
