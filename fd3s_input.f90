!*****************************************************************************
!***************** Subroutine to load parameters *****************************
!*****************************************************************************

SUBROUTINE fd3s_input
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'


    CHARACTER(len=15) :: junk    ! dummy variable 
    integer :: pp, dummy
    integer :: status(MPI_STATUS_SIZE)
    
    integer :: nx_min_loc, nx_max_loc
    integer :: ny_min_loc, ny_max_loc
    integer :: nz_min_loc, nz_max_loc
    
    integer :: ix_multi_loc, iy_multi_loc, iz_multi_loc
    
    real :: xmin_loc, xmax_loc
    real :: ymin_loc, ymax_loc
    real :: zmin_loc, zmax_loc


    	write(99,*) 'begin input'
    
	!=====================================================================
	! read box file, only process 0, and send to other processes
	!=====================================================================
	
	if (my_rank==0) then
		
		write(*,*) '----------------------------------------'
		write(*,*) 'begin input'
		write(*,*) 'process 0 reading boxfile' 
		write(99,*) 'process 0 reading boxfile'
		open(UNIT=15,FILE='MODELS/boxfile',STATUS='OLD',ACTION='READ')
	
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		read(15,*) junk
		
		read(15,*) pp
		write(*,*) pp
		read(15,*) px
		write(*,*) px
		read(15,*) py
		write(*,*) py
		read(15,*) pz
		write(*,*) pz
		
		if (pp/=p) then
			
			write(*,*) 'incorrect number of processes'
			write(*,*) 'inconsistency in boxfile and s_mpirun'
			write(99,*) 'incorrect number of processes'
			write(99,*) 'inconsistency in boxfile and s_mpirun'
			
		endif
		
		write(*,*) 'number of processes:', pp
		read(15,*) junk
		
		do k=1,p,1
			
			read(15,*) dummy
			read(15,*) ix_multi_loc, iy_multi_loc, iz_multi_loc
			read(15,*) nx_min_loc, nx_max_loc
			read(15,*) ny_min_loc, ny_max_loc
			read(15,*) nz_min_loc, nz_max_loc
			read(15,*) xmin_loc, xmax_loc
			read(15,*) ymin_loc, ymax_loc
			read(15,*) zmin_loc, zmax_loc
			read(15,*) junk
			
			if (k==1) then
			
				ix_multi=ix_multi_loc
				iy_multi=iy_multi_loc
				iz_multi=iz_multi_loc
				
				nx=nx_max_loc-nx_min_loc
				ny=ny_max_loc-ny_min_loc
				nz=nz_max_loc-nz_min_loc
				
				xmin=xmin_loc
				xmax=xmax_loc
				ymin=ymin_loc
				ymax=ymax_loc
				zmin=zmin_loc
				zmax=zmax_loc
				
			else
				call MPI_Send(ix_multi_loc, 1, MPI_INTEGER, k-1, 1, MPI_COMM_WORLD, ierr)
				call MPI_Send(iy_multi_loc, 1, MPI_INTEGER, k-1, 2, MPI_COMM_WORLD, ierr)
				call MPI_Send(iz_multi_loc, 1, MPI_INTEGER, k-1, 3, MPI_COMM_WORLD, ierr)
				call MPI_Send(nx_max_loc-nx_min_loc, 1, MPI_INTEGER, k-1, 4, MPI_COMM_WORLD, ierr)
				call MPI_Send(ny_max_loc-ny_min_loc, 1, MPI_INTEGER, k-1, 5, MPI_COMM_WORLD, ierr)
				call MPI_Send(nz_max_loc-nz_min_loc, 1, MPI_INTEGER, k-1, 6, MPI_COMM_WORLD, ierr)
				call MPI_Send(xmin_loc, 1, MPI_REAL, k-1, 7, MPI_COMM_WORLD, ierr)
				call MPI_Send(xmax_loc, 1, MPI_REAL, k-1, 8, MPI_COMM_WORLD, ierr)
				call MPI_Send(ymin_loc, 1, MPI_REAL, k-1, 9, MPI_COMM_WORLD, ierr)
				call MPI_Send(ymax_loc, 1, MPI_REAL, k-1, 10, MPI_COMM_WORLD, ierr)
				call MPI_Send(zmin_loc, 1, MPI_REAL, k-1, 11, MPI_COMM_WORLD, ierr)
				call MPI_Send(zmax_loc, 1, MPI_REAL, k-1, 12, MPI_COMM_WORLD, ierr)
				
			endif
			
		enddo
		
		close(UNIT=15)
		
	else
		
		call MPI_Recv(ix_multi, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(iy_multi, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(iz_multi, 1, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(nx, 1, MPI_INTEGER, 0, 4, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ny, 1, MPI_INTEGER, 0, 5, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(nz, 1, MPI_INTEGER, 0, 6, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(xmin, 1, MPI_REAL, 0, 7, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(xmax, 1, MPI_REAL, 0, 8, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ymin, 1, MPI_REAL, 0, 9, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(ymax, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(zmin, 1, MPI_REAL, 0, 11, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(zmax, 1, MPI_REAL, 0, 12, MPI_COMM_WORLD, status, ierr)
		
	endif
	
	call mpi_bcast(p,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(px,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(py,1,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(pz,1,mpi_integer,0,mpi_comm_world,ierr)
	
	write(99,*) '----------------------------------------'
	write(99,*) 'ix_multi=',ix_multi, 'iy_multi=',iy_multi, 'iz_multi=',iz_multi
	write(99,*) 'nx=',nx, 'ny=',ny, 'nz=',nz
	write(99,*) 'x_min=',xmin, 'xmax=',xmax
	write(99,*) 'ymin=',ymin, 'ymax=',ymax
	write(99,*) 'zmin=',zmin, 'zmax=',zmax
	write(99,*) '----------------------------------------'
	
	!=====================================================================
	! read Par file, only process 0
	!=====================================================================

	if (my_rank==0) then

		write(*,*) 'process 0 reading Par file'
		write(99,*) 'process 0 reading Par file'
		
		OPEN (UNIT=15,FILE='Par',STATUS='OLD',ACTION='READ') 

		READ(15,*)junk
		READ(15,1)seisfile
		!WRITE(*,*) 'seisfile=',seisfile
		READ(15,*)recfile
		!WRITE(*,1) 'recfile=',recfile
		READ(15,*)nt
		!WRITE(*,*) 'nt=', nt
		READ(15,*)dt
		!WRITE(*,*) 'dt=', dt
		!READ(15,*)izfree
		!WRITE(*,*) 'izfree=', izfree
		READ(15,*)junk
		READ(15,*)ibound
		!WRITE(*,*) 'ibound=', ibound
		READ(15,*)nabs
		!WRITE(*,*) 'nabs=', nabs
		READ(15,*)junk
		READ(15,*)xxs
		WRITE(*,*) 'xxs=',xxs
		READ(15,*)yys
		WRITE(*,*) 'yys=',yys
		READ(15,*)zzs
		WRITE(*,*) 'zzs=',zzs
		READ(15,*)source_type
		!WRITE(*,*) 'source_type=',source_type
		READ(15,*)dip
			dip=dip/180.*pi
		!WRITE(*,*) 'dip=',dip
		READ(15,*)rake
			rake=rake/180.*pi
		!WRITE(*,*) 'rake=',rake
		READ(15,*)strike
			strike=strike/180.*pi
		!WRITE(*,*) 'strike=',strike
		READ(15,*)junk
		READ(15,*)ssamp
		!WRITE(*,*) 'ssamp=',ssamp
		READ(15,*)iflush
		!WRITE(*,*) 'iflush=',iflush
		READ(15,*)iplanex
		!WRITE(*,*) 'iplanex=',iplanex
		READ(15,*)iplaney
		!WRITE(*,*) 'iplaney=',iplaney
		READ(15,*)iplanez
		!WRITE(*,*) 'iplanez=',iplanez
		READ(15,*)icheck
		!WRITE(*,*) 'icheck=',icheck
		READ(15,*)xmin_global
		!WRITE(*,*) 'xmin=',xmin
		READ(15,*)xmax_global
		!WRITE(*,*) 'xmax=',xmax
		READ(15,*)ymin_global
		!WRITE(*,*) 'ymin=',ymin
		READ(15,*)ymax_global
		!WRITE(*,*) 'ymax=',ymax
		READ(15,*)zmin_global
		!WRITE(*,*) 'zmin=',zmin
		READ(15,*)zmax_global
		!WRITE(*,*) 'zmax=',zmax
		!READ(15,*)rotdiv
		!WRITE(*,*) 'rotdiv=',rotdiv
		!READ(15,*)fs_model
		!WRITE(*,*) 'fs_model=',fs_model
		!READ(15,*)elast_inter
		!WRITE(*,*) 'elast_inter=',elast_inter
		READ(15,*)is_aniso
		READ(15,*)is_diss
		READ(15,*)is_homogeneous
		!write(*,*) is_homogeneous
		READ(15,*) is_formatted

		CLOSE(15)
		
 1  FORMAT(a,a)
 2  FORMAT(a,i8)
 3  FORMAT(a,e12.6)
 4  FORMAT(f5.0)
 5  FORMAT(i5)

		
	endif
	
	!=====================================================================
	! broadcast input parameters to the other processes
	!=====================================================================
	
	call mpi_bcast(seisfile, len(seisfile), mpi_character, 0, mpi_comm_world, ierr) 
	call mpi_bcast(recfile, len(recfile), mpi_character, 0, mpi_comm_world, ierr)
	call mpi_bcast(nt, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(dt, 1, mpi_real, 0, mpi_comm_world, ierr)
	!call mpi_bcast(izfree, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(ibound, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(nabs, 1, mpi_integer, 0, mpi_comm_world, ierr)
	
	call mpi_bcast(xxs, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(yys, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zzs, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(source_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(dip, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(rake, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(strike, 1, mpi_real, 0, mpi_comm_world, ierr)
	
	call mpi_bcast(ssamp, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(iflush, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(iplanex, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(iplaney, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(iplanez, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(icheck, 1, mpi_integer, 0, mpi_comm_world, ierr)
	
	call mpi_bcast(xmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(xmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(ymax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmin_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(zmax_global, 1, mpi_real, 0, mpi_comm_world, ierr)
	
	!call mpi_bcast(rotdiv, 1, mpi_integer, 0, mpi_comm_world, ierr)
	!call mpi_bcast(elast_inter, 1, mpi_integer, 0, mpi_comm_world, ierr)
	!call mpi_bcast(fs_model, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(is_aniso, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(is_diss, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(is_homogeneous, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(is_formatted, 1, mpi_integer, 0, mpi_comm_world, ierr)
	

	!=============================================================================
	! read source time function, adequate process only
	!=============================================================================

	if ((xmin .LE. xxs) .AND. (xmax .GE. xxs) .AND. (ymin .LE. yys) .AND. (ymax .GE. yys) .AND. &
		(zmin .LE. zmax_global-zzs) .AND. (zmax .GE. zmax_global-zzs)) then
		
		write(99,*) 'process ', my_rank, 'reading source time function'
		is_source=1.0
		
		OPEN(UNIT=15,FILE='MODELS/stf',STATUS='OLD',ACTION='READ')
	
		DO i=1, nt, 1
			READ(15,*) so(i)
		END DO
	
		CLOSE(UNIT=15)
    
	else
		
		write(99,*) 'source not located in this processor box'
		is_source=0.0
		
	endif
	
	call mpi_reduce(is_source*my_rank,source_processor,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
	
	if (my_rank==0) then
		
		write(*,*) 'process', int(source_processor), 'reading source time function'
	
	endif
	
	!======================================================================
	! read initial stress and velocity values, formatted or unformatted
	!======================================================================
	
	call int2str(my_rank,junk)
	
	if (my_rank==0) then
		write(*,*) "reading initial values"
	endif
	
	write(99,*) 'reading initial values'
	
	if (is_homogeneous==1) then
		
		if (my_rank==0) then
			write(*,*) "homogeneous"
		endif
		write(99,*) "homogeneous"
		
		w11(:,:,:)=0.0
		w12(:,:,:)=0.0
		w13(:,:,:)=0.0
		w14(:,:,:)=0.0
		w15(:,:,:)=0.0
		w16(:,:,:)=0.0
		w17(:,:,:)=0.0
		w18(:,:,:)=0.0
		w19(:,:,:)=0.0
		
	else
		
		if (my_rank==0) then
			write(*,*) "non-homogeneous"
		endif
		write(*,*) "non-homogeneous"
	
		if (is_formatted==1) then
		
			OPEN (UNIT=11,FILE='MODELS/v0_theta'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=12,FILE='MODELS/v0_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=13,FILE='MODELS/v0_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
	
			OPEN (UNIT=14,FILE='MODELS/s0_theta_theta'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=15,FILE='MODELS/s0_phi_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=16,FILE='MODELS/s0_r_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=17,FILE='MODELS/s0_theta_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=18,FILE='MODELS/s0_theta_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=19,FILE='MODELS/s0_phi_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			
		else
			
			OPEN (UNIT=11,FILE='MODELS/v0_theta'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=12,FILE='MODELS/v0_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=13,FILE='MODELS/v0_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
	
			OPEN (UNIT=14,FILE='MODELS/s0_theta_theta'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=15,FILE='MODELS/s0_phi_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=16,FILE='MODELS/s0_r_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=17,FILE='MODELS/s0_theta_phi'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=18,FILE='MODELS/s0_theta_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=19,FILE='MODELS/s0_phi_r'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			
		endif
	
		w11(:,:,:)=0
		w12(:,:,:)=0
		w13(:,:,:)=0
		w14(:,:,:)=0
		w15(:,:,:)=0
		w16(:,:,:)=0
		w17(:,:,:)=0
		w18(:,:,:)=0
		w19(:,:,:)=0
	
		if (is_formatted==1) then
		
			DO index_z=0, nz, 1
				DO index_y=0, ny, 1
					DO index_x=0, nx, 1
				
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! initial velocities
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
						READ(11,*) w11(index_x, index_y, index_z)
						READ(12,*) w12(index_x, index_y, index_z)
						READ(13,*) w13(index_x, index_y, index_z)
					
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! initial stresses
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
						READ(14,*) w14(index_x, index_y, index_z)
						READ(15,*) w15(index_x, index_y, index_z)
						READ(16,*) w16(index_x, index_y, index_z)
						READ(17,*) w17(index_x, index_y, index_z)
						READ(18,*) w18(index_x, index_y, index_z)
						READ(19,*) w19(index_x, index_y, index_z)
					
					END DO
				END DO
			END DO
			
		else
			
			read(11) w11(0:nx,0:ny,0:nz)
			read(12) w12(0:nx,0:ny,0:nz)
			read(13) w13(0:nx,0:ny,0:nz)
		
			read(14) w14(0:nx,0:ny,0:nz)
			read(15) w15(0:nx,0:ny,0:nz)
			read(16) w16(0:nx,0:ny,0:nz)
			read(17) w17(0:nx,0:ny,0:nz)
			read(18) w18(0:nx,0:ny,0:nz)
			read(19) w19(0:nx,0:ny,0:nz)
			
		endif
				
		CLOSE(11)
		CLOSE(12)
		CLOSE(13)
	
		CLOSE(14)
		CLOSE(15)
		CLOSE(16)
		CLOSE(17)
		CLOSE(18)
		CLOSE(19)
		
	endif
  
	!=============================================================================
	! Read structural information (rho, mu, lambda, A, B, C), formatted or unformatted
	!=============================================================================

	rhoinv(:,:,:)=1.0
	lam(:,:,:)=1.0
	mu(:,:,:)=1.0

	if (my_rank==0) then
		write(*,*) "reading structural information"
	endif
	
	write(99,*) 'reading structural information'
	
	if (my_rank==0) then
		
		IF (is_diss==1) THEN
			WRITE(*,*) 'Dissipation on'
			WRITE(99,*) 'Dissipation on'
		ELSE
			WRITE(*,*) 'Dissipation off'
			WRITE(99,*) 'Dissipation off'
		END IF
	
		IF (is_aniso==1) THEN
			WRITE(*,*) 'Anisotropy on'
			WRITE(99,*) 'Anisotropy on'
		ELSE
			WRITE(*,*) 'Anisotropy off'
			WRITE(99,*) 'Anisotropy off'
		END IF
		
	endif
	
	if (is_aniso==1) then
	
		A(:,:,:)=1.0
		B(:,:,:)=1.0
		C(:,:,:)=1.0
		
		if (is_formatted==1) then
		
			OPEN (UNIT=21,FILE='MODELS/A'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=22,FILE='MODELS/B'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=23,FILE='MODELS/C'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			
		else
			
			OPEN (UNIT=21,FILE='MODELS/A'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=22,FILE='MODELS/B'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=23,FILE='MODELS/C'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			
		endif
		
	endif
	
	if (is_diss==1) then
		
		tau_lambda(:,:,:)=1.0
		tau_mu(:,:,:)=1.0
		tau_s1_lambda(:,:,:)=1.0
		tau_s2_lambda(:,:,:)=1.0
		tau_s1_mu(:,:,:)=1.0
		tau_s2_mu(:,:,:)=1.0
		
		if (is_formatted==1) then
		
			OPEN (UNIT=31,FILE='MODELS/tau_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=32,FILE='MODELS/tau_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
	
			OPEN (UNIT=41,FILE='MODELS/tau_s1_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=42,FILE='MODELS/tau_s2_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=43,FILE='MODELS/tau_s1_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			OPEN (UNIT=44,FILE='MODELS/tau_s2_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
			
		else
			
			OPEN (UNIT=31,FILE='MODELS/tau_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=32,FILE='MODELS/tau_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
	
			OPEN (UNIT=41,FILE='MODELS/tau_s1_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=42,FILE='MODELS/tau_s2_lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=43,FILE='MODELS/tau_s1_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			OPEN (UNIT=44,FILE='MODELS/tau_s2_mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
			
		endif
		
	endif
	
	if (is_formatted==1) then
	
		OPEN (UNIT=15,FILE='MODELS/rhoinv'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
		OPEN (UNIT=16,FILE='MODELS/mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
		OPEN (UNIT=17,FILE='MODELS/lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ')
		
	else
		
		OPEN (UNIT=15,FILE='MODELS/rhoinv'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
		OPEN (UNIT=16,FILE='MODELS/mu'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
		OPEN (UNIT=17,FILE='MODELS/lambda'//junk(1:len_trim(junk)),STATUS='OLD',ACTION='READ',form='unformatted')
		
	endif

	if (is_formatted==1) then
	
		DO index_z=0, nz, 1
			DO index_y=0, ny, 1
				DO index_x=0, nx, 1
			
					READ(15,*) rhoinv(index_x, index_y, index_z)
					READ(16,*) mu(index_x, index_y, index_z)
					READ(17,*) lam(index_x, index_y, index_z)
				
					IF (is_aniso==1) THEN
						READ(21,*) A(index_x,index_y,index_z)
						READ(22,*) B(index_x,index_y,index_z)
						READ(23,*) C(index_x,index_y,index_z)
					END IF
					
					IF (is_diss==1) THEN
						READ(31,*) tau_lambda(index_x,index_y,index_z)
						READ(32,*) tau_mu(index_x,index_y,index_z)
				
						READ(41,*) tau_s1_lambda(index_x,index_y,index_z)
						READ(42,*) tau_s2_lambda(index_x,index_y,index_z)
						READ(43,*) tau_s1_mu(index_x,index_y,index_z)
						READ(44,*) tau_s2_mu(index_x,index_y,index_z)
					END IF
				END DO
			END DO
		END DO
		
	else
		
		read(15) rhoinv(0:nx,0:ny,0:nz)
		read(16) mu(0:nx,0:ny,0:nz)
		read(17) lam(0:nx,0:ny,0:nz)
		
		if (is_aniso==1) then
			
			read(21) A(0:nx,0:ny,0:nz)
			read(22) B(0:nx,0:ny,0:nz)
			read(23) C(0:nx,0:ny,0:nz)
			
		endif
		
		if (is_diss==1) then
			
			read(31) tau_lambda(0:nx,0:ny,0:nz)
			read(32) tau_mu(0:nx,0:ny,0:nz)
			
			read(41) tau_s1_lambda(0:nx,0:ny,0:nz)
			read(42) tau_s2_lambda(0:nx,0:ny,0:nz)
			read(43) tau_s1_mu(0:nx,0:ny,0:nz)
			read(44) tau_s2_mu(0:nx,0:ny,0:nz)
			
		endif
		
	endif
				
	CLOSE(15)
	CLOSE(16)
	CLOSE(17)

	
	if (is_aniso==1) then
	
		CLOSE(21)
		CLOSE(22)
		CLOSE(23)
		
	endif
	
	if (is_diss==1) then
	
		CLOSE(31)
		CLOSE(32)
	
		CLOSE(41)
		CLOSE(42)
		CLOSE(43)
		CLOSE(44)
		
	endif

	!======================================================================
	! read receiver locations, process  0 only
	!======================================================================
	
	if (my_rank==0) then
		
		write(*,*) 'process 0 reading receiver locations'
		write(99,*) 'process 0 reading receiver locations'
		
		open(unit=10, file=recfile)
		
		read(10,*) nr_global	!number of receivers
		
		do i=1,nr_global
			read(10,*) recloc_global(1,i), recloc_global(2,i)
			!write(*,*) recloc_global(1,i), recloc_global(2,i)
		enddo
		
		close(10)
		
	endif
	
	call mpi_bcast(recloc_global, 2*maxnr, mpi_real, 0, mpi_comm_world, ierr)
	call mpi_bcast(nr_global,1,mpi_integer,0,mpi_comm_world,ierr)
	
	write(99,*) 'end input'
	write(99,*) '----------------------------------------'
		
	return
        
END SUBROUTINE fd3s_input


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


