!=============================================================================
! Andreas Fichtner, 26. Mai 2006
!=============================================================================

PROGRAM generate_models
IMPLICIT NONE
include 'mpif.h'

!=============================================================================
! variables 
!=============================================================================

	character(len=20) :: junk
	INTEGER :: i, j, k, mi, is_aniso, is_formatted
	INTEGER :: model_type, status_rho, status_mu, status_lambda, status_A, status_B, status_C
	INTEGER :: nx, ny, nz, px, py, pz
	integer :: bnx, bny, bnz
	REAL :: rho_homogen, mu_homogen, lambda_homogen, A_homogen, B_homogen, C_homogen
	REAL :: z_max, z_min, y_max, y_min, x_max, x_min, dz, dy, dx
	REAL, ALLOCATABLE, DIMENSION(:) :: z, x, rho_1D, vp_1D, vs_1D
	REAL, ALLOCATABLE, DIMENSION(:,:,:) :: rhoinv, mu, lambda, A, B, C
	
	integer, allocatable, dimension(:,:) :: min_index, max_index		! minimum and maximum index boundaries
	integer :: status, ierr, p, my_rank


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
		read(10,*) mu_homogen
		!WRITE(*,*) mu_homogen
		read(10,*) lambda_homogen
		!WRITE(*,*) lambda_homogen
		read(10,*) rho_homogen
		!WRITE(*,*) rho_homogen
		read(10,*) A_homogen
		!WRITE(*,*) A_homogen
		read(10,*) B_homogen
		!WRITE(*,*) B_homogen
		read(10,*) C_homogen
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
	
	call MPI_BCAST(model_type,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mu_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lambda_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rho_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(A_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(B_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(C_homogen,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(is_aniso,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(is_formatted,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
	call MPI_BCAST(px,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(py,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	
!==============================================================================
! make boxfile containing the parameters for parallel computing
!==============================================================================

	if (my_rank==0) then

		write(*,*) 'writing boxfile'
	
	endif

	allocate(min_index(1:px*py*pz,1:3),STAT=status)
	allocate(max_index(1:px*py*pz,1:3),STAT=status)
	
	dz=(z_max-z_min)/nz
	dy=(y_max-y_min)/ny
	dx=(x_max-x_min)/nx
	
	bnx=floor(real((nx+1)/px))
	bny=floor(real((ny+1)/py))
	bnz=floor(real((nz+1)/pz))

	if (my_rank==0) then
	
		open(UNIT=11,FILE='boxfile',ACTION='WRITE')
		
		write(11,*) "- format description -"
		write(11,*) "total number of processes"
		write(11,*) "number of processes in x direction"
		write(11,*) "number of processes in y direction"
		write(11,*) "number of processes in z direction"
		write(11,*) "single index (npz-1)*py*px+(npy-1)*px+npx"
		write(11,*) "multi index (npx, npy, npz)"
		write(11,*) "index boundaries x (theta_min, theta_max)"
		write(11,*) "index boundaries y (phi_min, phi_max)"
		write(11,*) "index boundaries z (z_min, z_max)"
		write(11,*) "physical boundaries theta"
		write(11,*) "physical boundaries phi"
		write(11,*) "phyiscal boundaries r"
	
		write(11,*) '-------------------'
		write(11,*) px*py*pz
		write(11,*) px
		write(11,*) py
		write(11,*) pz
		write(11,*) '-------------------'
		
	endif
	
	do k=1,pz,1
		do j=1,py,1
			do i=1,px,1
				
				! indices (multiple and single) of the processor boxes
				mi=(k-1)*py*px+(j-1)*px+i	! long index of a processor box
		
				if (my_rank==0) then
				
					write(11,*) mi
					write(11,*) i,j,k
					
				endif
				
				! index boundaries of the processor boxes
				min_index(mi,1)=(i-1)*bnx+1-int(real(1/i))
				min_index(mi,2)=(j-1)*bny+1-int(real(1/j))
				min_index(mi,3)=(k-1)*bnz+1-int(real(1/k))
				
				if (i/=px) then
					max_index(mi,1)=i*bnx
				else
					max_index(mi,1)=nx
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
				
				if (my_rank==0) then
				
					write(11,*) min_index(mi,1), max_index(mi,1)
					write(11,*) min_index(mi,2), max_index(mi,2)
					write(11,*) min_index(mi,3), max_index(mi,3)
				
					! physical boundaries of the processor boxes
					write(11,*) x_min+min_index(mi,1)*dx, x_min+max_index(mi,1)*dx
					write(11,*) y_min+min_index(mi,2)*dy, y_min+max_index(mi,2)*dy
					write(11,*) z_min+min_index(mi,3)*dz, z_min+max_index(mi,3)*dz
					write(11,*) '-------------------'
					
				endif
			enddo
		enddo
	enddo
				
	if (my_rank==0) then
	
		close(UNIT=11)
	
	endif

!==============================================================================
! allocate memory
!==============================================================================


	ALLOCATE(rhoinv(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_rho)
	ALLOCATE(mu(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_mu)
	ALLOCATE(lambda(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_lambda)
	ALLOCATE(A(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_A)
	ALLOCATE(B(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_B)
	ALLOCATE(C(min_index(my_rank+1,1):max_index(my_rank+1,1), min_index(my_rank+1,2):max_index(my_rank+1,2), &
		(nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))),STAT=status_C)
	
	IF (status_A*status_B*status_C*status_rho*status_mu*status_lambda /= 0) THEN
		WRITE(*,*) 'Error: Processor ', my_rank, ': Could not allocate memory.'
	ELSE
		
		!=============================================================
		! homogeneous model
		!=============================================================
	
		IF (model_type==1) THEN
			
			rhoinv(:,:,:)=1.0/rho_homogen
			mu(:,:,:)=mu_homogen
			lambda(:,:,:)=lambda_homogen
			A(:,:,:)=A_homogen
			B(:,:,:)=B_homogen
			C(:,:,:)=C_homogen
			
		END IF
		
		!=============================================================
		! PREM isotropic (without ocean)
		!=============================================================
		
		IF (model_type==2) THEN
		
			ALLOCATE(z((nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))))
			ALLOCATE(x((nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))))
			ALLOCATE(rho_1D((nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))))
			ALLOCATE(vp_1D((nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))))
			ALLOCATE(vs_1D((nz-max_index(my_rank+1,3)):(nz-min_index(my_rank+1,3))))
		
			DO k=(nz-max_index(my_rank+1,3)),(nz-min_index(my_rank+1,3))
				z(k)=(z_max-k*dz)
				x(k)=z(k)/6371000
				
				! upper crust
				IF ((z(k)<=6371000) .AND. (z(k)>6356000)) THEN	
					rho_1D(k)=2600.0
					vp_1D(k)=5800.0
					vs_1D(k)=3200.0
				! lower crust
				ELSE IF ((z(k)<=6356000) .AND. (z(k)>6346600)) THEN
					rho_1D(k)=2900.0
					vp_1D(k)=6800.0
					vs_1D(k)=3900.0
				! seismic lithosphere
				ELSE IF ((z(k)<=6346600) .AND. (z(k)>6291000)) THEN
					rho_1D(k)=1000.0*(2.6910+0.6924*x(k))
					vp_1D(k)=1000.0*(4.1875+3.9382*x(k))
					vs_1D(k)=1000.0*(2.1519+2.3481*x(k))
				! low velocity zone
				ELSE IF ((z(k)<=6291000) .AND. (z(k)>6151000)) THEN
					rho_1D(k)=1000.0*(2.6910+0.6924*x(k))
					vp_1D(k)=1000.0*(4.1875+3.9382*x(k))
					vs_1D(k)=1000.0*(2.1519+2.3481*x(k))
				! transition zone 1 (depth from 220.0 km to 400.0 km)
				ELSE IF ((z(k)<=6151000) .AND. (z(k)>5971000)) THEN
					rho_1D(k)=1000*(7.1089-3.8045*x(k))
					vp_1D(k)=1000.0*(20.3926-12.2569*x(k))
					vs_1D(k)=1000.0*(8.9496-4.4597*x(k))
				! transition zone 2 (depth from 400.0 km to 600.0 km)
				ELSE IF ((z(k)<=5971000) .AND. (z(k)>5771000)) THEN
					rho_1D(k)=1000.0*(11.2494-8.0298*x(k))
					vp_1D(k)=1000.0*(39.7027-32.6166*x(k))
					vs_1D(k)=1000.0*(22.3512-18.5856*x(k))
				! transition zone 3 (depth from 600.0 km to 670.0 km
				ELSE IF ((z(k)<=5771000) .AND. (z(k)>5701000)) THEN
					rho_1D(k)=1000.0*(5.3197-1.4836*x(k))
					vp_1D(k)=1000.0*(19.0957-9.8672*x(k))
					vs_1D(k)=1000.0*(9.9839-4.9324*x(k))
				! lower mantle 1 (radius from 5600.0 km to 5701.0 km)
				ELSE IF ((z(k)<=5701000) .AND. (z(k)>5600000)) THEN
					rho_1D(k)=1000.0*(7.9565-6.4761*x(k)+5.5283*x(k)*x(k)-3.0807*x(k)*x(k)*x(k))
					vp_1D(k)=1000.0*(29.2766-23.6027*x(k)+5.5242*x(k)*x(k)-2.5514*x(k)*x(k)*x(k))
					vs_1D(k)=1000.0*(22.3459-17.2473*x(k)-2.0834*x(k)*x(k)+0.9783*x(k)*x(k)*x(k))
				! lower mantle 2 (radius from 5600.0 km to 5701.0 km)
				ELSE IF ((z(k)<=5600000) .AND. (z(k)>3630000)) THEN
					rho_1D(k)=1000.0*(7.9565-6.4761*x(k)+5.5283*x(k)*x(k)-3.0807*x(k)*x(k)*x(k))
					vp_1D(k)=1000.0*(24.9520-40.4673*x(k)+51.4832*x(k)*x(k)-26.6419*x(k)*x(k)*x(k))
					vs_1D(k)=1000.0*(11.1671-13.7818*x(k)+17.4575*x(k)*x(k)-9.2777*x(k)*x(k)*x(k))
				! lower mantle 3 (radius from 3480.0 km to 3630.0 km)
				ELSE IF ((z(k)<=3630000) .AND. (z(k)>3480000)) THEN
					rho_1D(k)=1000.0*(7.9565-6.4761*x(k)+5.5283*x(k)*x(k)-3.0807*x(k)*x(k)*x(k))
					vp_1D(k)=1000.0*(15.3891-5.3181*x(k)+5.5242*x(k)*x(k)-2.5514*x(k)*x(k)*x(k))
					vs_1D(k)=1000.0*(6.9254+1.4672*x(k)-2.0834*x(k)*x(k)+0.9783*x(k)*x(k)*x(k))
				! outer core (radius from 1221.5 km to 3480.0 km)
				ELSE IF ((z(k)<=3480000) .AND. (z(k)>1221500)) THEN
					rho_1D(k)=1000.0*(12.5815-1.2638*x(k)-3.6426*x(k)*x(k)-5.5281*x(k)*x(k)*x(k))
					vp_1D(k)=1000.0*(11.0487-4.0362*x(k)+4.8023*x(k)*x(k)-13.5732*x(k)*x(k)*x(k))
					vs_1D(k)=0
				! inner core
				ELSE IF (z(k)<=1221500) THEN
					rho_1D(k)=1000.0*(13.0885-8.8381*x(k)*x(k))
					vp_1D(k)=1000.0*(11.2622-6.3640*x(k)*x(k))
					vs_1D(k)=1000.0*(3.6678-4.4475*x(k)*x(k))
				END IF
				
				rhoinv(:,:,k)=1/rho_1D(k)
				mu(:,:,k)=rho_1D(k)*vs_1D(k)*vs_1D(k)
				lambda(:,:,k)=rho_1D(k)*vp_1D(k)*vp_1D(k)-2*mu(:,:,k)
				
				A(:,:,:)=0
				B(:,:,:)=0
				C(:,:,:)=0
			END DO
		
			DEALLOCATE(z)
			DEALLOCATE(x)
			DEALLOCATE(rho_1D)
			DEALLOCATE(vs_1D)
			DEALLOCATE(vp_1D)
			
		END IF
		
		
		
!==============================================================================
! save to files
!==============================================================================

		!==============================================================
		! open files, formatted or unformatted
		!==============================================================

		call int2str(my_rank,junk)

		if (is_formatted==1) then
		
			OPEN(UNIT=11,FILE='rhoinv'//junk(1:len_trim(junk)),ACTION='WRITE')
			OPEN(UNIT=12,FILE='mu'//junk(1:len_trim(junk)),ACTION='WRITE')
			OPEN(UNIT=13,FILE='lambda'//junk(1:len_trim(junk)),ACTION='WRITE')
			
		else
			
			OPEN(UNIT=11,FILE='rhoinv'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
			OPEN(UNIT=12,FILE='mu'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
			OPEN(UNIT=13,FILE='lambda'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
			
		endif
		
		if (is_aniso==1) then
			
			if (is_formatted==1) then
			
				OPEN(UNIT=14,FILE='A'//junk(1:len_trim(junk)),ACTION='WRITE')
				OPEN(UNIT=15,FILE='B'//junk(1:len_trim(junk)),ACTION='WRITE')
				OPEN(UNIT=16,FILE='C'//junk(1:len_trim(junk)),ACTION='WRITE')
				
			else
				
				OPEN(UNIT=14,FILE='A'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
				OPEN(UNIT=15,FILE='B'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
				OPEN(UNIT=16,FILE='C'//junk(1:len_trim(junk)),ACTION='WRITE',form='unformatted')
				
			endif
			
		endif
		
		!==============================================================
		! write parametres to files, formatted or unformatted
		!==============================================================
		
		if (is_formatted==1) then
		
			DO k=(nz-max_index(my_rank+1,3)),(nz-min_index(my_rank+1,3))
				DO j=min_index(my_rank+1,2), max_index(my_rank+1,2)
					DO i=min_index(my_rank+1,1), max_index(my_rank+1,1)
					
						WRITE(11,*) rhoinv(i,j,k)
						WRITE(12,*) mu(i,j,k)
						WRITE(13,*) lambda(i,j,k)
					
						if (is_aniso==1) then
					
							WRITE(14,*) A(i,j,k)
							WRITE(15,*) B(i,j,k)
							WRITE(16,*) C(i,j,k)
						
						endif
					
					END DO
				END DO
			END DO
			
		else
			
			write(11) rhoinv(:,:,:)
			write(12) mu(:,:,:)
			write(13) lambda(:,:,:)
			
			if (is_aniso==1) then
				
				write(14) A(:,:,:)
				write(15) B(:,:,:)
				write(16) C(:,:,:)
				
			endif
			
		endif
		
		CLOSE(UNIT=11)
		CLOSE(UNIT=12)
		CLOSE(UNIT=13)
		
		if (is_aniso==1) then
		
			CLOSE(UNIT=14)
			CLOSE(UNIT=15)
			CLOSE(UNIT=16)
			
		endif

!==============================================================================
! clean up
!==============================================================================
	
		DEALLOCATE(rhoinv)
		DEALLOCATE(mu)
		DEALLOCATE(lambda)
		
		DEALLOCATE(A)
		DEALLOCATE(B)
		DEALLOCATE(C)

	END IF
	
	deallocate(min_index)
	deallocate(max_index)
	
	call MPI_Finalize(ierr)

END PROGRAM generate_models



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
		
			!write(*,*) c
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
