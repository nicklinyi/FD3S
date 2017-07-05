!*****************************************************************************
!***************** Initialization of FD operator, model space,****************
!******************* boundaries, receiver location, source *******************
!*****************************************************************************

SUBROUTINE fd3s_init
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'

	! variables for receiver interpolation
    	REAL :: xxloc, yyloc, D
	REAL, DIMENSION(1:maxnr) :: hyx, hxx, hyy, hxy
	REAL, DIMENSION(1:maxnr) :: hy, hx
	
	! dummy variables
	real :: sxd, syd, szd
	character(len=20) :: dummy
	integer :: status(MPI_STATUS_SIZE)
    
	D=0.3
	
    	if (my_rank==0) then
    
	    write(*,*)'begin init '
    	    
	endif
	
	WRITE(99,*)'begin init '
	
	!======================================================================
	! Initialization of memory variables and modified elastic parameters
	!======================================================================
	
	IF (is_diss==1) THEN
	
		L1_11(:,:,:)=0
		L1_22(:,:,:)=0
		L1_33(:,:,:)=0
	
		L2_11(:,:,:)=0
		L2_22(:,:,:)=0
		L2_33(:,:,:)=0
	
		M1_11(:,:,:)=0
		M1_22(:,:,:)=0
		M1_33(:,:,:)=0
	
		M2_11(:,:,:)=0
		M2_22(:,:,:)=0
		M2_33(:,:,:)=0
	
		M1_12(:,:,:)=0
		M1_13(:,:,:)=0
		M1_23(:,:,:)=0
	
		M2_12(:,:,:)=0
		M2_13(:,:,:)=0
		M2_23(:,:,:)=0
		
	END IF
	
	lambda_tau=lam+lam*tau_lambda
	mu_tau=mu+mu*tau_mu
    
	!======================================================================
	! initialize source moment tensor
	!======================================================================

	!write(*,*) "test1"
	
	! DAS MUSS ICH SOWIESO NOCHMAL CHECKEN
	
	moment(4)= -SIN(dip)*COS(rake)*SIN(2*strike) &     ! M_theta,theta
     	        -SIN(2*dip)*SIN(rake)*SIN(strike)**2

	moment(5)=  SIN(dip)*COS(rake)*SIN(2*strike) &     ! M_phi,phi	
     		-SIN(2*dip)*SIN(rake)*COS(strike)**2

	moment(6)=  SIN(2*dip)*SIN(rake)                   ! M_r,r

	moment(7)=  SIN(dip)*COS(rake)*COS(2*strike) &     ! M_theta,phi
     		+.5*SIN(2*dip)*SIN(rake)*SIN(2*strike)

	moment(8)= -COS(dip)*COS(rake)*COS(strike)  &      ! M_theta,r
     		-COS(2*dip)*SIN(rake)*SIN(strike) 

	moment(9)= -COS(dip)*COS(rake)*SIN(strike)  &      ! M_phi,r
     		+COS(2*dip)*SIN(rake)*COS(strike)

	!======================================================================
	! Initialization of geometrical variables    
	!======================================================================

	!write(*,*) "test2"
	
	! convert to radians
	cc=pi/180.         
	ccinv=1./cc
	xmin=xmin*cc  ! minimal theta value, i.e. northermost
	xmax=xmax*cc
	ymin=ymin*cc  ! minimal phi value, i.e. westernmost
	ymax=ymax*cc

	! define model space range
	xrange=xmax-xmin   ! (rad)
	yrange=ymax-ymin   ! (rad)
	zrange=zmax-zmin   ! (m)

	! define grid spacing
	dx=xrange/nx     ! (rad)   
	dy=yrange/ny     ! (rad)   
	dz=zrange/nz     ! (m)  

	! inverse grid spacing for faster FD computation
	dxinv=1./dx  
	dyinv=1./dy  
	dzinv=1./dz  

	write(99,*) 'dx=',dx ,'dy=',dy ,'dz=',dz
	
	!======================================================================
	! coordinates
	!======================================================================

	!write(*,*) "test3"
	
	! in theta direction
	
	call int2str(my_rank, dummy)

	OPEN(unit=30,file='DATA/COORDINATES/xco_'//dummy(1:len_trim(dummy)),status='unknown')
   
	DO ix=-bndry,nx+bndry 
		xs(ix)=xmin+ix*dx;              ! theta coord. axis of model space (rad)
	ENDDO
	DO ix=0,nx
		WRITE(30,*) xs(ix)
	ENDDO

	CLOSE(30)
	
	! in phi direction

    	OPEN(unit=30,file='DATA/COORDINATES/yco_'//dummy(1:len_trim(dummy)),status='unknown')
     
	DO iy=-bndry,ny+bndry
		ys(iy)=ymin+iy*dy              ! phi coord. axis of model space (rad)
	ENDDO
	DO iy=0,ny
		WRITE(30,*) ys(iy)
	ENDDO
	
	CLOSE(30)

	! in vertical direction
	
	OPEN(unit=30,file='DATA/COORDINATES/zco_'//dummy(1:len_trim(dummy)),status='unknown')


    	DO iz=-bndry,nz+bndry 
		zs(iz)=zmax-(iz-0)*dz;	! vert. coord. axis of one node (m), conversion from radius to depth
	ENDDO
	DO iz=0,nz
        	WRITE(30,*)zs(iz)
	ENDDO
    
	CLOSE(30)

	!======================================================================
	! initialization of 3D spherical scaling factors for faster FD computations
	!======================================================================
 
	!write(*,*) "test4"
	
    	DO i=-bndry,nx+bndry
		sininv1(i,:,:)=1./SIN(xs(i))                    ! for sigma_(i,i),  
		cottheta1(i,:,:)=COS(xs(i))*sininv1(i,:,:)      !sigma_(phi,r),v_phi,v_r
		sininv2(i,:,:)=1./SIN(xs(i)+dx/2)               ! for sigma_(theta,r)  
		cottheta2(i,:,:)=COS(xs(i)+dx/2)*sininv2(i,:,:) !sigma_(theta,phi),v_theta
	ENDDO

	!======================================================================
	! Free surface, initialization of R1, R2, R3
	!======================================================================

	!write(*,*) "test5"
	
	R1(:,:)=1.0
	R2(:,:)=1.0
	R3(:,:)=1.0
	
	IF (is_aniso==1) THEN
		IF (is_diss==1) THEN
			R1=(lambda_tau(:,:,0)+2*mu_tau(:,:,0))*dzinv+(lambda_tau(:,:,0)+C(:,:,0))/(zmax_global+dz)
			R2=(lambda_tau(:,:,0)+2*mu_tau(:,:,0))*dzinv
			R3=-(lambda_tau(:,:,0)+C(:,:,0))/(zmax)
		ELSE
			R1=(lam(:,:,0)+2*mu(:,:,0))*dzinv+(lam(:,:,0)+C(:,:,0))/(zmax_global+dz)
			R2=(lam(:,:,0)+2*mu(:,:,0))*dzinv
			R3=-(lam(:,:,0)+C(:,:,0))/(zmax_global)
		END IF
	ELSE
		IF (is_diss==1) THEN
			R1=(lambda_tau(:,:,0)+2*mu_tau(:,:,0))*dzinv+(lambda_tau(:,:,0))/(zmax_global+dz)
			R2=(lambda_tau(:,:,0)+2*mu_tau(:,:,0))*dzinv
			R3=-(lambda_tau(:,:,0))/(zmax_global)
		ELSE
			R1=(lam(:,:,0)+2*mu(:,:,0))*dzinv+(lam(:,:,0))/(zmax_global+dz)
			R2=(lam(:,:,0)+2*mu(:,:,0))*dzinv
			R3=-(lam(:,:,0))/(zmax_global)
		END IF
	END IF


	!======================================================================
	! receiver interpolation
	!======================================================================
	
	recloc_global=recloc_global*cc
	
	k=1
	do i=1,nr_global
		
		if ((xmin<=recloc_global(1,i)) .and. (xmax>recloc_global(1,i)) .and. &
		    (ymin<=recloc_global(2,i)) .and. (ymax>recloc_global(2,i))) then
		    
		    recloc(1,k)=recloc_global(1,i)
		    recloc(2,k)=recloc_global(2,i)
		    k=k+1
		    
		endif
		
	enddo
	
	nr=k-1
	
	write(99,*) 'number of receivers', nr
	!write(*,*) 'processor ', my_rank, ', number of receivers', nr

	! calculate receiver location interpolation values  

       	do k=1,nr                  						! for each receiver

        	xxloc=recloc(1,k)					! theta coord. of receivers
		yyloc=recloc(2,k) 					! phi coord. of receivers
		!write(*,*) xxloc, xs(1)
		do i=0,nx-1               					! interpolation in theta direction
    
			!~~~~~~~~ theta velocity component ~~~~~~~~~~~~~~~~~~~~~
     
			if((xxloc>=(xs(i)+dx/2)) .and. (xxloc<(xs(i+1)+dx/2))) then

				ilowx(k)=i              			! neighboring grid point index to the north
				iupx(k)=i+1             			! neighboring grid point index to the south
				hxx(k)=xxloc-xs(i)-dx/2 			! distance receiver<->northern neighbor (rad), distance to the grid point of v_theta
			
				xqotx(k)=hxx(k)/dx    				! v_th: relat. distance to northern neighbor 
				xqotminx(k)=1-xqotx(k)  			! v_th: relat. distance to southern neighbor
				!write(*,*) i
		
			endif

			!~~~~~~~~ phi velocity component ~~~~~~~~~~~~~~~~~~~~~~~
			!~~~~~~~~ and vertical component ~~~~~~~~~~~~~~~~~~~~~~~
                                  
			if((xxloc .GE. xs(i)) .AND. (xxloc .LT. xs(i+1))) then 

				ilowy(k)=i              			! neighboring grid point index to the north
				iupy(k) =i+1            			! neighboring grid point index to the south
				hxy(k)=xxloc-xs(i)      			! distance receiver<->northern neighbor (rad)
			
				xqoty(k)=hxy(k)/dx    				! v_ph: relat. distance to northern neighbor
				xqotminy(k)=1-xqoty(k)       			! v_ph: relat. distance to southern neighbor
		
				
				ilowz(k)=i              			! neighboring grid point index to the north
				iupz(k) =i+1            			! neighboring grid point index to the south
				hx(k)=xxloc-xs(i)       			! distance receiver<->northern neighbor (rad)
			
				xqotz(k)=hx(k)/dx      				! v_r: relat. distance to northern neighbor
				xqotminz(k)=1-xqotz(k)  			! v_r: relat. distance to southern neighbor
			
				
			endif
       
		enddo  								! end search along theta direction

		
		do j=0,ny-1                  					! interpolation in phi direction

			!~~~~~~~~ theta velocity component ~~~~~~~~~~~~~~~~~~~~~
			!~~~~~~~~ and vertical component ~~~~~~~~~~~~~~~~~~~~~~~

         		if((yyloc .GE. ys(j)) .AND. (yyloc .LT. ys(j+1))) then

				jlowx(k)=j				        ! neighboring grid point index to the west
				jupx(k) =j+1            			! neighboring grid point index to the east

				hyx(k)=yyloc-ys(j)      			! distance receiver<->western neighbor (rad)
			
				yqotx(k)=hyx(k)/dy    				! v_th: relat. distance to western neighbor 
				yqotminx(k)=1-yqotx(k)  			! v_th: relat. distance to eastern neighbor 
				
				jlowz(k)=j              			! neighboring grid point index to the west
				jupz(k) =j+1            			! neighboring grid point index to the east
				hy(k)=yyloc-ys(j)       			! distance receiver<->western neighbor (rad)
			
				yqotz(k)=hy(k)/dy      				! v_r: relat. distance to western neighbor 
				yqotminz(k)=1-yqotz(k)       			! v_r: relat. distance to eastern neighbor 
		
			endif     

			!~~~~~~~~ phi velocity component ~~~~~~~~~~~~~~~~~~~~~~~ 
                                  
			if((yyloc .GE. ys(j)+dy/2) .AND. (yyloc .LT. ys(j+1)+dy/2)) then

				jlowy(k)=j              			! neighboring grid point index to the west
				jupy(k) =j+1            			! neighboring grid point index to the east
				hyy(k)=yyloc-ys(j)-dy/2 			! distance receiver<->western neighbor (rad)
			
				yqoty(k)=hyy(k)/dy    				! v_ph: relat. distance to western neighbor 
				yqotminy(k)=1-yqoty(k)       			! v_ph: relat. distance to eastern neighbor 

			endif
    
		enddo   							! end search along phi direction
                                 
	enddo  									! end loop over different receiver loactions
     

     	!=========================================================================
	! point source location
	!=========================================================================

	if (is_source==1) then
	
		isx=MAXVAL(MINLOC(ABS(xs-xxs/180.*pi)))-1-bndry        ! source theta index
		isy=MAXVAL(MINLOC(ABS(ys-yys/180.*pi)))-1-bndry        ! source phi index
		isz=MAXVAL(MINLOC(ABS(zmax_global-zzs-zs)))-1-bndry    ! source in z direction     
		
		write(99,*)'Source centered in grid point',isx,isy,isz
		write(99,*)'Source coordinates :',xs(isx)*180/pi,ys(isy)*180/pi,zs(isz)
		
		call mpi_send(isx,1,mpi_integer,0,1,mpi_comm_world,ierr)
		call mpi_send(isy,1,mpi_integer,0,2,mpi_comm_world,ierr)
		call mpi_send(isz,1,mpi_integer,0,3,mpi_comm_world,ierr)
		
		call mpi_send(xs(isx)*180/pi,1,mpi_real,0,4,mpi_comm_world,ierr)
		call mpi_send(ys(isy)*180/pi,1,mpi_real,0,5,mpi_comm_world,ierr)
		call mpi_send(zs(isz),1,mpi_real,0,6,mpi_comm_world,ierr)
		
	endif
	
	if (my_rank==0) then
		
		call mpi_recv(isx,1,mpi_integer,int(source_processor),1,mpi_comm_world,status,ierr)
		call mpi_recv(isy,1,mpi_integer,int(source_processor),2,mpi_comm_world,status,ierr)
		call mpi_recv(isz,1,mpi_integer,int(source_processor),3,mpi_comm_world,status,ierr)
		
		call mpi_recv(sxd,1,mpi_real,int(source_processor),4,mpi_comm_world,status,ierr)
		call mpi_recv(syd,1,mpi_real,int(source_processor),5,mpi_comm_world,status,ierr)
		call mpi_recv(szd,1,mpi_real,int(source_processor),6,mpi_comm_world,status,ierr)
	
		write(*,*)'Source centered in grid point',isx,isy,isz, 'of processor', int(source_processor)
		write(*,*)'Source coordinates :',sxd,syd,szd
		
	endif

	!======================================================================
	! damping windows
	!======================================================================
	
	! lower z boundary
	if (iz_multi==1) then
		window_z=1.0
		
		do i=0,nabs
			window_z(:,:,i)=window_z(:,:,i)*exp(-(D*i/nabs)*(D*i/nabs))
		enddo
	endif
	
	! left y boundary
	if (iy_multi==1) then
		window_y_left=1.0
	
		do i=0,nabs
			window_y_left(:,i,:)=window_y_left(:,i,:)*exp(-(D*(nabs-i)/nabs)*(D*(nabs-i)/nabs))
		enddo
	endif
	
	! right y boundary
	if (iy_multi==py) then
		window_y_right=1.0
	
		do i=0,nabs
			window_y_right(:,i,:)=window_y_right(:,i,:)*exp(-(D*i/nabs)*(D*i/nabs))
		enddo
	endif
	
	! left x boundary
	if (ix_multi==1) then
		window_x_left=1.0
	
		do i=0,nabs
			window_x_left(i,:,:)=window_x_left(i,:,:)*exp(-(D*(nabs-i)/nabs)*(D*(nabs-i)/nabs))
		enddo
	endif
	
	! right x boundary
	if (ix_multi==px) then
		window_x_right=1.0
	
		do i=0,nabs
			window_x_right(i,:,:)=window_x_right(i,:,:)*exp(-(D*i/nabs)*(D*i/nabs))
		enddo
	endif
	
	
	!======================================================================
	! clean up
	!======================================================================

	if (my_rank==0) then

		write(*,*) 'end init'
	
	endif
	
	write(99,*) 'end init'
	write(99,*) '----------------------------------------'

return
end subroutine fd3s_init
