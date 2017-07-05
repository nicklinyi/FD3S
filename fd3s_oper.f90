!*****************************************************************************
!***************** Interpolation & FD partial derivatives in 3D **************
!******************* for different operator lengths ****** *******************
!*****************************************************************************

! Interpolation 

SUBROUTINE inter(o,f,comp,dir)
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'


	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry), INTENT(in)  :: f
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry), INTENT(out) :: o

	CHARACTER(len=1) :: dir,comp

	o=0.

	IF ( comp == 'x' ) THEN

	    IF ( dir == '-' )THEN

		DO i=0,nx
		  o(i,:,:)=(f(i,:,:)+f(i-1,:,:))*.5
		ENDDO
	        !o(0,:,:)=0.

	    ELSEIF( dir == '+' )THEN

		DO i=0,nx
		  o(i,:,:)=(f(i+1,:,:)+f(i,:,:))*.5
		ENDDO
	        !o(nx,:,:)=0.

	    ENDIF

	ELSEIF( comp == 'y' ) THEN

	    IF ( dir == '-' ) THEN

		DO j=0,ny
		  o(:,j,:)=(f(:,j,:)+f(:,j-1,:))*.5
		ENDDO
		!o(:,0,:)=0.

	    ELSEIF( dir == '+' )THEN

		DO j=0,ny
		  o(:,j,:)=(f(:,j+1,:)+f(:,j,:))*.5
		ENDDO
		!o(:,ny,:)=0.

	    ENDIF

	ELSEIF( comp == 'z') THEN

	    IF ( dir == '-' ) THEN

		DO k=0,nz
      		  o(:,:,k)=(f(:,:,k)+f(:,:,k-1))*.5
		ENDDO
		!o(:,:,0)=0.

	    ELSEIF( dir == '+')THEN

		DO k=0,nz
		  o(:,:,k)=(f(:,:,k+1)+f(:,:,k))*.5
		ENDDO
		!o(:,:,nz)=0.

	    ENDIF

	ENDIF
	

END SUBROUTINE inter


!=============================================================================


! FD method: partial derivatives for different operator lengths (2<->8)

SUBROUTINE pder(h,f,comp,dir)
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'

	

	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry), INTENT(in) :: f
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry), INTENT(out) :: h
        
        CHARACTER(len=1) :: dir,comp
	
	!write(99,*) "pder debug 01"
	
	IF ( comp == 'x' ) THEN

		IF ( dir == '-' )THEN
			 
			!write(*,*) 'x -'
			
			DO i=0,nx
	       	  		h(i,:,:)=dxinv*(9/8)*(f(i,:,:)-f(i-1,:,:))-dxinv*(1/24)*(f(i+1,:,:)-f(i-2,:,:))
				!h(i,:,:)=17.5
			ENDDO
			!h(0:1,:,:)=0.
			!h(nx,:,:) =0.

		ELSEIF( dir == '+' )THEN

			DO i=0,nx
		  		h(i,:,:)=dxinv*(9/8)*(f(i+1,:,:)-f(i,:,:))-dxinv*(1/24)*(f(i+2,:,:)-f(i-1,:,:))
			ENDDO
			!h(0,:,:)=0.
			!h(nx-1:nx,:,:)=0.

		ENDIF

	!write(99,*) "pder debug 02"
		
	ELSEIF( comp == 'y' ) THEN

		IF ( dir == '-' ) THEN

			DO j=0,ny
				h(:,j,:)=dyinv*(9/8)*(f(:,j,:)-f(:,j-1,:))-dyinv*(1/24)*(f(:,j+1,:)-f(:,j-2,:))
			ENDDO
			!h(:,0:1,:)=0.
			!h(:,ny,:)=0.

		ELSEIF( dir == '+' )THEN

			DO j=0,ny
		  		h(:,j,:)=dyinv*(9/8)*(f(:,j+1,:)-f(:,j,:))-dyinv*(1/24)*(f(:,j+2,:)-f(:,j-1,:))
			ENDDO
			!h(:,0,:)=0.
			!h(:,ny-1:ny,:)=0.

		ENDIF

	!write(99,*) "pder debug 03"
		
	ELSEIF( comp == 'z') THEN

		IF ( dir == '-' ) THEN

			DO k=0,nz 
				!h(:,:,k)=dzinv*(9/8)*(f(:,:,k)-f(:,:,k-1))-dzinv*(1/24)*(f(:,:,k+1)-f(:,:,k-2))
				h(:,:,k)=dzinv*(9/8)*(f(:,:,k-1)-f(:,:,k))-dzinv*(1/24)*(f(:,:,k-2)-f(:,:,k+1))
			ENDDO
			!h(:,:,0:1)=0.
			!h(:,:,nz)=0.		

		ELSEIF( dir == '+')THEN

			DO k=0,nz
				!h(:,:,k)=dzinv*(9/8)*(f(:,:,k+1)-f(:,:,k))-dzinv*(1/24)*(f(:,:,k+2)-f(:,:,k-1))
				h(:,:,k)=dzinv*(9/8)*(f(:,:,k)-f(:,:,k+1))-dzinv*(1/24)*(f(:,:,k-1)-f(:,:,k+2))
			ENDDO
			!h(:,:,0)=0.
			!h(:,:,nz-1:nz)=0.

		ENDIF

	ENDIF
	
	!h(:,:,:)=17.5
	
END SUBROUTINE pder

!=============================================================================
