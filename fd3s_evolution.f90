! Noch zu checken:
!
! - Addition der Momententensorkomponenten
! - Die Interpolation von mu für die Nebendiagonalelemente von sigma war mir so schleierhaft, dass ich sie eliminiert habe! Wäre mal noch zu überprüfen.

!*****************************************************************************
!***************** FD time evolution and spatial derivatives******************
!*****************************************************************************

SUBROUTINE fd3s_evolution
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'

	!======================================================================
	! temporary variables, CAUTION: variables take different meanings during the routine
	!======================================================================
	
        REAL :: exx(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: eyy(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: ezz(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: exy(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: tmp1(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: tmp2(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: tmp3(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
        REAL :: tmp4(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
	REAL :: dummy(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
	
	character(len=20) :: message
	
	exx(:,:,:)=0.0
	eyy(:,:,:)=0.0
	ezz(:,:,:)=0.0
	exy(:,:,:)=0.0
	tmp1(:,:,:)=0.0
	tmp2(:,:,:)=0.0
	tmp3(:,:,:)=0.0
	tmp4(:,:,:)=0.0
	dummy(:,:,:)=0.0

	!======================================================================
	! communication and tapering
	!======================================================================
	
	call comm_tap_s

	!======================================================================
	! add source(s)
	!======================================================================
	
	if (is_source==1) then
	
		IF(source_type==4)THEN   ! Explosion
         		w14(isx,isy,isz)=w14(isx,isy,isz)+so(it)
			w15(isx,isy,isz)=w15(isx,isy,isz)+so(it)
			w16(isx,isy,isz)=w16(isx,isy,isz)+so(it)
		ELSEIF(source_type==10)THEN   ! Moment tensor source 
			w14(isx,isy,isz)=w14(isx,isy,isz)+moment(4)*so(it)        
			w15(isx,isy,isz)=w15(isx,isy,isz)+moment(5)*so(it)
			w16(isx,isy,isz)=w16(isx,isy,isz)+moment(6)*so(it)
			w17(isx,isy,isz)=w17(isx,isy,isz)+moment(7)*so(it)
			w18(isx,isy,isz)=w18(isx,isy,isz)+moment(8)*so(it)
			w19(isx,isy,isz)=w19(isx,isy,isz)+moment(9)*so(it)
		ENDIF
	endif
	
	!======================================================================
	! free surface boundary conditions for stress tensor elements:
	! antisymmetry condition for sigma_rr, sigma_(r,phi), sigma_(r,theta)
	! only if the processor box is located at the surface
	!======================================================================
	
	if (iz_multi==pz) then
	
		w16(:,:,0)=0.
		do i=1,2
        		w16(:,:,0-i)=-w16(:,:,0+i)		! sigma_(r,r)
			w18(:,:,0-i)=-w18(:,:,0-1+i)		! sigma_(theta,r)
			w19(:,:,0-i)=-w19(:,:,0-1+i)		! sigma_(phi,r)
		enddo
		
	endif
	
	!======================================================================
	! calculate acceleration from stresses (divergence of the stress tensor in spherical coordinates)
	!message="   acceleration"
	!call write_cpu_time_log(message)
	!======================================================================

	! d/dt(v_theta)*rho  (equal to w21 here)

    	CALL pder(exx,w14,'x','+')     ! derivatives d/dtheta of s_theta,theta
	CALL pder(eyy,w17,'y','-')     ! derivatives d/dphi of s_theta,phi
	CALL pder(ezz,w18,'z','-')     ! derivatives d/dr of s_theta,r
	CALL inter(tmp1,w14,'x','+')   ! interpolate s_(theta,theta) for v_theta
	CALL inter(tmp2,w15,'x','+')   ! interpolate s_(phi,phi) for v_theta
	CALL inter(tmp3,w18,'z','-')   ! interpolate s_(r,theta) for v_theta

	DO iz=0,nz	! w21=div[sigma]_(theta)
     		w21(0:nx,0:ny,iz) =	ezz(0:nx,0:ny,iz) + 1/zs(iz)*exx(0:nx,0:ny,iz) + &
					1/zs(iz)*sininv2(0:nx,0:ny,iz)*eyy(0:nx,0:ny,iz) + &
					1/zs(iz)*((tmp1(0:nx,0:ny,iz)-tmp2(0:nx,0:ny,iz))*cottheta2(0:nx,9:ny,iz)+ &
					3.*tmp3(0:nx,0:ny,iz)) 
	ENDDO

	! d/dt(v_phi)*rho (equal to w22 here)
        
    	CALL pder(exx,w17,'x','-')     ! derivatives d/dtheta of s_phi,theta
	CALL pder(eyy,w15,'y','+')     ! derivatives d/dphi of s_phi,phi
	CALL pder(ezz,w19,'z','-')     ! derivatives d/dr of s_phi,r
	CALL inter(tmp1,w19,'z','-')   ! interpolate s_r,phi for v_phi
	CALL inter(tmp2,w17,'x','-')   ! interpolate s_theta,phi for v_phi
        
	DO iz=0,nz     ! w22=div[sigma]_(phi)
     		w22(0:nx,0:ny,iz) = 	ezz(0:nx,0:ny,iz) + 1/zs(iz)*exx(0:nx,0:ny,iz) + &
					1/zs(iz)*sininv1(0:nx,0:ny,iz)*eyy(0:nx,0:ny,iz) + &
					1/zs(iz)*(3.*tmp1(0:nx,0:ny,iz) + &
					2.*cottheta1(0:nx,0:ny,iz)*tmp2(0:nx,0:ny,iz))
	ENDDO
	
	! d/dt(v_r)*rho (equal to w23 here)
        
	CALL pder(exx,w18,'x','-')     ! derivatives d/dtheta of s_r,theta
	CALL pder(eyy,w19,'y','-')     ! derivatives d/dphi of s_r,phi
	CALL pder(ezz,w16,'z','+')     ! derivatives d/dr of s_r,r
	CALL inter(tmp1,w14,'z','+')   ! interpolate s_theta,theta for v_r
	CALL inter(tmp2,w15,'z','+')   ! interpolate s_phi,phi for v_r
	CALL inter(tmp3,w16,'z','+')   ! interpolate s_r,r for v_r
	CALL inter(tmp4,w18,'x','-')   ! interpolate s_r,theta for v_r

	DO iz=0,nz 	! w23=div[sigma]_(r)
		w23(0:nx,0:ny,iz)  = 	ezz(0:nx,0:ny,iz)  + 1/(zs(iz)-dz/2)*exx(0:nx,0:ny,iz) + &
					1/(zs(iz)-dz/2)*sininv1(0:nx,0:ny,iz) *eyy(0:nx,0:ny,iz) + & 
					1/(zs(iz)-dz/2)*(2.*tmp3(0:nx,0:ny,iz) -tmp1(0:nx,0:ny,iz) - &
					tmp2(0:nx,0:ny,iz) + cottheta1(0:nx,0:ny,iz) *tmp4(0:nx,0:ny,iz) )
	ENDDO

	!======================================================================
	! add single forces
	!======================================================================
	
	if (is_source==1) then
	
		!write(*,*) 'source:', isx, isy, isz, so(it)
		
		IF(source_type==1)THEN                       ! Theta single force
          		w21(isx,isy,isz)=w21(isx,isy,isz)+so(it)
		ELSEIF(source_type==2)THEN                   ! Phi single force
          		w22(isx,isy,isz)=w22(isx,isy,isz)+so(it)
		ELSEIF(source_type==3)THEN                   ! Vertical single force
          		w23(isx,isy,isz)=w23(isx,isy,isz)+so(it)
		ENDIF
		
	endif
	
	!======================================================================
	! 2nd-order time extrapolation of velocities
	!======================================================================
	
	w11=w11+rhoinv*w21*dt
	w12=w12+rhoinv*w22*dt
	w13=w13+rhoinv*w23*dt

!------------------------------------------------------------------------------
	
	!======================================================================
	! communication and tapering
	!message="   comm_tap_v"
	!call write_cpu_time_log(message)
	!======================================================================
	
	call comm_tap_v
	
	!======================================================================
	! free surface velocity components
	!message="   graves"
	!call write_cpu_time_log(message)
	!======================================================================
	
	if (iz_multi==pz) then
	
		call fs_graves(w11(:,:,0-1),w12(:,:,0-1),w13(:,:,0-1))
	
	endif
	
	!======================================================================
	! strain rate from spatial derivatives of new velocities
	!message="   strain rate"
	!call write_cpu_time_log(message)
	!======================================================================
	
	! d/dt(e_theta,theta) (equal to exx here)

    	call pder(w21,w11,'x','-')
	call inter(tmp1,w13,'z','-')
 
	do iz=0,nz
		exx(0:nx,0:ny,iz)=1/zs(iz)*(w21(0:nx,0:ny,iz)+tmp1(0:nx,0:ny,iz))
	enddo

	! d/dt(e_phi,phi) (equal to eyy here)

    	call pder(w22,w12,'y','-')
	call inter(tmp2,w11,'x','-')

	do iz=0,nz 
		eyy(0:nx,0:ny,iz)=	1/zs(iz)*(sininv1(0:nx,0:ny,iz)*w22(0:nx,0:ny,iz)+ &
					tmp1(0:nx,0:ny,iz)+cottheta1(0:nx,0:ny,iz)*tmp2(0:nx,0:ny,iz))
	enddo

	! d/dt(e_r,r) (equal to ezz here) 

    	call pder(ezz,w13,'z','-')						

	! ezz in free surface calculated via 2nd order derivative
	
	if (iz_multi==pz) then
		ezz(0:nx,0:ny,0)=(w13(0:nx,0:ny,0)-w13(0:nx,0:ny,0-1))*dzinv
	endif

	! d/dt(e_theta,phi) (equal to exy here)

    	call pder(w21,w11,'y','+')
	call pder(w22,w12,'x','+')
	call inter(tmp1,w12,'x','+')

	do iz=0,nz 
		exy(0:nx,0:ny,iz)=	0.5*1/zs(iz)*(sininv2(0:nx,0:ny,iz)*w21(0:nx,0:ny,iz)+ &
					w22(0:nx,0:ny,iz)- cottheta2(0:nx,0:ny,iz)*tmp1(0:nx,0:ny,iz))  
	enddo

	! d/dt(e_theta,r) first part

    	call pder(w23,w11,'z','+')
	call pder(tmp3,w13,'x','+')   
	call inter(tmp1,w11,'z','+')

	! Continuation d/dt(e_theta,r) (equal to tmp3 here)
  
    	do iz=0,nz  
		tmp3(0:nx,0:ny,iz)=	0.5*(1/(zs(iz)-dz/2)*tmp3(0:nx,0:ny,iz)+w23(0:nx,0:ny,iz)- &
                   			1/(zs(iz)-dz/2)*tmp1(0:nx,0:ny,iz))      
	enddo

	! d/dt(e_phi,r) (equal to tmp4 here)

	call pder(tmp2,w12,'z','+')   ! For e_phi,r already (also for curl)
	call pder(tmp4,w13,'y','+')   ! For e_phi,r already (also for curl)
    	call inter(tmp1,w12,'z','+')

	do iz=0,nz
		tmp4(0:nx,0:ny,iz)=	0.5*(1/(zs(iz)-dz/2)*sininv1(0:nx,0:ny,iz)*tmp4(0:nx,0:ny,iz)+ &
                   			tmp2(0:nx,0:ny,iz)-1/(zs(iz)-dz/2)*tmp1(0:nx,0:ny,iz))
	enddo
	
	!======================================================================
	! Hookes Relation & final second order time extrapolation (stress)
	! Normal stresses (Located in same grid points as Lame's Parameters)
	!message="   Hooke"
	!call write_cpu_time_log(message)
	!======================================================================
	
	if (is_diss==1) then
		dummy=L1_11+L2_11+L1_22+L2_22+L1_33+L2_33
		
		if (is_aniso==1) then

			w14=w14+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*exx+C*ezz+A*(eyy+exx)+2*(M1_11+M2_11)+dummy)*dt
			w15=w15+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*eyy+C*ezz+A*(eyy+exx)+2*(M1_22+M2_22)+dummy)*dt
			w16=w16+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*ezz+C*(eyy+exx)+2*(M1_22+M2_22)+dummy)*dt
		
		else
			w14=w14+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*exx+2*(M1_11+M2_11)+dummy)*dt
			w15=w15+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*eyy+2*(M1_22+M2_22)+dummy)*dt
			w16=w16+(lambda_tau*(exx+eyy+ezz)+2*mu_tau*ezz+2*(M1_22+M2_22)+dummy)*dt

		endif
		
	else
		if (is_aniso==1) then
				
			w14=w14+(lam*(exx+eyy+ezz)+2*mu*exx+C*ezz+A*(eyy+exx))*dt
			w15=w15+(lam*(exx+eyy+ezz)+2*mu*eyy+C*ezz+A*(eyy+exx))*dt
			w16=w16+(lam*(exx+eyy+ezz)+2*mu*ezz+C*(eyy+exx))*dt
			
		else
	
			w14=w14+(lam*(exx+eyy+ezz)+2*mu*exx)*dt
			w15=w15+(lam*(exx+eyy+ezz)+2*mu*eyy)*dt
			w16=w16+(lam*(exx+eyy+ezz)+2*mu*ezz)*dt
			
		endif
	endif
		
	!======================================================================
	! Forward time step of diagonal memory veriables
	!======================================================================

	if (is_diss==1) then
		L1_11=L1_11-dt*(0.5*lam*tau_lambda*exx+L1_11)/tau_s1_lambda
		L1_22=L1_22-dt*(0.5*lam*tau_lambda*eyy+L1_22)/tau_s1_lambda
		L1_33=L1_33-dt*(0.5*lam*tau_lambda*ezz+L1_33)/tau_s1_lambda
		L2_11=L2_11-dt*(0.5*lam*tau_lambda*exx+L2_11)/tau_s2_lambda
		L2_22=L2_22-dt*(0.5*lam*tau_lambda*eyy+L2_22)/tau_s2_lambda
		L2_33=L2_33-dt*(0.5*lam*tau_lambda*ezz+L2_33)/tau_s2_lambda
		
		M1_11=M1_11-dt*(0.5*mu*tau_mu*exx+M1_11)/tau_s1_mu
		M1_22=M1_22-dt*(0.5*mu*tau_mu*eyy+M1_22)/tau_s1_mu
		M1_33=M1_33-dt*(0.5*mu*tau_mu*ezz+M1_33)/tau_s1_mu
		M2_11=M2_11-dt*(0.5*mu*tau_mu*exx+M2_11)/tau_s2_mu
		M2_22=M2_22-dt*(0.5*mu*tau_mu*eyy+M2_22)/tau_s2_mu
		M2_33=M2_33-dt*(0.5*mu*tau_mu*ezz+M2_33)/tau_s2_mu
		
		write (*,*) 'memory forward 1'
	endif

	!======================================================================
	! Hookes Relation & final second order time extrapolation (stress)
	! Shear stresses (Located in same grid points as Lame's Parameters)
	!======================================================================

	if (is_diss==1) then
		
		if (is_aniso==1) then
			
			w19=w19+(2*mu_tau*tmp4+2*B*tmp4+2*(M1_23+M2_23))*dt		! sigma(r,phi)
			w18=w18+(2*mu_tau*tmp3+2*B*tmp3+2*(M1_13+M2_13))*dt		! sigma(r,theta)
			w17=w17+(2*mu_tau*exy+2*(M1_12+M2_12))*dt			! sigma(phi,theta)
			
		else
			
			w19=w19+(2*mu_tau*tmp4+2*(M1_23+M2_23))*dt		! sigma(r,phi)
			w18=w18+(2*mu_tau*tmp3+2*(M1_13+M2_13))*dt		! sigma(r,theta)
			w17=w17+(2*mu_tau*exy+2*(M1_12+M2_12))*dt		! sigma(phi,theta)
			
		endif
		
	else
		
		if (is_aniso==1) then
	
			w19=w19+(2*mu*tmp4+2*B*tmp4)*dt		! sigma(r,phi)
			w18=w18+(2*mu*tmp3+2*B*tmp3)*dt		! sigma(r,theta)
			w17=w17+(2*mu*exy)*dt			! sigma(phi,theta)
			
		else
			
			w19=w19+(2*mu*tmp4)*dt		! sigma(r,phi)
			w18=w18+(2*mu*tmp3)*dt		! sigma(r,theta)
			w17=w17+(2*mu*exy)*dt		! sigma(phi,theta)
			
		endif
		
	endif
	
	!======================================================================
	! Forward time step for off-diagonal memory variables
	!======================================================================
	
	if (is_diss==1) then
		
		M1_12=M1_12-dt*(0.5*mu*tau_mu*exy+M1_12)/tau_s1_mu
		M1_13=M1_13-dt*(0.5*mu*tau_mu*tmp3+M1_13)/tau_s1_mu
		M1_23=M1_23-dt*(0.5*mu*tau_mu*tmp4+M1_23)/tau_s1_mu
		M2_12=M2_12-dt*(0.5*mu*tau_mu*exy+M2_12)/tau_s2_mu
		M2_13=M2_13-dt*(0.5*mu*tau_mu*tmp3+M2_13)/tau_s2_mu
		M2_23=M2_23-dt*(0.5*mu*tau_mu*tmp4+M2_23)/tau_s2_mu
		
	endif

	!call write_cpu_time(proc)
	
	!======================================================================
	! Free surface of normal stresses
	!======================================================================
       
      	!w14(:,:,0:izfree-1)=0.
	!w15(:,:,0:izfree-1)=0.
	
	!message="   finish"
	!call write_cpu_time_log(message)
	
	if (my_rank==0) then
		write(*,*) 'end evolution'
	endif

END SUBROUTINE fd3s_evolution



!########### Additional subroutines ##########################################
! CALL fs_graves(w11(:,:,izfree-1),w12(:,:,izfree-1),w13(:,:,izfree-1)) 

subroutine fs_graves(fs_x,fs_y,fs_z)
use parameters
use variables
implicit none
include 'mpif.h'
 
    	real :: fs_x(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
	real :: fs_y(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
	real :: fs_z(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
	real :: out1(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
	real :: out2(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
	real :: A1, A2
    
	fs_x(:,:)=0.0
	fs_y(:,:)=0.0
	fs_z(:,:)=0.0
	out1(:,:)=0.0
	out2(:,:)=0.0
	
	A1=dzinv-1/(2*(zmax_global+dz))
	A2=dzinv-1/(2*(zmax_global-dz))

	!======================================================================
	! Vertical Component dz/2 above free surface
	!======================================================================
	
	call pder2(out1,w11(:,:,0),'-','x')	! d_theta v_theta
	call pder2(out2,w12(:,:,0),'-','y')	! d_phi v_phi

	if (is_diss==1) then
		fs_z=R2(0:nx,0:ny)*w13(0:nx,0:ny,0)/R1 + &
    			R3(0:nx,0:ny)*(sininv1(0:nx,0:ny,0)*out2(0:nx,0:ny)+w13(0:nx,0:ny,0)+out1(0:nx,0:ny)+ &
			cottheta1(0:nx,0:ny,0)*w11(0:nx,0:ny,0))/R1(0:nx,0:ny) - &
			(L1_11(0:nx,0:ny,0)+L1_22(0:nx,0:ny,0)+L1_33(0:nx,0:ny,0)+ & 
			L2_11(0:nx,0:ny,0)+L2_22(0:nx,0:ny,0))/R1(0:nx,0:ny) + &
			(L2_33(0:nx,0:ny,0)+2*M1_33(0:nx,0:ny,0)+2*M2_33(0:nx,0:ny,0))/R1(0:nx,0:ny)
	else
		fs_z(0:nx,0:ny)=R2(0:nx,0:ny)*w13(0:nx,0:ny,0)/R1(0:nx,0:ny) + &
    			R3(0:nx,0:ny)*(sininv1(0:nx,0:ny,0)*out2(0:nx,0:ny)+w13(0:nx,0:ny,0)+ &
			out1(0:nx,0:ny)+cottheta1(0:nx,0:ny,0)*w11(0:nx,0:ny,0))/R1(0:nx,0:ny)
	endif
    
	!======================================================================
	! Theta-Component dz above free surface
	!======================================================================

	call pder2(out1,w13(:,:,0),'+','x')	! d_theta v_r at the surface
	call pder2(out2,fs_z,'+','x')		! d_theta v_r above the surfac
    
	if (is_diss==1) then
	    fs_x=(w11(:,:,0)/zmax-out1/zmax-out2/zmax)/A1+A2*w11(:,:,0+1)/A1 - &
    		 4*(M1_13(:,:,0)+M2_13(:,:,0))/(A1*(mu_tau(:,:,0)+B(:,:,0)))
	else
	    fs_x=(w11(:,:,0)/zmax-out1/zmax-out2/zmax)/A1+A2*w11(:,:,0+1)/A1
	endif
	
	!======================================================================
	! Phi-Component dz above free surface
	!======================================================================

	call pder2(out1,w13(:,:,0),'+','y')	! d_phi v_r at the surface
	call pder2(out2,fs_z,'+','y')		! d_phi v_r above the surface
    
	if (is_diss==1) then
	    fs_y=(w12(:,:,0)/zmax-out1*sininv1(:,:,0)/zmax-out2*sininv1(:,:,0)/zmax)/A1 + &
    		 A2*w12(:,:,0+1)/A1 - &
		 4*(M1_12(:,:,0)+M2_12(:,:,0))/(A1*(mu_tau(:,:,0)+B(:,:,0)))
	else
	    fs_y=(w12(:,:,0)/zmax-out1*sininv1(:,:,0)/zmax-out2*sininv1(:,:,0)/zmax)/A1 + &
    		 A2*w12(:,:,0+1)/A1
	endif
	
	return
end subroutine fs_graves


!#############################################################################
! subroutine for partial FD derivatives in 2D , particularly for the surface


SUBROUTINE pder2(h2,f2,dir,comp)
USE parameters
USE variables
IMPLICIT NONE
include 'mpif.h'

	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry), INTENT(in) :: f2
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry), INTENT(out) :: h2        
	CHARACTER*1 dir,comp

	IF ( comp == 'x' ) THEN         ! spatial dimension of derivative 

		IF ( dir == '-' )THEN        ! direction of difference operator 

			DO i=0,nx
				h2(i,:)=((9/8)*(f2(i,:)-f2(i-1,:))-(1/24)*(f2(i+1,:)-f2(i-2,:)))*dxinv;
			ENDDO
			!h2(0:1,:)=0.
			!h2(nx,:) =0.
		
		ELSEIF( dir == '+' )THEN

			DO i=0,nx
				h2(i,:)=((9/8)*(f2(i+1,:)-f2(i,:))-(1/24)*(f2(i+2,:)-f2(i-1,:)))*dxinv;
			ENDDO
			!h2(0,:)=0.
			!h2(nx-1:nx,:)=0.
			
		ENDIF		

	ELSEIF( comp == 'y' ) THEN

		IF ( dir == '-' ) THEN

			DO j=0,ny
				h2(:,j)=((9/8)*(f2(:,j)-f2(:,j-1))-(1/24)*(f2(:,j+1)-f2(:,j-2)))*dyinv
			ENDDO
			!h2(:,0:1)=0.
			!h2(:,ny)=0.
			
		ELSEIF( dir == '+' )THEN
	
			DO j=0,ny
				h2(:,j)=((9/8)*(f2(:,j+1)-f2(:,j))-(1/24)*(f2(:,j+2)-f2(:,j-1)))*dyinv
			ENDDO
			!h2(:,0)=0.
			!h2(:,ny-1:ny)=0.
			
		ENDIF		
			
	ENDIF

END SUBROUTINE pder2

!#############################################################################
