!*****************************************************************************
!***************** Output of snapshot files and seismograms ******************
!*****************************************************************************

subroutine fd3s_output(n)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! local variables
	!======================================================================
	
	integer, intent(in) :: n
	character(len=10) :: dummy, ns
	character(len=50) :: fn
	
	!======================================================================
	! write velocity slices
	!======================================================================
	
	call int2str(my_rank,dummy)
	call int2str(n,ns)
	
	! plane theta=const

	fn='DATA/OUTPUT/v_theta_thetaplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=10,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_phi_thetaplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=11,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_r_thetaplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=12,file=fn,action='WRITE')
	
	do k=0,nz
		do j=0,ny

			write(10,*) w11(iplanex,j,k)
			write(11,*) w12(iplanex,j,k)
			write(12,*) w13(iplanex,j,k)
			
		enddo
	enddo
	
	close(unit=10)
	close(unit=11)
	close(unit=12)
	
	! plane phi=const
	
	fn='DATA/OUTPUT/v_theta_phiplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=10,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_phi_phiplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=11,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_r_phiplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=12,file=fn,action='WRITE')
	
	do k=0,nz
		do i=0,nx

			write(10,*) w11(i,iplaney,k)
			write(11,*) w12(i,iplaney,k)
			write(12,*) w13(i,iplaney,k)
			
		enddo
	enddo
	
	close(unit=10)
	close(unit=11)
	close(unit=12)
	
	! plane z=const
	
	fn='DATA/OUTPUT/v_theta_rplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=10,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_phi_rplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=11,file=fn,action='WRITE')
	fn='DATA/OUTPUT/v_r_rplane'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
	open(unit=12,file=fn,action='WRITE')
	
	do j=0,ny
		do i=0,nx

			write(10,*) w11(i,j,iplanez)
			write(11,*) w12(i,j,iplanez)
			write(12,*) w13(i,j,iplanez)
			
		enddo
	enddo
	
	close(unit=10)
	close(unit=11)
	close(unit=12)

	!======================================================================
	! write seismograms
	!======================================================================
	
	if ((iz_multi==pz) .and. (nr>0)) then
	
		open(unit=51,file='DATA/OUTPUT/seismogram_theta'//dummy(1:len_trim(dummy)),action='WRITE')
		open(unit=52,file='DATA/OUTPUT/seismogram_phi'//dummy(1:len_trim(dummy)),action='WRITE')
		open(unit=53,file='DATA/OUTPUT/seismogram_r'//dummy(1:len_trim(dummy)),action='WRITE')
	
		write(51,*) 'theta component seismograms'
		write(52,*) 'phi component seismograms'
		write(53,*) 'r component seismograms'
	
		write(51,*) 'nr=', nr
		write(51,*) 'nt=', nt
		write(51,*) 'dt=', dt
	
		write(52,*) 'nr=', nr
		write(52,*) 'nt=', nt
		write(52,*) 'dt=', dt
	
		write(53,*) 'nr=', nr
		write(53,*) 'nt=', nt
		write(53,*) 'dt=', dt
	
		do i=1,nr
	
			write(51,*) 'theta=', recloc(1,i), 'phi=', recloc(2,i)
			write(52,*) 'theta=', recloc(1,i), 'phi=', recloc(2,i)
			write(53,*) 'theta=', recloc(1,i), 'phi=', recloc(2,i)
		
			do j=1,nt
		
				write(51,*) s_x(i,j)
				write(52,*) s_y(i,j)
				write(53,*) s_z(i,j)
			
			enddo
		
		enddo
	
		close(unit=51)
		close(unit=52)
		close(unit=53)
		
	endif

end subroutine fd3s_output


!==============================================================================
! record seismograms
!==============================================================================

subroutine record_seismograms
use parameters
use variables
implicit none
include 'mpif.h'

	integer :: depth
	
	depth=0;

	do i=1,nr
		
		! theta velocity component interpolation
	
		!s_x(i,it)=w11(ilowx(i),jlowx(i),0)
		s_x(i,it)=(w11(ilowx(i),jlowx(i),depth)*xqotminx(i) + &
			   w11(iupx(i),jlowx(i),depth)*xqotx(i))*yqotminx(i) + &
			  (w11(ilowx(i),jupx(i),depth)*xqotminx(i) + &
			   w11(iupx(i),jupx(i),depth)*xqotx(i))*yqotx(i)
			    
		! phi velocity component interpolation
			  
		s_y(i,it)=(w12(ilowy(i),jlowy(i),depth)*xqotminy(i) + &
			   w12(iupy(i),jlowy(i),depth)*xqoty(i))*yqotminy(i) + &
			  (w12(ilowy(i),jupy(i),depth)*xqotminy(i) + &
			   w12(iupy(i),jupy(i),depth)*xqoty(i))*yqoty(i)
			   
		! z velocity component interpolation
		
		s_z(i,it)=(0.5*(w13(ilowz(i),jlowz(i),depth)+w13(ilowz(i),jlowz(i),depth-1))*xqotminz(i) + &
			   0.5*(w13(iupz(i),jlowz(i),depth)+w13(iupz(i),jlowz(i),depth-1))*xqotz(i))*yqotminz(i) + &
			  (0.5*(w13(ilowz(i),jupz(i),depth)+w13(ilowz(i),jupz(i),depth-1))*xqotminz(i) + &
			   0.5*(w13(iupz(i),jupz(i),depth)+w13(iupz(i),jupz(i),depth-1))*xqotz(i))*yqotz(i)
			   
	enddo
		
end subroutine record_seismograms


!==============================================================================
! save velocity blocks, unformatted
!==============================================================================

subroutine save_velocity(n)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! some variables
	!======================================================================

	integer, intent(in) :: n
	character(len=10) :: dummy, ns
	
	!======================================================================
	! write velocity boxes
	!======================================================================
	
	call int2str(my_rank,dummy)
	call int2str(n,ns)

	open(unit=10,file='DATA/OUTPUT/v_theta_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=11,file='DATA/OUTPUT/v_phi_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=12,file='DATA/OUTPUT/v_r_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	
	write(10) w11(0:nx,0:ny,0:nz)
	write(11) w12(0:nx,0:ny,0:nz)
	write(12) w13(0:nx,0:ny,0:nz)
	
	close(unit=10)
	close(unit=11)
	close(unit=12)

end subroutine save_velocity


!==============================================================================
! save stress blocks, unformatted
!==============================================================================

subroutine save_stress(n)
use parameters
use variables
implicit none
include 'mpif.h'

	!======================================================================
	! some variables
	!======================================================================

	integer, intent(in) :: n
	character(len=10) :: dummy, ns
	
	!======================================================================
	! write velocity boxes
	!======================================================================
	
	call int2str(my_rank,dummy)
	call int2str(n,ns)

	open(unit=10,file='DATA/OUTPUT/s_theta_theta_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=11,file='DATA/OUTPUT/s_theta_phi_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=12,file='DATA/OUTPUT/s_theta_r_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=13,file='DATA/OUTPUT/s_phi_phi_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=14,file='DATA/OUTPUT/s_phi_r_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	open(unit=15,file='DATA/OUTPUT/s_r_r_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns)),action='WRITE',form='unformatted')
	
	write(10) w14(0:nx,0:ny,0:nz)
	write(11) w17(0:nx,0:ny,0:nz)
	write(12) w18(0:nx,0:ny,0:nz)
	write(13) w15(0:nx,0:ny,0:nz)
	write(14) w19(0:nx,0:ny,0:nz)
	write(15) w16(0:nx,0:ny,0:nz)
	
	close(unit=10)
	close(unit=11)
	close(unit=12)
	close(unit=13)
	close(unit=14)
	close(unit=15)

end subroutine save_stress
