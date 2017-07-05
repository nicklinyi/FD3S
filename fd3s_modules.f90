!*****************************************************************************
!***************** Declaration of all global variables ***********************
!*****************************************************************************

! Predefined parameters

MODULE parameters
IMPLICIT NONE

	INTEGER, PARAMETER :: nx_max=106
	INTEGER, PARAMETER :: ny_max=106
	INTEGER, PARAMETER :: nz_max=106
	
	INTEGER, PARAMETER :: maxnt=13000	! maximum number of time steps
	INTEGER, PARAMETER :: maxnr=500		! maximum number of receivers
	INTEGER, PARAMETER :: nabs_max=50	! maximum width of the damping region
	
	INTEGER, PARAMETER :: bndry=2		! width of the overlap region
	
	REAL, PARAMETER :: pi=3.1415926535898

END MODULE parameters
    
!==============================================================================

! Initial values

MODULE variables
USE parameters
IMPLICIT NONE

	integer :: nx, ny, nz		! local box sizes
	integer :: ix_multi, iy_multi, iz_multi

	INTEGER :: index_x, index_y, index_z
	INTEGER :: p, my_rank, ierr, px, py, pz

!==============================================================================

! File Names
 
    CHARACTER(LEN=50) :: recfile,seisfile,outfile

!==============================================================================

! Miscellaneous Variables

    INTEGER :: i,j,k,ix,iy,iz   ! index variables
    INTEGER :: icheck,isamp     ! output increments (see Par)

!==============================================================================
 
! Time parameters
 
    INTEGER :: it,nt            ! time index, number of time steps
    REAL :: dt                  ! time increment 

!==============================================================================

! Source Variables 

    REAL :: dip,rake,strike
    INTEGER :: isx,isy,isz,zsrc                   ! indices for source location
    INTEGER :: ij(3,3),srcproc,source_type,it0    ! node of source location
    real :: is_source, source_processor

    REAL, DIMENSION(1:maxnt+1) :: so ! source time function
    REAL :: moment(4:9)         ! moment tensor elements
    REAL :: xxs,yys,zzs         ! source coord. (see Par)
    REAL :: zzzs                ! source radius

!==============================================================================

! Geometrical variables

    REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: sininv1,sininv2,cottheta1,cottheta2

    REAL :: xs(-bndry:nx_max+bndry),ys(-bndry:ny_max+bndry),zs(-bndry:nz_max+bndry)     ! coordinates (rad,rad,m)
    REAL :: dx,dy,dz,dxinv,dyinv,dzinv     ! (inverse) grid spacing 
    REAL :: cc,ccinv                       ! conversion factors rad-deg 
    REAL :: xmin,xmax,ymin,ymax,zmin,zmax  ! (see Par)
    REAL :: xmin_global, xmax_global, ymin_global, ymax_global, zmin_global, zmax_global
    REAL :: xrange,yrange,zrange           ! model space range

!==============================================================================

! Free Surface Variables

    !INTEGER :: izfree,fs_model             ! (see Par)

! Graves surface implementation 
    
    real :: R1(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
    real :: R2(-bndry:nx_max+bndry,-bndry:ny_max+bndry)
    real :: R3(-bndry:nx_max+bndry,-bndry:ny_max+bndry)

!==============================================================================

! Absorbing boundaries

    INTEGER :: ibound,nabs                        ! (see Par)
    REAL :: window_z(-bndry:nx_max+bndry,-bndry:ny_max+bndry,0:nabs_max)
    REAL :: window_y_left(-bndry:nx_max+bndry,0:nabs_max,-bndry:nz_max+bndry)
    REAL :: window_y_right(-bndry:nx_max+bndry,0:nabs_max,-bndry:nz_max+bndry)
    REAL :: window_x_left(0:nabs_max,-bndry:ny_max+bndry,-bndry:nz_max+bndry)
    REAL :: window_x_right(0:nabs_max,-bndry:ny_max+bndry,-bndry:nz_max+bndry)

!==============================================================================

! Elastic Parameters & Model variables

    REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: rhoinv,mu,lam      ! 3D elastic parameters
    REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: mu_tau, lambda_tau ! mu_tau=mu*(1+tau_mu), lambda_tau=lam*(1+tau_lambda) 
    !INTEGER :: elast_inter
    
!==============================================================================

! Anelasticity and Anisotropy Parameters

	INTEGER :: is_aniso, is_diss, is_homogeneous, is_formatted
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: A, B, C				! anisotropy parameters
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: tau_lambda, tau_mu		! relaxation times
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: tau_s1_lambda, tau_s2_lambda	! relaxation times
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: tau_s1_mu, tau_s2_mu		! relaxation times
	
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: L1_11, L1_22, L1_33	! 1st memory variables, lambda
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: L2_11, L2_22, L2_33	! 2nd memory variables, lambda
	
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: M1_11, M1_22, M1_33, M1_12, M1_13, M1_23	! 1st memory variables, mu
	REAL, DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: M2_11, M2_22, M2_33, M2_12, M2_13, M2_23	! 1st memory variables, mu

!==============================================================================

! Velocity & Stress fields

    REAL,DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: w11,w12,w13 ! vel.comp. theta(t), phi(p), r
    REAL,DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: w14,w15,w16 ! norm.stress(t,t),(p,p),(r,r)
    REAL,DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: w17,w18,w19 ! shear stress(t,p),(t,r),(p,r)
    REAL,DIMENSION(-bndry:nx_max+bndry,-bndry:ny_max+bndry,-bndry:nz_max+bndry) :: w21,w22,w23 ! acceleration theta,phi,r

    REAL :: rotsurf(-bndry:nx_max+bndry,-bndry:ny_max+bndry),divsurf(-bndry:nx_max+bndry,-bndry:ny_max+bndry)  ! curl & divergence, surface 
    REAL :: rotx(-bndry:ny_max+bndry,-bndry:nz_max+bndry),divx(-bndry:ny_max+bndry,-bndry:nz_max+bndry)        ! curl & divergence, t=const 
    REAL :: roty(-bndry:nx_max+bndry,-bndry:nz_max+bndry),divy(-bndry:nx_max+bndry,-bndry:nz_max+bndry)        ! curl & divergence, p=const

!==============================================================================

! Receivers 
       
    INTEGER :: iflush,rotdiv        				! (see Par)
    INTEGER :: nr_global, nr                   			! number of receivers (first line of recfile)
    REAL :: recloc_global(1:2,1:maxnr), recloc(1:2,1:maxnr)     ! receiver location theta(deg),phi(deg) 
    
    REAL, DIMENSION(1:maxnr,1:maxnt) :: s_x, s_y, s_z		! seismograms

! Receiver interpolation variables, see fd3s_init.f90

    INTEGER :: ilowx(maxnr),iupx(maxnr),jlowx(maxnr),jupx(maxnr)
    INTEGER :: ilowy(maxnr),iupy(maxnr),jlowy(maxnr),jupy(maxnr)
    INTEGER :: ilowz(maxnr),iupz(maxnr),jlowz(maxnr),jupz(maxnr)   
    REAL :: xqotx(maxnr),xqotminx(maxnr),yqotx(maxnr),yqotminx(maxnr)
    REAL :: xqoty(maxnr),xqotminy(maxnr),yqoty(maxnr),yqotminy(maxnr)
    REAL :: xqotz(maxnr),xqotminz(maxnr),yqotz(maxnr),yqotminz(maxnr) 
 
!==============================================================================

! Snapshots

    INTEGER :: ssamp,iplanex,iplaney,iplanez ! (see Par)
    

!==============================================================================
! adjoint calculations
!==============================================================================

	character(len=50) :: bndryfile
	integer, dimension(1:maxnt) :: is_bndry_x_left, is_bndry_x_right
	integer, dimension(1:maxnt) :: is_bndry_y_left, is_bndry_y_right
	integer, dimension(1:maxnt) :: is_bndry_z_lower
	integer :: zero_ratio
	real :: m

END MODULE variables
!==============================================================================
