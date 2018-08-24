MODULE LAGRANGE_MOD
USE PARAMS_MOD
USE GRID_2D
IMPLICIT NONE
INTEGER, PARAMETER :: d_per_s = 86400   ! How many seconds of one day
INTEGER, PARAMETER :: h_per_s = 3600    ! How many seconds of one hour
INTEGER            :: ct   = 0, ctd=0   ! Current time and time-of-the-day

! Summer noon above surface irradiance (W m-2) around Oshima
REAL,    PARAMETER :: PAR0 = 300.    
REAL,    PARAMETER :: kext = 0.06   ! Attenuation coefficient (m-1) of water
REAL,    PARAMETER :: dayL = 50400. ! Day length (14 hours)

! The shallowest vertical location for zoo.
REAL,    PARAMETER :: Z_MAX     =  -10.
REAL,    PARAMETER :: pi        = 3.1415926535897932384633

! Time to save
INTEGER, PARAMETER :: nsave     =    h_per_s/int(dt)

INTEGER, PARAMETER :: Nstep     =    d_per_s*NDays/int(dt)

! Current surface PAR in terms of diel cycle
REAL               :: SPAR      =    PAR0  

TYPE Particle

    ! the mth group
    integer :: m

    ! the nth particle in the mth group
    integer :: n

    ! Grid indices for particles
    integer :: ix, iz
    
    ! coordinates for particles
    real    :: rx, rz
    
    ! associated u,w,and kap (vertical diffusivity)
    real    :: u, Uw, kap

    ! swimming speed (m/s)
    real    :: w = 0.

    ! Maximal swimming speed (m/s)
    real    :: wm= 3.5D-2  ! m/s (Batchelder et al. PIO 2002)

    ! Prefered optimal light (isolume)
    real    :: iso= 1D2 ! W m-2
END TYPE Particle

! Positions of particles
TYPE (Particle)  :: Zp(TNPAR)

CONTAINS

SUBROUTINE Init_Particles
IMPLICIT NONE
INTEGER :: i,j,k,n
REAL    :: rnd
REAL    :: bom ! A scratch variable for bottom depth
REAL    :: z2  ! A scratch z variable


!Minimal and maximal X grid
INTEGER, parameter :: x0= 2050, x1= 3000

!Minimal and maximal depths
REAL,    parameter :: z0=-200., z1=-10.

k = 1

! Assign particles (randomly initialize Zp and Ip)
DO i = 1, NG
   DO j = 1, N_PAR
      Zp(k)%m = i
      Zp(k)%n = j

      !Set optimal light (isolume)
      SELECT CASE(i)
      CASE(1)  !No DVM (control)
        Zp(k)%wm  = 0d0
      CASE(2)  !Prefer high light
        Zp(k)%iso = 15.
      CASE(3)  !Prefer around 200 m
        Zp(k)%iso = 0.04
      CASE(4)  !Prefer around 500 m
        Zp(k)%iso = 2D-11
      CASE DEFAULT
        stop "There is no such kind of particle group!"
      END SELECT

      IF (i .EQ. 1) THEN

         CALL RANDOM_NUMBER(rnd)

         !Randomly distributed in one grid
         !Different groups start from the same position

         Zp(k)%ix= INT(rnd * DBLE(x1-x0) + DBLE(x0))

         ! Compute rx
         Zp(k)%rx= X_r(Zp(k)%ix)

         CALL RANDOM_NUMBER(rnd)

         ! Vertical position must be above bottom depth
         ! Find corresponding bottom depth
         bom     = Bot_dep(Zp(k)%ix, 2)
         z2      = max(z0, bom+5.)
         Zp(k)%rz= rnd * (z1-z2) + z2

         ! Find iz
         DO n=1, NZ
            IF (Z_w(n-1) < Zp(k)%rz .AND. Z_w(n) >= Zp(k)%rz) THEN
               Zp(k)%iz = n
               EXIT
            ENDIF
         END DO

      ELSE
         ! Track the same particle in the 1st group
         Zp(k)%iz= Zp(j)%iz
         Zp(k)%rz= Zp(j)%rz
         Zp(k)%ix= Zp(j)%ix
         Zp(k)%rx= Zp(j)%rx
      ENDIF

      ! Calculate zooplankton vertical swimming speed
      Zp(k)%w = swim_speed(SPAR,Zp(k)%rz, Zp(k)%iso, bom)&
              * Zp(k)%wm

      k = k + 1
   ENDDO
ENDDO
RETURN
END SUBROUTINE Init_Particles

! !ROUTINE: Lagrangian particle random walk at 2D space (Horizontal and vertical)
SUBROUTINE lagrange2D(nz,nx,zlev,xlev,mask,kvz,kvx,w,u,npar,zoo)
!
! !DESCRIPTION:
!
! Here a Lagrangian particle random walk for spatially
! inhomogeneous turbulence according to \cite{Visser1997} is implemented.
! With the random walk, the particle $i$ is moved from the vertical
! position $z_i^n$ to $z_i^{n+1}$ according to the following algorithm:
! \begin{equation}
! \begin{array}{rcl}
! z_i^{n+1} &=&
! z^n_i + \partial_z \nu_t (z^n_i)\Delta t \\ \\
! &+&
! R \left\{2 r^{-1} \nu_t (z^n_i + \frac12  \partial_z \nu_t (z^n_i)\Delta t)
! \Delta t\right\}^{1/2},
! \end{array}
! \end{equation}
! where $R$ is a random process with $\langle R \rangle =0$ (zero mean) and
! and the variance $\langle R^2 \rangle=r$.
! Set {\tt visc\_corr=.true.} for
! evaluating eddy viscosity in a semi-implicit way. A background viscosity
! ({\tt visc\_back}) may be set. The variance $r$ of the random walk scheme
! ({\tt rnd\_var}) has to be set manually as well here.
!
! !USES:
!
IMPLICIT NONE
!
! !INPUT PARAMETERS:
! Number of vertical points
INTEGER, intent(in)                :: nz

! Number of horizontal points
INTEGER, intent(in)                :: nx

! Vertical z coordinates from bottom to surface
REAL,    intent(in)                :: zlev(0:nz)

! X coordinates from left to right 
REAL,    intent(in)                :: xlev(0:nx)

! mask of each grid
INTEGER, intent(in)                :: mask(nx,nz)

! Vertical Eddy diffusivities 
REAL,    intent(in)                :: kvz(nx,0:nz)

! Horizontal Eddy diffusivities 
REAL,    intent(in)                :: kvx(0:nx,nz)

! Vertical velocity of water
REAL,    intent(in)                :: w(1:nx,0:nz)

! Horizontal velocity of water
REAL,    intent(in)                :: u(0:nx,1:nz)

! Number of particles
INTEGER, intent(in)                :: npar

! INPUT/OUTPUT PARAMETERS:
! Index for particle positions
TYPE (Particle), intent(inout)     :: ZOO(npar)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
INTEGER              :: i,n,k
REAL, ALLOCATABLE    :: vel(:), nu(:) 
INTEGER, ALLOCATABLE :: mask_(:)
!EOP
!-----------------------------------------------------------------------
!BOC

! Save u, w, and kappat
DO n = 1, npar
   zoo(n)%u   =   u(zoo(n)%ix, zoo(n)%iz)
   zoo(n)%Uw  =   w(zoo(n)%ix, zoo(n)%iz)
   zoo(n)%kap = kvz(zoo(n)%ix, zoo(n)%iz)

   !Check whether zoo is in a fake grid (mask ==0)
   if (mask(zoo(n)%ix, zoo(n)%iz) .eq. 0) then
       write(6,*) 'ix = ', zoo(n)%ix
       write(6,*) 'iz = ', zoo(n)%iz
       write(6,*) 'Bottom depth = ', bot_dep(zoo(n)%ix,2)
       stop "The particle is in a fake grid!"
   endif
ENDDO

! Vertical movement

ALLOCATE(vel(0:nz))
vel(:) = 0d0

ALLOCATE(mask_(nz))
mask_(:) = 0

ALLOCATE(nu(0:nz))
nu(:) = 0.

DO i = 1, nx
   DO n = 1,npar
      IF (zoo(n)%ix == i) THEN

        ! Vertical velocity including swimming/sinking
        do k = 1, (nz-1)
           vel(k) = zoo(n)%w + w(i,k)
            nu(k) = kvz(i,k)
        enddo
        do k = 1, nz
           mask_(k) = mask(i, k)
        enddo

        CALL LAGRANGE(nz,zlev,mask_,nu,vel,zoo(n)%iz, zoo(n)%rz)

      ENDIF
   ENDDO
ENDDO
DEALLOCATE(vel)
DEALLOCATE(nu)
DEALLOCATE(mask_)

! Horizontal movement
ALLOCATE(vel(0:nx))
vel(:) = 0d0
ALLOCATE(mask_(nx))
mask_(:)= 0
ALLOCATE(nu(0:nx))
nu(:)= 0.


DO i = 1, nz
   DO n = 1,npar
      IF (zoo(n)%iz == i) THEN

        ! Horizontal velocity 
        do k = 1, (nx-1)
           vel(k) = u(k,i)
            nu(k) = kvx(k,i)
        enddo

        do k = 1, nx
           mask_(k) = mask(k,i)
        enddo

        CALL LAGRANGE(nx,xlev,mask_,nu,vel, zoo(n)%ix, zoo(n)%rx)

      ENDIF
   ENDDO
ENDDO
DEALLOCATE(vel)
DEALLOCATE(nu)
DEALLOCATE(mask_)

RETURN
END SUBROUTINE LAGRANGE2D

! Lagrange particle tracking on one particle
! !INTERFACE:
SUBROUTINE LAGRANGE(nlev,zlev,mask,nuh,w,zi,zp)
!
IMPLICIT NONE
!
! !INPUT PARAMETERS:
INTEGER, intent(in)                :: nlev
REAL,    intent(in)                :: zlev(0:nlev)
INTEGER, intent(in)                :: mask(nlev)
REAL,    intent(in)                :: nuh(0:nlev)
REAL,    intent(in)                :: w(0:nlev)

! INPUT/OUTPUT PARAMETERS:
! the ith grid where the particle currently resides
INTEGER, intent(inout)             :: zi

! Z position of particles
REAL,    intent(inout)             :: zp
!
! !LOCAL VARIABLES:
INTEGER            :: i
REAL               :: rnd, rnd_var_inv
REAL,parameter     :: visc_back=0.e-6,rnd_var=0.333333333
REAL               :: dz(nlev),dzn(nlev),step,zp_old
REAL               :: visc,rat,dt_inv,zloc
logical,parameter  :: visc_corr=.TRUE.
!EOP
!-----------------------------------------------------------------------
!BOC

!Judge whether the input is correct
IF (zp < zlev(0) .or. zp > zlev(nlev)) THEN
   stop "The particle out of the domain!"
ENDIF

dt_inv      = 1./dt
rnd_var_inv = 1./rnd_var

CALL RANDOM_NUMBER(rnd)

rnd = (2.*rnd-1.)

do i = 1,nlev
   
   ! Depth of each grid
   dz(i) = zlev(i)-zlev(i-1)

   ! gradient of Kv
   dzn(i)= (nuh(i)-nuh(i-1))/dz(i)
end do


!  local viscosity calculation
IF (visc_corr) THEN ! correction suggested by Visser [1997]

    ! The new position after dt
   zloc = zp + 0.5*(dzn(zi)+w(zi))*dt 

   step = zloc - zp   ! The distance traveled

   ! Judge whether the particle overshoots through a grid. 
   ! If yes, the time step needs to be shortened
   if (zi .lt. nlev .and. zi .gt. 1) then
     if (zloc > zlev(zi+1) .or. zloc < zlev(zi-2)) then
       stop "The particle overshoots! Reduce time step!"
     endif
   endif

   ! First judge whether the particle jumps out of the domain
   if (zloc > zlev(nlev)) then
      zloc = 2.*zlev(nlev) - zloc 
      zi   = nlev
    ! Judge whether after bouncing the particle overshoots through a grid. 
      if (zloc < zlev(nlev-1)) then 
          stop "The particle overshoots! Reduce time step!"
      endif
   elseif (zloc < zlev(0)) then     
      zloc = 2d0*zlev(0) - zloc      
      zi   = 1
      if (zloc > zlev(1)) then 
          stop "The particle overshoots! Reduce time step!"
      endif
   endif      

   ! To judge whether the particle enters a fake grid (mask = 0)
   ! If yes, then bounce back
   if (step > 0.) then
      if (zloc > zlev(zi)) then
         if (mask(zi + 1) == 0) then
             zloc = 2d0*zlev(zi) - zloc
        ! Judge whether after bouncing the particle overshoots through a grid. 
             if (zloc < zlev(zi-1)) then 
                 stop "The particle overshoots! Reduce time step!"
             endif
         endif
      endif
   elseif (step < 0.) then
      if (zloc < zlev(zi-1)) then
         if (mask(zi - 1) == 0) then
             zloc = 2d0*zlev(zi-1) - zloc

        ! Judge whether after bouncing the particle overshoots through a grid. 
             if (zloc > zlev(zi)) then 
                 stop "The particle overshoots! Reduce time step!"
             endif
         endif
      endif
   endif

   step = zloc-zp     ! The distance traveled

! The following is to find the index for the grid where the particle currently resides
   if (step > 0.) then ! search new index above old index
      do i=zi, nlev
         if (zlev(i) .gt. zloc) then
             zi = i
             EXIT
         endif
      end do
   else                ! search new index below old index
      do i=zi,1,-1
         if (zlev(i-1) .lt. zloc) then
             zi=i
             EXIT
         endif
      end do
   end if
ELSE
   i   =zi
   zloc=zp
END IF

rat  = (zloc-zlev(i-1))/dz(i)
visc = rat*nuh(i)+(1.-rat)*nuh(i-1)  ! interpolate Kv

IF (visc.lt.visc_back) visc=visc_back   !Kv background

zp_old = zp

step = dt*( sqrt(2d0*rnd_var_inv*dt_inv*visc) * rnd + w(i) + dzn(i) )

zp = zp + step

! Judge whether the particle overshoots through a grid. 
! If yes, the time step needs to be shortened
if (i .lt. nlev .and. i .gt. 1) then
  if (zp > zlev(i+1) .or. zp < zlev(i-2)) then
    stop "The particle overshoots! Reduce time step!"
  endif
endif

! First judge whether the particle jumps out of the domain
if (zp > zlev(nlev)) then
   zp = 2.*zlev(nlev) - zp 
 ! Judge whether after bouncing the particle overshoots through a grid. 
   if (zp < zlev(nlev-1)) then 
       stop "The particle overshoots! Reduce time step!"
   endif
elseif (zp < zlev(0)) then     
   zp = 2d0*zlev(0) - zp      
   if (zp > zlev(1)) then 
       stop "The particle overshoots! Reduce time step!"
   endif
endif      

! To judge whether the particle enters a fake grid (mask = 0)
! If yes, then bounce back
if (step > 0.) then
   if (zp > zlev(i)) then
      if (mask(i + 1) == 0) then
          zp = 2.*zlev(i) - zp
     ! Judge whether after bouncing the particle overshoots through a grid. 
          if (zp < zlev(i-1)) then 
              stop "The particle overshoots! Reduce time step!"
          endif
      endif
   endif
elseif (step < 0.) then
   if (zp < zlev(i-1)) then
      if (mask(i - 1) == 0) then
          zp = 2.*zlev(i-1) - zp

     ! Judge whether after bouncing the particle overshoots through a grid. 
          if (zp > zlev(i)) then 
              stop "The particle overshoots! Reduce time step!"
          endif
      endif
   endif
endif

step = zp - zp_old

if (step.gt.0) then ! search new index above old index
   do i=zi,nlev
      if (zlev(i) .gt. zp) EXIT
   end do
else                ! search new index below old index
   do i=zi,1,-1
      if (zlev(i-1) .lt. zp) EXIT
   end do
end if
zi = i

return
END SUBROUTINE LAGRANGE

REAL FUNCTION surface_par(timeofday)
IMPLICIT NONE
INTEGER, INTENT(IN) :: timeofday
real,    parameter  :: Trise= 21600. ! Timing of sun rise
real,    parameter  :: Tset = 72000. ! Timing of sun set


! Initialize surface PAR based on diel cycle
IF (timeofday < Trise .or. timeofday > Tset) then
    surface_par = 0.
else
    surface_par = sin(pi*dble(timeofday-Trise)/dayL)*PAR0
endif

END FUNCTION surface_par

!Relative swimming speed (-1, 1)
REAL FUNCTION swim_speed(PAR_0,z, isolume, bot_dep)
IMPLICIT NONE
REAL, INTENT(IN) :: PAR_0, isolume, z, bot_dep
REAL             :: PAR_, swimdirection, z_des, diff, dis

! Half-saturation constant for swimming speed ~ distance between current position and target position
REAL, parameter  :: kd = 100.

! Half-saturation constant for swimming speed ~ distance to surface or bottom
REAL, parameter  :: ks = 3.

IF (z > 0. .OR. z < bot_dep) then
    write(6,*) 'ZOO current z = ', z
    write(6,*) 'Bottom depth = ', bot_dep
    STOP "z position incorrect!"
ENDIF

! Calculate current light at z
PAR_ = PAR_0*exp(z*kext)

! Calculate desired location
if (PAR_0 > 0.) then
   z_des = log(isolume/PAR_0)/kext
   z_des = min(z_des, Z_MAX)
else
   z_des = Z_MAX
endif

! Distance between current location and desired location
diff = abs(z-z_des)

! swimming speed increases with the distance between current location and desired location
diff = diff/(diff+kd)  

! Distance from surface or bottom
! Retrieve bottom depth for current particle

dis  = min(-z, z-bot_dep)

if (dis < 0.) then
   write(6,*) 'Current vertical position: ', z
   write(6,*) 'Current Bottom depth: ', bot_dep
   stop "The vertical position is incorrect!"
else
   dis  = dis/(dis + ks)
endif


IF (z .GE. z_des) THEN
   swimdirection = -1. ! swim downward
ELSE
   swimdirection = 1.  ! swim upward
ENDIF
swim_speed = diff*dis*swimdirection
END FUNCTION swim_speed
!
!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
END MODULE LAGRANGE_MOD
