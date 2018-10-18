MODULE LAGRANGE_MOD
USE PARAMS_MOD
USE GRID_3D
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
INTEGER, PARAMETER :: nsave     =    h_per_s/int(dt)  ! Save every hour

INTEGER, PARAMETER :: NSTEP     =    D_PER_S*NDAYS/INT(DT)

! Current surface PAR in terms of diel cycle
REAL               :: SPAR      =    PAR0  

TYPE Particle

    ! the mth group
    INTEGER :: m

    ! the nth particle in the mth group
    INTEGER :: n

    ! Grid indices for particles
    INTEGER :: ix, iy, iz
    
    ! coordinates for particles
    REAL    :: rx,ry,rz
    
    ! associated u,w,and kap (vertical diffusivity)
    REAL    :: up,vp,wp,kapp

    ! swimming speed (m/s)
    REAL    :: w = 0.

    ! Maximal swimming speed (m/s)
    REAL    :: wm= 3.5D-2  ! m/s (Batchelder et al. PIO 2002)

    ! Prefered optimal light (isolume)
    REAL    :: iso= 1D2 ! W m-2
END TYPE Particle

! Positions of particles
TYPE (Particle)  :: Zp(TNPAR)

CONTAINS

SUBROUTINE INIT_PARTICLES
IMPLICIT NONE
INTEGER :: i,j,k,n
REAL    :: rnd

!Minimal and maximal X positions
REAL,    parameter :: x0= -5D4,  x1= -30D3
REAL,    parameter :: y0= -150D3,y1= -13D4

!Minimal and maximal depths
REAL,    parameter :: z0=-200.
REAL               :: bom ! A scratch variable for bottom depth
REAL               :: z2  ! A scratch z variable


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

         SELECT CASE(initial_condition)

         CASE(uniform)

           CALL RANDOM_NUMBER(rnd)

           !Randomly distributed in one grid
           !Different groups start from the same position

           Zp(k)%rx= rnd * (x1-x0) + x0

           CALL RANDOM_NUMBER(rnd)
           Zp(k)%ry= rnd * (y1-y0) + y0

          CASE(point)

           Zp(k)%rx= -37D3 
           Zp(k)%ry= -146.5D3 

         CASE DEFAULT
           STOP "INITIAL OPTIONS INCORRECT!"
         END SELECT

         ! Compute ix
         DO n = 1, NX
            IF (X_w(n) > Zp(k)%rx .AND. X_w(n-1) <= Zp(k)%rx) THEN
               Zp(k)%ix = n
               EXIT
            ENDIF
         ENDDO

         ! Compute iy
         DO n = 1, NY
            IF (Y_w(n) > Zp(k)%ry .AND. Y_w(n-1) <= Zp(k)%ry) THEN
               Zp(k)%iy = n
               EXIT
            ENDIF
         ENDDO

         ! Vertically, the particles uniformly randomly distributed 
         ! Bottom depth
         bom = BOT(Zp(k)%ix, Zp(k)%iy)

         CALL RANDOM_NUMBER(rnd)

         ! Vertical position must be above bottom depth
         ! Find corresponding bottom depth
         z2      = max(z0, bom)
         Zp(k)%rz= rnd * (Z_MAX-z2) + z2

         ! Find iz
         DO n=1, NZ
            IF (Z_w(n-1) < Zp(k)%rz .AND. Z_w(n) >= Zp(k)%rz) THEN
               Zp(k)%iz = n
               EXIT
            ENDIF
         END DO

         !Check if in a invalid grid
         DO WHILE (mask(Zp(k)%ix, Zp(k)%iy, Zp(k)%iz) == 0)
            IF (mask(Zp(k)%ix, Zp(k)%iy, NZ) == 1) THEN
               ! The particle is too deep
               Zp(k)%rz=Zp(k)%rz+.1
               DO n=Zp(k)%iz,NZ
                  if (Z_w(n-1) < Zp(k)%rz .and. Z_w(n) > Zp(k)%rz) then
                     Zp(k)%iz=n
                     EXIT
                  endif
               END DO
            ELSE
               ! The particle locates on the land
               Zp(k)%ry=Zp(k)%ry-.1
               DO n=Zp(k)%ix,NX
                  if (X_w(n-1) < Zp(k)%rx .and. X_w(n) > Zp(k)%rx) then
                     Zp(k)%ix=n
                     EXIT
                  endif
               END DO
            ENDIF
         ENDDO

      ELSE
         ! Track the same particle in the 1st group
         Zp(k)%iz= Zp(j)%iz
         Zp(k)%rz= Zp(j)%rz
         Zp(k)%ix= Zp(j)%ix
         Zp(k)%iy= Zp(j)%iy
         Zp(k)%rx= Zp(j)%rx
         Zp(k)%ry= Zp(j)%ry
      ENDIF

      ! Calculate zooplankton vertical swimming speed
      Zp(k)%w = swim_speed(SPAR,Zp(k)) * Zp(k)%wm

      k = k + 1
   ENDDO
ENDDO
RETURN
END SUBROUTINE Init_Particles

! !ROUTINE: Lagrangian particle random walk at 3D space (Horizontal and vertical)
SUBROUTINE lagrange3D
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
INTEGER              :: n,k
REAL, ALLOCATABLE    :: vel(:), nu(:) 
INTEGER, ALLOCATABLE :: mask_(:)
!EOP
!-----------------------------------------------------------------------
!BOC

! Save u, w, and kappat

DO n = 1, TNPAR
   Zp(n)%up  =   u(Zp(n)%ix, Zp(n)%iy, Zp(n)%iz)

   Zp(n)%vp  =   v(Zp(n)%ix, Zp(n)%iy, Zp(n)%iz)
   Zp(n)%wp  =   w(Zp(n)%ix, Zp(n)%iy, Zp(n)%iz)
   Zp(n)%kapp= kap(Zp(n)%ix, Zp(n)%iy, Zp(n)%iz)

   !Check whether zoo is in a fake grid (mask ==0)
   IF (mask(Zp(n)%ix,Zp(n)%iy, Zp(n)%iz) .eq. 0) THEN
       write(6,*) 'ix = ', Zp(n)%ix
       write(6,*) 'iy = ', Zp(n)%iy
       write(6,*) 'iz = ', Zp(n)%iz
       write(6,*) 'Bottom depth = ', BOT(Zp(n)%ix,Zp(n)%iy)
       stop "The particle is in a fake grid!"
   ENDIF
ENDDO

! Vertical movement
ALLOCATE(vel(0:nz))
vel(:) = 0d0

ALLOCATE(mask_(nz))
mask_(:) = 0

ALLOCATE(nu(0:nz))
nu(:) = 0.

DO n = 1,TNPAR
   ! Vertical velocity including swimming/sinking
   do k = 1, (nz-1)
      vel(k) = Zp(n)%w + w(Zp(n)%ix,Zp(n)%iy, k)
       nu(k) = kap(Zp(n)%ix,Zp(n)%iy,k)
   enddo
   do k = 1, nz
      mask_(k) = mask(Zp(n)%ix, Zp(n)%iy, k)
   enddo

   CALL LAGRANGE(NZ, Z_w, mask_,nu,vel,Zp(n)%iz, Zp(n)%rz)
ENDDO

DEALLOCATE(vel)
DEALLOCATE(nu)
DEALLOCATE(mask_)

! Horizontal movement (u)
ALLOCATE(vel(0:nx))
vel(:) = 0d0
ALLOCATE(mask_(nx))
mask_(:)= 0
ALLOCATE(nu(0:nx))
nu(:)= 0.1

DO n = 1,TNPAR
   ! Horizontal velocity 
   do k = 1, (NX-1)
      vel(k) =  u(k,Zp(n)%iy,Zp(n)%iz)
   enddo
   do k = 1, NX
      mask_(k) = mask(k,Zp(n)%iy,Zp(n)%iz)
   enddo

   CALL LAGRANGE(NX,X_w,mask_,nu,vel, Zp(n)%ix, Zp(n)%rx)

ENDDO

DEALLOCATE(vel)
DEALLOCATE(nu)
DEALLOCATE(mask_)

! Horizontal movement (v)
ALLOCATE(vel(0:NY))
vel(:) = 0d0
ALLOCATE(mask_(NY))
mask_(:)= 0
ALLOCATE(nu(0:NY))
nu(:)= 0.1

DO n = 1, TNPAR
   ! v
   do k = 1, (NY-1)
      vel(k) =  u(Zp(n)%ix, k, Zp(n)%iz)
   enddo

   do k = 1, NY
      mask_(k) = mask(Zp(n)%ix,k,Zp(n)%iz)
   enddo

   CALL LAGRANGE(NY, Y_w, mask_, nu, vel, Zp(n)%iy, Zp(n)%ry)

ENDDO
RETURN
END SUBROUTINE LAGRANGE3D

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
REAL FUNCTION swim_speed(PAR_0,z)
IMPLICIT NONE
REAL,           INTENT(IN) :: PAR_0
TYPE(PARTICLE), INTENT(IN) :: Z
REAL             :: PAR_, swimdirection, z_des, diff, dis
REAL             :: bot_dep

! Half-saturation constant for swimming speed ~ distance between current position and target position
REAL, PARAMETER  :: kd = 100.

! Half-saturation constant for swimming speed ~ distance to surface or bottom
REAL, PARAMETER  :: ks = 3.
INTEGER          :: i

bot_dep=BOT(z%ix, z%iy)
IF (z%rz > 0. .OR. z%rz < bot_dep) then
    write(6,*) 'ZOO current ix = ', z%ix
    write(6,*) 'ZOO current iy = ', z%iy
    write(6,*) 'ZOO current iz = ', z%iz
    write(6,*) 'ZOO current rz = ', z%rz
    write(6,*) 'Bottom depth = ', bot_dep
    do i=1,NZ
       write(6,*) 'i   = ',i
       write(6,*) 'Z_w = ',Z_w(i)
       write(6,*) 'mask= ',mask(z%ix,z%iy,i)
    enddo

    STOP "z position incorrect!"
ENDIF

! Calculate current light at z
PAR_ = PAR_0*exp(z%rz*kext)

! Calculate desired location
if (PAR_0 > 0.) then
   z_des = log(z%iso/PAR_0)/kext
   z_des = min(z_des, Z_MAX)
else
   z_des = Z_MAX
endif

! Distance between current location and desired location
diff = abs(z%rz - z_des)

! swimming speed increases with the distance between current location and desired location
diff = diff/(diff+kd)  

! Distance from surface or bottom
! Retrieve bottom depth for current particle

dis  = min(-z%rz, z%rz - bot_dep)

if (dis < 0.) then
   write(6,*) 'Current vertical position: ', z
   write(6,*) 'Current Bottom depth: ', bot_dep
   stop "The vertical position is incorrect!"
else
   dis  = dis/(dis + ks)
endif

IF (z%rz .GE. z_des) THEN
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
