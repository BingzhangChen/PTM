PROGRAM IBM_2D
USE LAGRANGE_MOD
USE Read_u
IMPLICIT NONE

! Kv along X
real     :: Kvx(0:NX, NZ)       =    0d0

! Scratch indices:
integer  :: i,  k,  isave
real     :: start, finish, cday, bom


CALL cpu_time(start)

CALL RANDOM_SEED()

CALL initialize_grid

! Current time
ct    = 0
step0 = 1

! Initialize u, w, kap
CALL get_u(step0,   u1, w1, kap1) 
CALL get_u(step0+1, u2, w2, kap2) 

u   = u1
w   = w1
kap = kap1

! Current time of the day
ctd = mod(ct, d_per_s)

! Surface PAR
SPAR= surface_par(ctd)


CALL Init_Particles
CALL Create_particlefile
CALL Create_Flowfile

! Write into initial condition
CALL write_particlefile(1, 0.0)
CALL     write_Flowfile(1, 0.0)

! Check whether simulated time period exceeds offline SUNTANS physics data
IF (dt*Nstep > dt0*Nsteps) STOP "Simulated time period exceeds offline SUNTANS physics data!"

! Start main loop:

! Indeces of records for saving results:
isave = 1
Do i = 1, Nstep

   ! Update current time (in seconds)
   ct = ct + dt 

   ! Current day
   cday = dble(ct)/dble(d_per_s)

   ! Interpolate u,  w, and kap
   ! Determine whether needs to read new data
   IF (ct > step0*period) THEN

      step0= step0 + 1
   
      ! Update u1 and u2
      u1   = u2
      w1   = w2
      kap1 = kap2
   
      ! Read external data of the next time step
      CALL get_u(step0, u2, w2, kap2) 
   END IF

   CALL interp_u(ct, step0) 

   ! Current time of the day
   ctd = MOD(ct, d_per_s)
   
   ! Surface PAR
   SPAR= surface_par(ctd)

   ! Calculate zooplankton vertical swimming speed
   DO k = 1, TNPAR
      ! Find corresponding bottom depth
      bom     = Bot_dep(Zp(k)%ix, 2)
      Zp(k)%w = swim_speed(SPAR, Zp(k)%rz, Zp(k)%iso, bom) * Zp(k)%wm
   ENDDO

   ! Transport the particles
   CALL lagrange2D(nz,nx,Z_w,X_w,mask,kap,kvx,w,u, TNPAR, Zp)
   
   ! save data (positions of each particle) every hour
   IF (MOD(i,nsave) == 0) THEN

     isave = isave+1
     CALL write_particlefile(isave, cday)
     CALL     write_Flowfile(isave, cday)

   ENDIF
ENDDO

CALL CPU_TIME(finish)
PRINT '("Time = ",f8.3," hours.")', (finish-start)/3600.0 

END PROGRAM
