PROGRAM IBM_3D
USE LAGRANGE_MOD
USE Read_u
IMPLICIT NONE

! Scratch indices:
INTEGER  :: i, j, k, isave, LASTTIME
REAL     :: CDAY
REAL(4)  :: start, finish
REAL(4)  :: time_initial, t1, t2, t3, t4
REAL(4)  :: time_interp
REAL(4)  :: time_PTM, time_save

CALL CPU_TIME(START)

CALL RANDOM_SEED()

CALL initialize_grid

! Current time
ct    = 0
step0 = 1

! Initialize u, w, kap
CALL get_u(step0,   u1, v1, w1, kap1) 
CALL get_u(step0+1, u2, v2, w2, kap2) 

u   = u1
v   = v1
w   = w1
kap = kap1

! Current time of the day
CTD = MOD(CT, D_PER_S)

! Surface PAR
SPAR= surface_par(ctd)

IF (.NOT. RESTART) THEN
    CALL Init_Particles
    CALL Create_particlefile
    LASTTIME=0
ELSE
    CALL READ_RESTART(RESTARTNAME, LASTTIME)
ENDIF

CT=CT+LASTTIME
CALL Create_Flowfile

! Write into initial condition
CALL WRITE_PARTICLEFILE(1, 0.0)
CALL     WRITE_FLOWFILE(1, 0.0)

CALL CPU_TIME(time_initial)

PRINT '("Initialization takes ",f8.3," hours.")', (time_initial - start)/3600.0 

! Check whether simulated time period exceeds offline SUNTANS physics data
IF (dt*Nstep > dt0*Nsteps) STOP "Simulated time period exceeds offline SUNTANS physics data!"

! Start main loop:
TIME_INTERP = 0.
TIME_PTM    = 0.
TIME_SAVE   = 0.

! Indeces of records for saving results:
isave = 1
Do i = 1, Nstep

   ! Update current time (in seconds)
   ct = ct + dt 

   ! Current day
   CDAY = DBLE(CT)/DBLE(D_PER_S)

   ! Interpolate u,  w, and kap
   ! Determine whether needs to read new data
   CALL CPU_TIME(T1)
   IF (ct > step0*period) THEN

      step0= step0 + 1
   
      ! Update u1 and u2
      u1   = u2
      v1   = v2
      w1   = w2
      kap1 = kap2
   
      ! Read external data of the next time step
      CALL get_u(step0+1, u2, v2,  w2, kap2) 
   END IF

   CALL interp_u(ct, step0) 

   CALL CPU_TIME(T2)

   TIME_INTERP = TIME_INTERP + T2 - T1

   ! Current time of the day
   ctd = MOD(ct, d_per_s)
   
   ! Surface PAR
   SPAR= surface_par(ctd)

   ! Calculate zooplankton vertical swimming speed
   DO k = 1, TNPAR
      Zp(k)%w = swim_speed(SPAR, Zp(k)) * Zp(k)%wm
   ENDDO

   ! Transport the particles
   CALL lagrange3D

   CALL CPU_TIME(T3)
   TIME_PTM = TIME_PTM + T3 - T2

   ! save data (positions of each particle) every hour
   IF (MOD(i,nsave) == 0) THEN
     write(6,1804) 'Save results at day', cday
     isave = isave+1
     CALL write_particlefile(isave, cday)
     CALL     write_Flowfile(isave, cday)
   ENDIF

   CALL CPU_TIME(T4)
   TIME_save = TIME_save + T4 - T3
ENDDO

CALL CPU_TIME(finish)
PRINT '("Total time = ",f8.3," hours.")', (finish-start)/3600.0 
PRINT '("Physics interpolation time = ",f8.3," hours.")', time_interp/3600.0 
PRINT '("Particle tracking time = ",f8.3," hours.")', time_PTM/3600.0 
PRINT '("Data saving time = ",f8.3," hours.")', time_save/3600.0 
1804 FORMAT(A19,1X, F12.2)
END PROGRAM
