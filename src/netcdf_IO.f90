MODULE NETCDF_IO
USE NETCDF
USE LAGRANGE_MOD, only : NG, N_PAR, Zp, d_per_s
USE Grid_3D,      only : NX, NY, NZ, X_r, Y_r,  Z_r, X_w, Y_w, Z_w, MASK
USE Grid_3D,      only : u, v, w, kap, BOT
IMPLICIT NONE

! settings in offline SUNTANS data
INTEGER, PARAMETER :: dt0    = 10        !Time step in seconds in SUNTANS model
INTEGER, PARAMETER :: Nsteps = 86400
INTEGER, PARAMETER :: Ntout  = 360

! Seconds between each time step in offline SUNTANS data
integer, parameter :: period = Ntout*dt0

integer            :: step0  !Current step in offline SUNTANS data

character (len=20) :: ufile  = 'u.nc'
character (len=20) :: kapfile= 'kappat.nc'

character (len=20) :: Particle_FNAME = 'Particles.nc'
character (len=20) ::     Flow_FNAME = 'Flow.nc'

!Number of Dimensions for particles and u
INTEGER,            PARAMETER  ::     NDIMS =  4
CHARACTER (LEN=5),  PARAMETER  ::   IX_NAME = 'IX'
CHARACTER (LEN=5),  PARAMETER  ::   IY_NAME = 'IY'
CHARACTER (LEN=5),  PARAMETER  ::   IZ_NAME = 'IZ'
CHARACTER (LEN=5),  PARAMETER  ::    X_NAME = 'X'
CHARACTER (LEN=5),  PARAMETER  ::    Y_NAME = 'Y'
CHARACTER (LEN=5),  PARAMETER  ::    Z_NAME = 'Z'
CHARACTER (LEN=5),  PARAMETER  ::   Xr_NAME = 'X_r'
CHARACTER (LEN=5),  PARAMETER  ::   Xw_NAME = 'X_w'
CHARACTER (LEN=5),  PARAMETER  ::   Yr_NAME = 'Y_r'
CHARACTER (LEN=5),  PARAMETER  ::   Yw_NAME = 'Y_w'
CHARACTER (LEN=5),  PARAMETER  ::   Zr_NAME = 'Z_r'
CHARACTER (LEN=5),  PARAMETER  ::   Zw_NAME = 'Z_w'
CHARACTER (LEN=5),  PARAMETER  :: MASK_NAME = 'Mask'
CHARACTER (LEN=10), PARAMETER  ::  BOT_NAME = 'Bathymetry'
CHARACTER (LEN=5),  PARAMETER  ::  DAY_NAME = 'Day'
CHARACTER (LEN=5),  PARAMETER  ::    u_NAME = 'u'
CHARACTER (LEN=5),  PARAMETER  ::    v_NAME = 'v'
CHARACTER (LEN=5),  PARAMETER  ::    w_NAME = 'w'
CHARACTER (LEN=5),  PARAMETER  ::  KAP_NAME = 'kap'
CHARACTER (LEN=5),  PARAMETER  ::     UNITS = 'units'

! Variable name for DVM vertical swimming speed
CHARACTER (LEN=10), PARAMETER ::   SPEED_NAME = 'DVMSPD' 

LOGICAL :: FLOWFILEEXISTS = .FALSE.

CONTAINS

!Obtain current particles position
SUBROUTINE READ_RESTART(RESTARTFILE, TIME_SEC)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: RESTARTFILE
INTEGER,          INTENT(OUT):: TIME_SEC  ! CURRENT TIME IN SECONDS
INTEGER :: iz, i, j
INTEGER :: ncid, X_varid,Y_varid, Z_varid, DAY_varid
INTEGER :: IX_varid, IY_varid, IZ_varid
INTEGER :: u_varid,  v_varid, w_varid, kap_varid, spd_varid
INTEGER :: TimeDimID, NRECORDS, start(3)
REAL    :: day
LOGICAL :: file_exists
CHARACTER(len=10) :: TimeDimName

! Check whether the restart file exists
INQUIRE(FILE=restartfile, EXIST=file_exists)

if (file_exists) then
   CALL CHECK(NF90_OPEN(RESTARTFILE, NF90_NOWRITE, NCID) )
else
   stop "Restart file does not exist!"
endif

! Inquire the number of time records
! What is the name of the unlimited dimension, how many records are there?
call check(nf90_inq_dimid(ncid, "Time", TimeDimID))
call check(nf90_inquire_dimension(ncid, TimeDimID, &
               name = TimeDimName, len = NRECORDS))

! get variable IDs
CALL check(NF90_INQ_VARID(ncid, X_NAME,     X_varid))    
CALL check(NF90_INQ_VARID(ncid, Y_NAME,     Y_varid))    
CALL check(NF90_INQ_VARID(ncid, Z_NAME,     Z_varid))     
CALL check(NF90_INQ_VARID(ncid, IX_NAME,    IX_varid))    
CALL check(NF90_INQ_VARID(ncid, IY_NAME,    IY_varid))    
CALL check(NF90_INQ_VARID(ncid, IZ_NAME,    IZ_varid))    
CALL check(NF90_INQ_VARID(ncid, u_NAME,     u_varid))      
CALL check(NF90_INQ_VARID(ncid, v_NAME,     v_varid))      
CALL check(NF90_INQ_VARID(ncid, w_NAME,     w_varid))      
CALL check(NF90_INQ_VARID(ncid, Speed_NAME, spd_varid))      
CALL check(NF90_INQ_VARID(ncid, kap_NAME,   kap_varid))      
CALL check(NF90_INQ_VARID(ncid, DAY_NAME,   DAY_varid))

! Read the final record from the restart file.
! step is the NO. of the record

! Read running time
CALL CHECK( NF90_GET_VAR(ncid, DAY_varid, day, start = (/Nrecords/)))
TIME_SEC=int(day*dble(d_per_s))

iz = 0  !Index for Zp
DO i = 1, NG
  DO j = 1, N_PAR
     iz    = iz + 1
     start = (/i, j, NRECORDS/)
     CALL check(NF90_GET_VAR(ncid, IX_varid, Zp(iz)%ix,  start=start))
     CALL check(NF90_GET_VAR(ncid, IY_varid, Zp(iz)%iy,  start=start))
     CALL check(NF90_GET_VAR(ncid, IZ_varid, Zp(iz)%iz,  start=start))
     CALL check(NF90_GET_VAR(ncid, X_varid,  Zp(iz)%rx,  start=start))
     CALL check(NF90_GET_VAR(ncid, Y_varid,  Zp(iz)%ry,  start=start))
     CALL check(NF90_GET_VAR(ncid, Z_varid,  Zp(iz)%rz,  start=start))
     CALL check(NF90_GET_VAR(ncid, u_varid,  Zp(iz)%up,  start=start))
     CALL check(NF90_GET_VAR(ncid, v_varid,  Zp(iz)%vp,  start=start))
     CALL check(NF90_GET_VAR(ncid, w_varid,  Zp(iz)%wp,  start=start))
     CALL check(NF90_GET_VAR(ncid, spd_varid,Zp(iz)%w,   start=start))
     CALL check(NF90_GET_VAR(ncid, kap_varid,Zp(iz)%kapp,start=start))
  ENDDO
ENDDO

CALL CHECK(NF90_CLOSE(NCID) )

END SUBROUTINE READ_RESTART

SUBROUTINE CREATE_PARTICLEFILE
IMPLICIT NONE
CHARACTER (len=10), parameter ::     GRP_NAME = 'Group'
CHARACTER (len=10), parameter ::     PAT_NAME = 'ParticleID'
CHARACTER (len=10), parameter ::     REC_NAME = 'Time'
CHARACTER (len=10), parameter ::    UNIT_dist = 'm'
CHARACTER (len=8)             ::    date

INTEGER :: ncid, rec_dimid, GRP_dimid, PAT_dimid
INTEGER :: X_varid, Y_varid, Z_varid, DAY_varid
INTEGER :: IX_varid,IY_varid, IZ_varid
INTEGER :: u_varid, v_varid, w_varid, kap_varid, spd_varid
INTEGER :: dimids(3)

! Create particle file (nc file)
CALL CHECK(NF90_CREATE(PARTICLE_FNAME, NF90_CLOBBER, NCID) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
CALL CHECK(NF90_DEF_DIM(NCID, GRP_NAME, NG,    GRP_DIMID) )
CALL CHECK(NF90_DEF_DIM(NCID, PAT_NAME, N_PAR, PAT_DIMID) )
CALL CHECK(NF90_DEF_DIM(NCID, REC_NAME, NF90_UNLIMITED, REC_DIMID) )

! Define the variables of particle group and ID. 
! Ordinarily we would need to provide
! an array of dimension IDs for each variable's dimensions, but
! since group and particle variables only have one dimension, we can
! simply provide the address of that dimension ID
CALL CHECK(NF90_DEF_VAR(NCID, DAY_NAME, NF90_REAL, REC_DIMID, DAY_VARID) )

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.
dimids = (/ GRP_dimid, PAT_dimid, rec_dimid /)

! Define the netCDF variables for the position data.
CALL check( nf90_def_var(ncid, IX_NAME,   NF90_INT,  dimids, IX_varid) )
CALL check( nf90_def_var(ncid, IY_NAME,   NF90_INT,  dimids, IY_varid) )
CALL check( nf90_def_var(ncid, IZ_NAME,   NF90_INT,  dimids, IZ_varid) )
CALL check( nf90_def_var(ncid, X_NAME,    NF90_REAL, dimids, X_varid)  )
CALL check( nf90_def_var(ncid, Y_NAME,    NF90_REAL, dimids, Y_varid)  )
CALL check( nf90_def_var(ncid, Z_NAME,    NF90_REAL, dimids, Z_varid)  )

CALL check( nf90_def_var(ncid, u_NAME,    NF90_REAL, dimids, u_varid)   )
CALL check( nf90_def_var(ncid, v_NAME,    NF90_REAL, dimids, v_varid)   )
CALL check( nf90_def_var(ncid, w_NAME,    NF90_REAL, dimids, w_varid)   )
CALL check( nf90_def_var(ncid, kap_NAME,  NF90_REAL, dimids, kap_varid) )
CALL check( nf90_def_var(ncid, Speed_NAME,NF90_REAL, dimids, spd_varid) )

! Assign units attributes to the netCDF variables.
! Get the date of the nc file
CALL DATE_AND_TIME(DATE=DATE)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))
CALL check( nf90_put_att(ncid, X_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, Y_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, Z_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, u_varid,   UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, v_varid,   UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, w_varid,   UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, spd_varid, UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, kap_varid, UNITS, 'm^2 s-1'))
CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))

! End define mode.
CALL check( nf90_enddef(ncid) )
CALL check( nf90_close(ncid) )
RETURN
END SUBROUTINE Create_particlefile
  
SUBROUTINE WRITE_PARTICLEFILE(REC, DAY)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: rec  !The time index to be written
REAL,    INTENT(IN)  :: day 
INTEGER              :: start(3)
integer :: ncid, X_varid,Y_varid, Z_varid, DAY_varid
integer :: IX_varid, IY_varid, IZ_varid
integer :: u_varid,  v_varid, w_varid, kap_varid, spd_varid
integer :: iz, i, j

! These settings tell netcdf to write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(Particle_FNAME, NF90_WRITE, ncid))

! get variable IDs
CALL check(NF90_INQ_VARID(ncid, X_NAME,     X_varid))    
CALL check(NF90_INQ_VARID(ncid, Y_NAME,     Y_varid))    
CALL check(NF90_INQ_VARID(ncid, Z_NAME,     Z_varid))     
CALL check(NF90_INQ_VARID(ncid, IX_NAME,    IX_varid))    
CALL check(NF90_INQ_VARID(ncid, IY_NAME,    IY_varid))    
CALL check(NF90_INQ_VARID(ncid, IZ_NAME,    IZ_varid))    
CALL check(NF90_INQ_VARID(ncid, u_NAME,     u_varid))      
CALL check(NF90_INQ_VARID(ncid, v_NAME,     v_varid))      
CALL check(NF90_INQ_VARID(ncid, w_NAME,     w_varid))      
CALL check(NF90_INQ_VARID(ncid, Speed_NAME, spd_varid))      
CALL check(NF90_INQ_VARID(ncid, kap_NAME,   kap_varid))      
CALL check(NF90_INQ_VARID(ncid, DAY_NAME,   DAY_varid))

CALL check(NF90_PUT_VAR(ncid, DAY_varid, day, start=(/rec/)))

iz = 0  !Index for Zp
DO i = 1, NG
  DO j = 1, N_PAR
     iz    = iz + 1
     start = (/i, j, rec/)
     CALL check(NF90_PUT_VAR(ncid, IX_varid, Zp(iz)%ix,  start=start))
     CALL check(NF90_PUT_VAR(ncid, IY_varid, Zp(iz)%iy,  start=start))
     CALL check(NF90_PUT_VAR(ncid, IZ_varid, Zp(iz)%iz,  start=start))
     CALL check(NF90_PUT_VAR(ncid, X_varid,  Zp(iz)%rx,  start=start))
     CALL check(NF90_PUT_VAR(ncid, Y_varid,  Zp(iz)%ry,  start=start))
     CALL check(NF90_PUT_VAR(ncid, Z_varid,  Zp(iz)%rz,  start=start))
     CALL check(NF90_PUT_VAR(ncid, u_varid,  Zp(iz)%up,  start=start))
     CALL check(NF90_PUT_VAR(ncid, v_varid,  Zp(iz)%vp,  start=start))
     CALL check(NF90_PUT_VAR(ncid, w_varid,  Zp(iz)%wp,  start=start))
     CALL check(NF90_PUT_VAR(ncid, spd_varid,Zp(iz)%w,   start=start))
     CALL check(NF90_PUT_VAR(ncid, kap_varid,Zp(iz)%kapp,start=start))
  ENDDO
ENDDO

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))
return
END SUBROUTINE write_particlefile

SUBROUTINE CREATE_FLOWFILE
IMPLICIT NONE
CHARACTER (LEN=10), PARAMETER ::     REC_NAME = 'Time'
CHARACTER (LEN=10), PARAMETER ::    UNIT_dist = 'm'

INTEGER :: ncid, rec_dimid, NX_dimid, NY_dimid, NZ_dimid, NX1_dimid, NY1_dimid,NZ1_dimid
INTEGER :: Xr_varid,Yr_varid, Zr_varid, Xw_varid,Yw_varid, Zw_varid,DAY_varid
INTEGER :: u_varid,v_varid, w_varid, kap_varid
INTEGER :: MASK_varid, bot_varid
INTEGER :: i,j,k
INTEGER :: dimids1(NDIMS),  dimids3(NDIMS)

IF (FLOWFILEEXISTS) THEN
    RETURN
ELSE
    ! Create flow file (nc file)
    CALL CHECK(NF90_CREATE(FLOW_FNAME, NF90_CLOBBER, NCID) )
      
    ! Define the dimensions. The record dimension is defined to have
    ! unlimited length - it can grow as needed.
    CALL check(nf90_def_dim(ncid, Xr_NAME,  NX,    NX_dimid))
    CALL check(nf90_def_dim(ncid, Xw_NAME,  NX+1, NX1_dimid))
    CALL check(nf90_def_dim(ncid, Yr_NAME,  NY,    NY_dimid))
    CALL check(nf90_def_dim(ncid, Yw_NAME,  NY+1, NY1_dimid))
    CALL check(nf90_def_dim(ncid, Zr_NAME,  NZ,    NZ_dimid))
    CALL check(nf90_def_dim(ncid, Zw_NAME,  NZ+1, NZ1_dimid))
    CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    
    ! Define the variables of particle group and ID. 
    ! Ordinarily we would need to provide
    ! an array of dimension IDs for each variable's dimensions, but
    ! since group and particle variables only have one dimension, we can
    ! simply provide the address of that dimension ID
    CALL check(nf90_def_var(ncid, Xr_NAME,  NF90_real,  NX_dimid, Xr_varid) )
    CALL check(nf90_def_var(ncid, Xw_NAME,  NF90_real, NX1_dimid, Xw_varid) )
    CALL check(nf90_def_var(ncid, Yr_NAME,  NF90_real,  NY_dimid, Yr_varid) )
    CALL check(nf90_def_var(ncid, Yw_NAME,  NF90_real, NY1_dimid, Yw_varid) )
    CALL check(nf90_def_var(ncid, Zr_NAME,  NF90_real,  NZ_dimid, Zr_varid) )
    CALL check(nf90_def_var(ncid, Zw_NAME,  NF90_real, NZ1_dimid, Zw_varid) )
    CALL check(nf90_def_var(ncid, DAY_NAME, NF90_real, rec_dimid, DAY_varid) )
    
    CALL check(nf90_def_var(ncid, MASK_NAME,NF90_short, [NX_dimid, NY_dimid,NZ_dimid], MASK_varid) )
    CALL check(nf90_def_var(ncid,  BOT_NAME,NF90_real,[NX_dimid, NY_dimid], BOT_varid) )
    
    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids1 = (/  NX_dimid,  NY_dimid,  NZ_dimid, rec_dimid /)
    dimids3 = (/  NX_dimid,  NY_dimid, NZ1_dimid, rec_dimid /)
    
    ! Define the netCDF variables for u, w, and kappa
    CALL check( nf90_def_var(ncid, u_NAME,   NF90_real, dimids1,   u_varid) )
    CALL check( nf90_def_var(ncid, v_NAME,   NF90_real, dimids1,   v_varid) )
    CALL check( nf90_def_var(ncid, w_NAME,   NF90_real, dimids3,   w_varid) )
    CALL check( nf90_def_var(ncid, kap_NAME, NF90_real, dimids3, kap_varid) )
    
    ! Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid,   Xr_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid,   Yr_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid,   Xw_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid,   Yw_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid,   Zr_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid,   Zw_varid, UNITS, UNIT_dist))
    CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))
    CALL check( nf90_put_att(ncid,   u_varid, UNITS, 'm/s'))
    CALL check( nf90_put_att(ncid,   w_varid, UNITS, 'm/s'))
    CALL check( nf90_put_att(ncid, kap_varid, UNITS, 'm^2/s'))
    
    ! End define mode.
    CALL CHECK( NF90_ENDDEF(NCID) )
    
    ! Write the X and Z coordinate data
    CALL check( nf90_put_var(ncid, Xr_varid,   X_r) )
    CALL check( nf90_put_var(ncid, Yr_varid,   Y_r) )
    CALL check( nf90_put_var(ncid, Xw_varid,   X_w) )
    CALL check( nf90_put_var(ncid, Yw_varid,   Y_w) )
    CALL check( nf90_put_var(ncid, Zr_varid,   Z_r) )
    CALL check( nf90_put_var(ncid, Zw_varid,   Z_w) )
    CALL check( nf90_put_var(ncid, BOT_varid,  BOT))
    do i = 1, NX
       do j = 1, NY
          do k = 1, NZ
             CALL check(nf90_put_var(ncid, mask_varid, mask(i,j,k),start=[i,j,k]))
          enddo
       enddo
    enddo
    CALL CHECK( NF90_CLOSE(NCID) )
ENDIF
RETURN
END SUBROUTINE Create_Flowfile

SUBROUTINE WRITE_FLOWFILE(REC, DAY)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: rec  !The time index to be written
REAL,    INTENT(IN)  :: day 
INTEGER              :: i, j, k
INTEGER              :: ncid, u_varid,v_varid, w_varid, kap_varid, DAY_varid
REAL,    ALLOCATABLE :: cff(:,:,:)

IF (FLOWFILEEXISTS) RETURN

! write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(Flow_FNAME, NF90_WRITE, ncid))

CALL check(NF90_INQ_VARID(ncid,   u_NAME, U_VARID))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid,   v_NAME, V_VARID))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid,   w_NAME, W_VARID))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, kap_NAME, KAP_VARID))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DAY_NAME, DAY_VARID))  ! get variable IDs

CALL CHECK(NF90_PUT_VAR(NCID, DAY_VARID, DAY, START = (/REC/)))

ALLOCATE(cff(NX,NY,NZ))
!Interpolate to r points
cff(:,:,:)=0.
cff(1,:,:)=u(1,:,:)
cff(NX,:,:)=u(NX-1,:,:)
DO i = 2, NX-1
   do j = 1, NY
      do k = 1, NZ
         cff(i,j,k)=(u(i-1,j,k) + u(i,j,k))/2.
      enddo
   enddo
ENDDO

CALL CHECK(NF90_PUT_VAR(NCID, U_VARID, REAL(CFF),START=[1,1,1,REC],&
                                                 COUNT=[NX,NY,NZ,1]))

cff(:,:,:)=0.
cff(:,1,:)=v(i,1,k)
cff(:,NY,:)=v(i,NY-1,k)
do i = 1, NX
   do j = 2, NY-1
      do k = 1, NZ
         cff(i,j,k)=(v(i,j-1,k) + v(i,j,k))/2.
      enddo
   enddo
enddo

CALL CHECK(NF90_PUT_VAR(NCID, V_VARID, REAL(CFF),START=[1,1,1,REC],&
                                                 COUNT=[NX,NY,NZ,1]))


DEALLOCATE(cff)

ALLOCATE(cff(NX,NY,NZ+1))
cff(:,:,:)=0.
do i = 1, NX
   do j = 1, NY
      do k = 0, NZ
         cff(i,j,k+1)=w(i,j,k)
      enddo
   enddo
enddo
CALL CHECK(NF90_PUT_VAR(NCID, W_VARID, real(CFF), &
           START=[1,1,1,REC], COUNT=[NX,NY,NZ+1,1]))

DO i = 1, NX
   do j = 1, NY
      do k = 0, NZ
         cff(i,j,k+1)=kap(i,j,k)
      enddo
   enddo
ENDDO
CALL CHECK(NF90_PUT_VAR(NCID, kap_VARID, real(CFF), &
           START=[1,1,1,REC], COUNT=[NX,NY,NZ+1,1]))


DEALLOCATE(cff)

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL CHECK(NF90_CLOSE(NCID))
RETURN
END SUBROUTINE write_Flowfile

SUBROUTINE NCREAD_3D(fname, vname, N_X, N_Z, step, dat)
IMPLICIT NONE
! This is the name of the data file we will read. 
CHARACTER(LEN=*), INTENT(IN)  :: fname, vname
INTEGER,          INTENT(IN)  :: N_X, N_Z, step
REAL,             INTENT(OUT) :: dat(N_X, N_Z)

! This will be the netCDF ID for the file and data variable.
integer :: ncid, varid, i,j

real    :: val

if (step > Nsteps/Ntout + 1) stop "You are retrieving data beyond the time range!"

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.
CALL check(NF90_OPEN(FNAME, NF90_NOWRITE, NCID) )

! Get the varid of the data variable, based on its name.
CALL check( NF90_INQ_VARID(ncid, vname, varid) )

! Read the data.
! step is the NO. of the record
do i = 1, N_X
    do j = 1, N_Z
        CALL CHECK(NF90_GET_VAR(ncid, varid, val, start = (/i, j, step/)))
        dat(i,j) = val
    enddo
enddo

! Close the file, freeing all resources.
CALL CHECK( NF90_CLOSE(NCID) )
RETURN

END SUBROUTINE NCREAD_3D

SUBROUTINE CHECK(STATUS)
INTEGER, INTENT (IN) :: STATUS
IF(STATUS /= NF90_NOERR) THEN
  PRINT *, TRIM(NF90_STRERROR(STATUS))
  STOP "STOPPED"
END IF
END SUBROUTINE CHECK
END MODULE NETCDF_IO
