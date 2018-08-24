MODULE NETCDF_IO
USE NETCDF
USE LAGRANGE_MOD, only : NG, N_PAR, Zp
USE Grid_2D,      only : NX, NZ, X_r, Z_r, X_w, Z_w, MASK
USE Grid_2D,      only : u, w, kap
IMPLICIT NONE

! settings in offline SUNTANS data
integer, parameter :: dt0    = 10        !Time step in seconds in SUNTANS model
integer, parameter :: Nsteps = 356400
integer, parameter :: Ntout  = 360


! Seconds between each time step in offline SUNTANS data
integer, parameter :: period = Ntout*dt0

integer            :: step0  !Current step in offline SUNTANS data

character (len=20) :: ufile  = 'u.nc'
character (len=20) :: wfile  = 'w.nc'
character (len=20) :: kapfile= 'kappat.nc'

character (len=20) :: Particle_FNAME = 'Particles.nc'
character (len=20) ::     Flow_FNAME = 'Flow.nc'

!Number of Dimensions for particles and u
integer,            parameter  ::     NDIMS =  3
character (len=5),  parameter  ::    X_NAME = 'X'
character (len=5),  parameter  ::    Z_NAME = 'Z'
character (len=5),  parameter  ::   Xr_NAME = 'X_r'
character (len=5),  parameter  ::   Xw_NAME = 'X_w'
character (len=5),  parameter  ::   Zr_NAME = 'Z_r'
character (len=5),  parameter  ::   Zw_NAME = 'Z_w'
character (len=5),  parameter  :: MASK_NAME = 'Mask'
character (len=5),  parameter  ::  DAY_NAME = 'Day'
character (len=5),  parameter  ::    u_NAME = 'u'
character (len=5),  parameter  ::    w_NAME = 'w'
character (len=5),  parameter  ::  kap_NAME = 'kap'
character (len=10), parameter  ::     UNITS = 'units'

! Variable name for DVM vertical swimming speed
character (len=10), parameter ::   Speed_NAME = 'DVMSpd' 


CONTAINS

SUBROUTINE Create_particlefile

IMPLICIT NONE
character (len=10), parameter ::     GRP_NAME = 'Group'
character (len=10), parameter ::     PAT_NAME = 'ParticleID'
character (len=10), parameter ::     REC_NAME = 'Time'
character (len=10), parameter ::    UNIT_dist = 'm'
character (len=8)             ::    date

integer :: ncid, rec_dimid, GRP_dimid, PAT_dimid
integer :: X_varid, Z_varid, DAY_varid
integer :: u_varid, w_varid, kap_varid, spd_varid
integer :: dimids(NDIMS)

! Create particle file (nc file)
CALL check(nf90_create(Particle_FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
CALL check(nf90_def_dim(ncid, GRP_NAME, NG, GRP_dimid) )
CALL check(nf90_def_dim(ncid, PAT_NAME, N_PAR, PAT_dimid) )
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

! Define the variables of particle group and ID. 
! Ordinarily we would need to provide
! an array of dimension IDs for each variable's dimensions, but
! since group and particle variables only have one dimension, we can
! simply provide the address of that dimension ID
CALL check(nf90_def_var(ncid, DAY_NAME, NF90_REAL,rec_dimid, DAY_varid) )

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.
dimids = (/ GRP_dimid, PAT_dimid, rec_dimid /)

! Define the netCDF variables for the position data.
CALL check( nf90_def_var(ncid, X_NAME, NF90_REAL, dimids, X_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_REAL, dimids, Z_varid) )

CALL check( nf90_def_var(ncid, u_NAME, NF90_REAL, dimids, u_varid) )
CALL check( nf90_def_var(ncid, w_NAME, NF90_REAL, dimids, w_varid) )
CALL check( nf90_def_var(ncid, kap_NAME,NF90_REAL,dimids, kap_varid) )
CALL check( nf90_def_var(ncid, Speed_NAME,NF90_REAL,dimids, spd_varid) )

! Assign units attributes to the netCDF variables.
! Get the date of the nc file
CALL date_and_time(DATE=date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))

CALL check( nf90_put_att(ncid, X_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, Z_varid,   UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, u_varid,   UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, w_varid,   UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, spd_varid, UNITS, 'm s-1'))
CALL check( nf90_put_att(ncid, kap_varid, UNITS, 'm^2 s-1'))
CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))

! End define mode.
CALL check( nf90_enddef(ncid) )
  
CALL check( nf90_close(ncid) )
RETURN
END SUBROUTINE Create_particlefile

  
SUBROUTINE write_particlefile(rec, day)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: rec  !The time index to be written
REAL,    INTENT(IN)  :: day 
INTEGER              :: start(NDIMS)
integer :: ncid, X_varid, Z_varid, DAY_varid
integer :: u_varid, w_varid, kap_varid, spd_varid
integer :: iz, i, j

! These settings tell netcdf to write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(Particle_FNAME, NF90_WRITE, ncid))

! get variable IDs
CALL check(NF90_INQ_VARID(ncid, X_NAME,     X_varid))    
CALL check(NF90_INQ_VARID(ncid, Z_NAME,     Z_varid))     
CALL check(NF90_INQ_VARID(ncid, u_NAME,     u_varid))      
CALL check(NF90_INQ_VARID(ncid, w_NAME,     w_varid))      
CALL check(NF90_INQ_VARID(ncid, Speed_NAME, spd_varid))      
CALL check(NF90_INQ_VARID(ncid, kap_NAME,   kap_varid))      
CALL check(NF90_INQ_VARID(ncid, DAY_NAME,   DAY_varid))

CALL check(NF90_PUT_VAR(ncid, DAY_varid, REAL(day), start=(/rec/)))

iz = 0  !Index for Zp
DO i = 1, NG
  DO j = 1, N_PAR
     iz    = iz + 1
     start = (/i, j, rec/)
     CALL check(NF90_PUT_VAR(ncid, X_varid, REAL(Zp(iz)%rx), start=start))
     CALL check(NF90_PUT_VAR(ncid, Z_varid, REAL(Zp(iz)%rz), start=start))
     CALL check(NF90_PUT_VAR(ncid, u_varid, REAL(Zp(iz)%u), start=start))
     CALL check(NF90_PUT_VAR(ncid, w_varid, REAL(Zp(iz)%Uw), start=start))
     CALL check(NF90_PUT_VAR(ncid, spd_varid,REAL(Zp(iz)%w), start=start))
     CALL check(NF90_PUT_VAR(ncid, kap_varid,REAL(Zp(iz)%kap), start=start))
  ENDDO
ENDDO

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))
return
END SUBROUTINE write_particlefile

SUBROUTINE Create_Flowfile

IMPLICIT NONE
CHARACTER (LEN=10), PARAMETER ::     REC_NAME = 'Time'
CHARACTER (LEN=10), PARAMETER ::    UNIT_dist = 'm'

INTEGER :: ncid, rec_dimid, NX_dimid, NZ_dimid, NX1_dimid, NZ1_dimid
INTEGER :: Xr_varid, Zr_varid, Xw_varid, Zw_varid,DAY_varid, u_varid, w_varid, kap_varid
INTEGER :: MASK_varid
INTEGER :: dimids1(NDIMS), dimids2(NDIMS)

! Create flow file (nc file)
CALL check(nf90_create(Flow_FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed.
CALL check(nf90_def_dim(ncid, Xr_NAME,  NX,    NX_dimid))
CALL check(nf90_def_dim(ncid, Xw_NAME,  NX+1, NX1_dimid))
CALL check(nf90_def_dim(ncid, Zr_NAME,  NZ,    NZ_dimid))
CALL check(nf90_def_dim(ncid, Zw_NAME,  NZ+1, NZ1_dimid))
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

! Define the variables of particle group and ID. 
! Ordinarily we would need to provide
! an array of dimension IDs for each variable's dimensions, but
! since group and particle variables only have one dimension, we can
! simply provide the address of that dimension ID
CALL check(nf90_def_var(ncid, Xr_NAME,  NF90_REAL,  NX_dimid, Xr_varid) )
CALL check(nf90_def_var(ncid, Xw_NAME,  NF90_REAL, NX1_dimid, Xw_varid) )
CALL check(nf90_def_var(ncid, Zr_NAME,  NF90_REAL,  NZ_dimid, Zr_varid) )
CALL check(nf90_def_var(ncid, Zw_NAME,  NF90_REAL, NZ1_dimid, Zw_varid) )
CALL check(nf90_def_var(ncid, DAY_NAME, NF90_REAL, rec_dimid, DAY_varid) )

CALL check(nf90_def_var(ncid, MASK_NAME,NF90_INT, [NX_dimid, NZ_dimid], MASK_varid) )

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.
dimids1 = (/ NX1_dimid,  NZ_dimid, rec_dimid /)
dimids2 = (/  NX_dimid, NZ1_dimid, rec_dimid /)

! Define the netCDF variables for u, w, and kappa
CALL check( nf90_def_var(ncid, u_NAME,   NF90_REAL, dimids1, u_varid) )
CALL check( nf90_def_var(ncid, w_NAME,   NF90_REAL, dimids2, w_varid) )
CALL check( nf90_def_var(ncid, kap_NAME, NF90_REAL, dimids2, kap_varid) )

! Assign units attributes to the netCDF variables.
CALL check( nf90_put_att(ncid,   Xr_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid,   Xw_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid,   Zr_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid,   Zw_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))
CALL check( nf90_put_att(ncid,   u_varid, UNITS, 'm/s'))
CALL check( nf90_put_att(ncid,   w_varid, UNITS, 'm/s'))
CALL check( nf90_put_att(ncid, kap_varid, UNITS, 'm^2/s'))

! End define mode.
CALL check( nf90_enddef(ncid) )
  
! Write the X and Z coordinate data
CALL check( nf90_put_var(ncid, Xr_varid, X_r) )
CALL check( nf90_put_var(ncid, Xw_varid, X_w) )
CALL check( nf90_put_var(ncid, Zr_varid, Z_r) )
CALL check( nf90_put_var(ncid, Zw_varid, Z_w) )
CALL check( nf90_put_var(ncid, MASK_varid, MASK) )
CALL check( nf90_close(ncid) )
return

END SUBROUTINE Create_Flowfile

SUBROUTINE write_Flowfile(rec, day)
implicit none
integer, intent(in)  :: rec  !The time index to be written
real,    intent(in)  :: day 
integer              :: i, j
integer :: ncid, u_varid, w_varid, kap_varid, DAY_varid
real(4), allocatable :: cff(:,:)

! write one timestep of data.

!Open the nc file for writing
CALL check(NF90_OPEN(Flow_FNAME, NF90_WRITE, ncid))

CALL check(NF90_INQ_VARID(ncid, u_NAME, u_varid))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, w_NAME, w_varid))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, kap_NAME, kap_varid))      ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DAY_NAME, DAY_varid))  ! get variable IDs

CALL check(NF90_PUT_VAR(ncid, DAY_varid, real(day), start = (/rec/)))
allocate(cff(NX+1,NZ))
cff(:,:)=0.
do i = 0, NX
   do j = 1, NZ
      cff(i+1,j)=real(u(i,j))
   enddo
enddo
CALL check(NF90_PUT_VAR(ncid, u_varid, cff,start=[1,1,rec],&
                                           count=[NX+1,NZ,1]))
DEALLOCATE(cff)
ALLOCATE(cff(NX,NZ+1))
cff(:,:)=0.
do i = 1, NX
   do j = 0, NZ
      cff(i,j+1)=real(w(i,j))
   enddo
enddo

CALL check(NF90_PUT_VAR(ncid, w_varid, cff,start=[1,1,rec],&
                                           count=[NX,NZ+1,1]))

DO i = 1, NX
   do j = 0, NZ
      cff(i,j+1)=real(kap(i,j))
   enddo
ENDDO

CALL check(NF90_PUT_VAR(ncid,kap_varid, cff,start=[1,1,rec],&
                                            count=[NX,NZ+1,1]))

DEALLOCATE(cff)

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))
return
end subroutine write_Flowfile

SUBROUTINE NCREAD_3D(fname, vname, NX, NZ, step, dat)
IMPLICIT NONE
! This is the name of the data file we will read. 
character (len=*), intent(in) :: fname, vname

! We are reading 2D data, a 6 x 12 grid. 
integer, intent(in)  :: NX, NZ, step
real,    intent(out) :: dat(NX, NZ)

! This will be the netCDF ID for the file and data variable.
integer :: ncid, varid

! Loop indexes, and error handling.
integer :: start(3), count(3)

if (step > Nsteps/Ntout + 1) stop "You are retrieving data beyond the time range!"

! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
! the file.
CALL check(nf90_open(FNAME, NF90_NOWRITE, NCID) )

! Get the varid of the data variable, based on its name.
CALL check( nf90_inq_varid(ncid, vname, varid) )

! Read the data.
! step is the NO. of the record
start = [1, 1, step]
count = [NX,NZ,1]

call check( nf90_get_var(ncid, varid, dat, start = start, &
                                           count = count) )

! Close the file, freeing all resources.
call check( nf90_close(ncid) )
return

END SUBROUTINE NCREAD_3D

subroutine check(status)
integer, intent (in) :: status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check
END MODULE
