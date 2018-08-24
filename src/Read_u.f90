Module Read_u
USE NETCDF_IO
USE grid_2D
IMPLICIT NONE

CONTAINS

SUBROUTINE interp_u(t, step) 
implicit none

! Current model time
integer, intent(in)   :: t


! Current step in offline SUNTANS data
integer, intent(in)   :: step

integer               :: i,j
integer               :: t1, t2

! Calculate two time points of step and step + 1
t1 = (step-1)*period
t2 =  step   *period

if (t < t1 .or. t > t2) then
    stop "Current time NOT within the time interval! Something wrong!"
else
    do i = 0, NX
       do j = 1, NZ
          u(i,j) = (u1(i,j)*(t2 - t) + u2(i,j)*(t - t1))/(t2 - t1)
       enddo
    enddo
    do i = 1, NX
       do j = 0, NZ
            w(i,j)=  (w1(i,j)*(t2 - t) +   w2(i,j)*(t - t1))/(t2 - t1)
          kap(i,j)=(kap1(i,j)*(t2 - t) + kap2(i,j)*(t - t1))/(t2 - t1)
       enddo
    enddo
end if

return
END SUBROUTINE


SUBROUTINE get_u(rec, u0, w0, kap0) 
IMPLICIT NONE

! The record to be read from nc files
INTEGER, intent(in)   :: rec

REAL,    intent(out)  :: u0(0:NX,NZ), w0(NX,0:NZ), kap0(NX,0:NZ)

INTEGER               :: i,j

! Scratch matrix to store temporary data read from nc file
REAL                  :: cff(NX,NZ)= 0d0
REAL                  :: cff1(NX,NZ)= 0d0

! Filename for u and w
character(len=10)     :: ufile     = 'u.nc'
character(len=10)     :: wfile     = 'w.nc'
character(len=10)     :: kappafile = 'kappat.nc'
real,    parameter    :: NAN       = 9D4

! Read the external u,  w, kappa to set up initial condition

CALL NCREAD_3D(ufile, 'u', NX, NZ, rec, cff)

do j = 1, NZ
   do i = 1, NX
      ! Correct NA values
      if (cff(i,j) > NAN) then
          cff(i,j) = 0.
      endif
   enddo
enddo


! Interpolate cff ==> u
u0(0, :) = 0.
u0(NX,:) = 0.
do j = 1, NZ
   do i = 1, NX-1
      u0(i,NZ-j+1) = (cff(i,j)*dist(i+1) + cff(i+1,j)*dist(i)) &
                   / (dist(i) + dist(i+1))
   enddo
enddo

CALL NCREAD_3D(wfile, 'w', NX, NZ, 1, cff)

do j = 1, NZ
   do i = 1, NX
      ! Correct NA values
      if (cff(i,j) > NAN) cff(i,j) = 0.
   enddo
enddo

! Reverse cff
do j = 1, NZ
   do i = 1, NX
      cff1(i,j) = cff(i, NZ+1-j)
   enddo
enddo

! Interpolate cff1 ==> w
w0(:, 0)  = 0.
w0(:, NZ) = 0.
do i = 1, NX
   do j = 1, NZ-1
      w0(i,j) = (cff1(i,j)*RHz(j+1) + cff1(i,j+1)*RHz(j)) &
              / (RHz(j)+RHz(j+1))
   enddo
enddo

CALL NCREAD_3D(kappafile, 'kappat', NX, NZ, 1, cff)

do j = 1, NZ
   do i = 1, NX
      ! Correct NA values
      if (cff(i,j) > NAN) cff(i,j) = 0.
   enddo
enddo


! Reverse cff
do j = 1, NZ
   do i = 1, NX
      cff1(i,j) = cff(i, NZ+1-j)
   enddo
enddo

! Interpolate cff1 ==> kap
kap0(:, 0)  = 0.
kap0(:, NZ) = 0.
do i = 1, NX
   do j = 1, NZ-1
      kap0(i,j) = (cff1(i,j)*RHz(j+1) + cff1(i,j+1)*RHz(j)) &
                / (RHz(j)+RHz(j+1))
   enddo
enddo

return
END SUBROUTINE get_u
END MODULE
