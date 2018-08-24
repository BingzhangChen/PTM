MODULE GRID_2D
IMPLICIT NONE

! Number of grids along the X axis
integer, parameter :: NX      =    5098       

! Number of grids for u points at the interface
integer, parameter :: NX1     =    NX + 1       

! Number of grids along the Z axis
integer, parameter :: NZ      =    150       

! Number of grids for w points at the interface
integer, parameter :: NZ1     =    NZ + 1       

! Bottom depth
real     ::  Bot_dep(NX, 2)   =    -3000.

! Length along X for each grid
real     ::  dist(NX)         =    1D3     !Unit: m

! Horizontal index for horizontal grids (centered at middle points)
real     ::  X_r0(NX, 1)  =    0d0
real     ::  X_r(NX)      =    0d0
real     ::  X_w0(NX1,1)  =    0d0
real     ::  X_w(0:NX)    =    0d0

! Depth (m) of each grid
real     ::  Z_r(NZ)      =    0d0   
real     ::  Z_w(0:NZ)    =    0d0   


! Mask of each grid
integer  :: mask(NX, NZ)  =    0

! The index of the deepest grid with mask == 1
integer  ::   di(NX)      =    1

! Vertical distance of each grid
real     ::  Hz(NZ,1)     =    0d0   
real     ::  RHz(NZ)      =    0d0   

! U and w
real     :: u(0:NX,  NZ)      =    0d0
real     :: u1(0:NX, NZ)      =    0d0 ! For temporal interpolation
real     :: u2(0:NX, NZ)      =    0d0 ! For temporal interpolation

real     :: w(NX,  0:NZ)      =    0d0
real     :: w1(NX, 0:NZ)      =    0d0
real     :: w2(NX, 0:NZ)      =    0d0

real     :: kap(NX,  0:NZ)    =    0d0
real     :: kap1(NX, 0:NZ)    =    0d0
real     :: kap2(NX, 0:NZ)    =    0d0

CONTAINS

SUBROUTINE initialize_grid
IMPLICIT NONE
character(len=20) :: Xrfile    = 'X_r.dat'
character(len=20) :: Xwfile    = 'X_w.dat'
character(len=20) :: Depthfile = 'Depth.dat'
character(len=20) :: vertspace = 'vertspace.dat'
integer           :: i,k,j

! Read horizontal grid
CALL Readcsv(Xrfile, NX, 1, X_r0)
CALL Readcsv(Xwfile, NX1,1, X_w0)

! Write into X_r and X_w
do i = 1, NX
   X_r(i)   = X_r0(i,1)
   
   X_w(i-1) = X_w0(i,1)
enddo
X_w(NX) = X_w0(NX1,1)

! Calculate horizontal distance of each grid
do i = 1, NX
   dist(i)  = X_w(i) - X_w(i-1)
enddo

! Obtain bottom depth
CALL Readcsv(Depthfile, NX, 2, Bot_Dep)

! Obtain vertical spacing
CALL Readcsv(vertspace, NZ, 1, Hz)

! Reverse Hz ==> RHz
do k = 1, NZ
   RHz(k) = Hz(NZ-k+1, 1)
enddo

! Compute Z_r and Z_w
Z_w(0) =-3000.

DO k = 1, NZ
  Z_w(k) = Z_w(k-1) + RHz(k)
  Z_r(k) = Z_w(k-1) + RHz(k)/2.
ENDDO

! Determine the mask of each grid (deepest grid that is above the sea floor)
mask(:,:) = 1

DO j = 1, NX
    do i = NZ, 1, -1
        IF (Z_w(i-1) < (Bot_dep(j,2)+1D0) ) THEN
            di(j) = i+1
            do k = 1, i
               mask(j,k) = 0
            enddo
            EXIT
        ENDIF
    enddo
ENDDO

!!Check mask
!do j = 1, NX
!   write(6,*) 'Depth of the deepest valid grid :', Z_w(di(j)-1)
!   write(6,*) 'Bottom depth :', Bot_dep(j,2)
!   if (Z_w(di(j)-1) .lt. Bot_dep(j,2)) stop
!enddo
!stop
RETURN
END SUBROUTINE initialize_grid
END MODULE grid_2D

