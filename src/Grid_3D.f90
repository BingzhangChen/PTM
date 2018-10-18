MODULE GRID_3D
IMPLICIT NONE

! Number of Delaunay points 
integer, parameter :: numpoints = 21595

! Number of cells
integer, parameter :: NUMCELLS  = 41856

! Number of grids along the Z axis
integer, parameter :: NZ        = 100      

! Number of grids for w points at the interface
INTEGER, PARAMETER :: NZ1       = NZ + 1       

! Number of grids along the X axis
integer, parameter :: NX        = 300      

! Number of grids for u points at the interface
integer, parameter :: NX1       = NX + 1       

! Number of grids along the Y axis
integer, parameter :: NY        = 150      

! Number of grids for v points at the interface
integer, parameter :: NY1       = NY + 1       

! Bottom depth of original triangle domain
REAL     ::  Bot_dep(NUMCELLS,3)= -3000.

! Bottom depth of rectangle grids
REAL     ::      BOT(NX,NY)     = -3000.

! Depth (m) of each grid
REAL     ::  Z_r(NZ)            = 0d0   
REAL     ::  Z_w(0:NZ)          = 0d0   
REAL     ::  X_r(NX)            = 0d0   
REAL     ::  X_w(0:NX)          = 0d0   
REAL     ::  Y_r(NY)            = 0d0   
REAL     ::  Y_w(0:NY)          = 0d0   

! Mask of each grid (whole domain)
INTEGER  :: MASK(NX, NY, NZ)    = 0

! The index of the deepest grid with mask == 1
integer  ::    di(NX, NY)       = 1

! Vertical distance of each grid
real     ::  Hz(NZ)             = 0d0   

real     :: points(numpoints, 2)= 0d0
real     ::  cells(numcells,  5)= 0d0  !Cells of the whole domain

!SET COMPUTATION BOUNDARIES (IMPORTANT!!)
REAL, PARAMETER :: XMIN = -80D3, XMAX=0D0
REAL, PARAMETER :: YMIN = -160D3, YMAX=-120D3

REAL,    allocatable  :: nndist(:,:)
INTEGER, allocatable  ::  nnind(:,:)

! U, v, w and kappat
REAL     ::  u(0:NX, NY, NZ)  = 0d0
REAL     :: u1(0:NX, NY, NZ)  = 0d0 ! For temporal interpolation
REAL     :: u2(0:NX, NY, NZ)  = 0d0 ! For temporal interpolation

REAL     ::  v(NX, 0:NY, NZ)  = 0d0
REAL     :: v1(NX, 0:NY, NZ)  = 0d0 ! For temporal interpolation
REAL     :: v2(NX, 0:NY, NZ)  = 0d0 ! For temporal interpolation

REAL     ::  w(NX,NY, 0:NZ)   = 0d0
REAL     :: w1(NX,NY, 0:NZ)   = 0d0
REAL     :: w2(NX,NY, 0:NZ)   = 0d0

REAL     ::  kap(NX, NY,0:NZ) = 0d0
REAL     :: kap1(NX, NY,0:NZ) = 0d0
REAL     :: kap2(NX, NY,0:NZ) = 0d0
REAL     ::          DX       = 0D0
REAL     ::          DY       = 0D0


REAL, PARAMETER :: eps  = 1d-5   !A very small number

REAL, PARAMETER :: NAN  = 9D4    !Not a number in original output

CONTAINS

SUBROUTINE initialize_grid
IMPLICIT NONE
CHARACTER(LEN=8),  PARAMETER :: RUNDATA   = 'rundata/'
CHARACTER(LEN=8),  PARAMETER :: DATA      = 'data/'
CHARACTER(LEN=25), PARAMETER :: POINTFILE = trim(rundata)//'points.dat'
CHARACTER(LEN=25), PARAMETER :: EDGEFILE  = trim(rundata)//'edges.dat'
CHARACTER(LEN=25), PARAMETER :: CELLFILE  = trim(rundata)//'cells.dat'
CHARACTER(LEN=25), PARAMETER :: NCELLFILE = 'cellsNEW.dat'
CHARACTER(LEN=30), PARAMETER :: DEPTHFILE = trim(data)//'depth.dat-voro'
CHARACTER(LEN=25), PARAMETER :: VERTSPACE = trim(data)//'vertspace.dat'
CHARACTER(LEN=6 ), PARAMETER :: DISTFILE  = 'NNDIST'
CHARACTER(LEN=6 ), PARAMETER ::  INDFILE  = 'NNIND'
INTEGER,           PARAMETER :: DISTUNIT  = 31
INTEGER,           PARAMETER ::  INDUNIT  = 32
INTEGER,           PARAMETER :: CELLUNIT  = 33

INTEGER                      :: i,k,j,ig, L,cff2, cff3
REAL, ALLOCATABLE            :: cff(:,:), cff1(:)
REAL                         :: dist
LOGICAL                      :: THERE1, THERE2
REAL                         :: x1, y1, x2, y2, x3, y3

! Determine the Delaunay points and edges of each grid
! Read point file
ALLOCATE(cff(numpoints, 3))
CALL  READCSV(POINTFILE, NUMPOINTS, 3, .FALSE., CFF)
DO i=1,numpoints
   do j=1,2
      points(i,j)=cff(i,j)
   enddo
ENDDO
DEALLOCATE(cff)

! Read cells data
! Inquire whether the new cell file already exists
INQUIRE(FILE=NCELLFILE, EXIST=THERE1)

IF ((.NOT. THERE1) ) THEN
   WRITE(6,*) 'Cannot find the file ', trim(ncellfile),'. Create new file of ', trim(ncellfile)

   ALLOCATE(cff(numcells, 8))
   CALL  READCSV(CELLFILE, NUMCELLS, 8, .FALSE., CFF)

   ALLOCATE(nnind(numcells, 3))
   NNIND(:,:) = 0
   
   ALLOCATE(nndist(numcells, 3))
   NNDIST(:,:) = 0.

   DO i=1,NUMCELLS
      do j=1,5
         cells(i,j)=cff(i,j)
      enddo
   
      !Find the three points of each triangle
      do L=1, 3
         nndist(i, L) = 1d10
         SELECT CASE(L)
         CASE(1)
            cff2 = -999999
            cff3 = -999999
         CASE(2)
            cff2 = nnind(i, L-1)
            cff3 = -999999
         CASE(3)
            cff2 = nnind(i, L-1)
            cff3 = nnind(i, L-2)
         CASE DEFAULT
            stop "No such L!"
         END SELECT
   
         DO k=1, numpoints
            IF (k .ne. cff2 .and. k .ne. cff3) THEN
              DIST = SQRT( (points(k,1) - cells(i,1))**2 + (points(k,2) - cells(i,2))**2 )
              IF (DIST < nndist(i,L)) THEN  ! THe Lth nearest distance
                  nndist(i,L)   = DIST
                   nnind(i,L)   = k
                   cells(i,L+2) = real(k)
              ENDIF
            ENDIF
         ENDDO  ! END OF K
      ENDDO   ! END OF L
   ENDDO
   IF (ALLOCATED(cff))    DEALLOCATE(cff)
   IF (ALLOCATED(nnind))  DEALLOCATE(nnind)
   IF (ALLOCATED(nndist)) DEALLOCATE(nndist)

   ! Write out new cells file
   OPEN(CELLUNIT, FILE=NCELLFILE, STATUS='REPLACE', ACTION='WRITE')
 
   do i = 1, NUMCELLS
    ! Write into file
      WRITE(CELLUNIT, 1802) ( cells(i,L),  L = 1, 5)
   enddo
   
   CLOSE(CELLUNIT)
ELSE
   WRITE(6,*) 'READ FILE :', trim(NCELLFILE)
   CALL READCSV(NcellFILE, numcells, 5, .FALSE., cells)
ENDIF

! Obtain bottom depth
CALL READCSV(DEPTHFILE, NUMCELLS, 3, .FALSE., BOT_DEP)

! Obtain vertical spacing
IF (.NOT. (ALLOCATED(cff))) ALLOCATE(cff(NZ, 1))
CALL READCSV(VERTSPACE, NZ, 1, .FALSE., CFF)

! Reverse cff ==> Hz
DO k = 1, NZ
   Hz(k) = cff(NZ-k+1, 1)
ENDDO
DEALLOCATE(cff)

! Compute Z_r and Z_w
Z_w(0) = -MAXVAL(BOT_DEP(:,3))

DO k = 1, NZ
  Z_w(k) = Z_w(k-1) + Hz(k)
  Z_r(k) = Z_w(k-1) + Hz(k)/2.
ENDDO
Z_w(NZ)  = 0d0

! Compute horizontal distance of each rectangle grid
DX=(XMAX-XMIN)/DBLE(NX)
DY=(YMAX-YMIN)/DBLE(NY)

! Compute X_r and X_w
X_w(0) = XMIN

DO k = 1, NX
  X_W(k) = X_W(k-1) + DX
  X_R(k) = X_W(k-1) + DX/2.
ENDDO

Y_w(0) = YMIN

DO k = 1, NY
  Y_w(k) = Y_w(k-1) + DY
  Y_r(k) = Y_w(k-1) + DY/2.
ENDDO

! The first step of determining mask: to check if the grid has any point in it
mask(:,:,:) = 0
DO i = 1, NX
   DO j = 1, NY
      do k = 1, numcells
         ! Obtain the coordinates of each cell
         x1 = POINTS(INT(CELLS(k, 3)), 1)
         y1 = POINTS(INT(CELLS(k, 3)), 2)
         x2 = POINTS(INT(CELLS(k, 4)), 1)
         y2 = POINTS(INT(CELLS(k, 4)), 2)
         x3 = POINTS(INT(CELLS(k, 5)), 1)
         y3 = POINTS(INT(CELLS(k, 5)), 2)
         IF (POINTINTRIANGLE(x1, y1, x2, y2, x3, y3, X_r(i), Y_r(j))) THEN
             mask(i,j,:) = 1
             EXIT
         ENDIF
      enddo
   ENDDO
ENDDO

! Construct the database of Three nearest distances and cell 
IF (ALLOCATED(NNDIST)) DEALLOCATE(NNDIST)
IF (ALLOCATED(NNIND))  DEALLOCATE(NNIND)

ALLOCATE(nndist(NX*NY, 3))
nndist(:,:) = 0.

ALLOCATE(nnind(NX*NY, 3))
 nnind(:,:) = 0

! Inquire whether the files of nndist and nnind already exist
INQUIRE(FILE=DISTFILE, EXIST=THERE1)
INQUIRE(FILE= INDFILE, EXIST=THERE2)

IF ((.NOT. THERE1) .or. (.NOT. THERE2)) THEN
   WRITE(6,*) 'Cannot find the file ',distfile,'. Create new file of ', distfile, ' and ', indfile
   OPEN(DISTUNIT, FILE=DISTFILE, STATUS='REPLACE',ACTION='WRITE')
   OPEN(INDUNIT,  FILE=INDFILE, STATUS='REPLACE',ACTION='WRITE')
   do i=1,NX
      do j=1,NY
         ! select the three nearest points with mask == 1
         ig = (i-1)*NY + j
         do L=1, 3
            nndist(ig, L) = 1d10
            select case(L)
            case(1)
               cff2 = -999999
               cff3 = -999999
            case(2)
               cff2 = nnind(ig, L-1)
               cff3 = -999999
            case(3)
               cff2 = nnind(ig, L-1)
               cff3 = nnind(ig, L-2)
            case default
               stop "No such L!"
            end select

            do k=1,numcells
               IF (k .ne. cff2 .and. k .ne. cff3) THEN
                 DIST = SQRT( (X_r(i) - cells(k,1))**2 + (Y_r(j) - cells(k,2))**2 )
                 IF (DIST < nndist(ig,L)) THEN  ! THe Lth nearest distance
                     nndist(ig,L) = DIST
                      nnind(ig,L) = k
                 ENDIF
               ENDIF
            enddo  ! end of k
         enddo   ! end of L

         ! Write into file
         WRITE(DISTUNIT, 1800) (    nndist(ig,L),  L = 1, 3)
         WRITE(INDUNIT,  1801) (REAL(NNIND(ig,L)), L = 1, 3)
      enddo
   enddo

   CLOSE(DISTUNIT)
   CLOSE(INDUNIT)

ELSE
   WRITE(6,*) 'Read files ',distfile,' and ', indfile
   CALL READCSV(DISTFILE, NX*NY, 3, .FALSE., NNDIST)

   IF (ALLOCATED(CFF)) DEALLOCATE(CFF)

   ALLOCATE(cff(NX*NY,3))
   cff(:,:)=0.
   CALL READCSV( INDFILE, NX*NY, 3, .FALSE., cff)

   do k=1, NX*NY
      do i=1,3
         nnind(k,i) = int(cff(k,i))
      enddo
   enddo
   DEALLOCATE(cff)
ENDIF

! Interpolate bottom depth
IF (.NOT. ALLOCATED(cff1)) ALLOCATE(cff1(numcells))
cff1(:)=Bot_Dep(:,3)
call KNN3(cff1, BOT)
BOT = -BOT

! Determine the vertical mask of each grid (deepest grid that is above the sea floor)

DO i = 1, NX
  DO j = 1, NY
     IF (mask(i,j,NZ) == 1) THEN
        DO k = NZ, 1, -1
            IF (Z_w(k-1) < BOT(i,j) ) THEN
                di(i,j) = k+1
                do ig = 1, k
                   mask(i,j,ig) = 0
                enddo
                EXIT
            ENDIF
        ENDDO

        IF (BOT(i,j) > -10.) THEN  ! Remove very shallow grids
           mask(i,j,:) = 0
        ENDIF
     ENDIF
  ENDDO
ENDDO
IF (ALLOCATED(cff1))   DEALLOCATE(cff1)

RETURN
1800 FORMAT(3(F12.3, 2X))
1801 FORMAT(3(F12.0, 2X))
1802 FORMAT(2(F12.2, 2X), 3(F12.0, 2X))
END SUBROUTINE initialize_grid

!Copied from http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
LOGICAL FUNCTION POINTINTRIANGLE(x1, y1, x2, y2, x3, y3, x, y)
IMPLICIT NONE
real, intent(in) :: x1, y1, x2, y2, x3, y3, x, y
real             :: denominator, a, b, c
real, parameter  :: z1 = -eps
real, parameter  :: z2 = 1.+eps

denominator = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))
a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denominator
b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denominator
c = 1. - a - b
if (a .ge. z1 .and. a .le. z2 .and. b .ge. z1 .and. b .le. z2 .and. c .ge. z1 .and. c .le. z2) then
  pointInTriangle = .true.
else
  pointInTriangle = .false.
endif
END FUNCTION POINTINTRIANGLE

SUBROUTINE KNN3(X, XOUT)
IMPLICIT NONE

REAL, INTENT(IN)     ::     X(NUMCELLS)  ! INPUT VECTOR CONTAINING VALUES FOR EACH CELL
REAL, INTENT(OUT)    ::  XOUT(NX, NY)

INTEGER              :: i, j, k, ig

! Whether already found a very close point?
LOGICAL              :: interp    = .FALSE.
REAL                 :: suminvw   = 0.
REAL                 :: suminvwd  = 0.

IF (.NOT. ALLOCATED(NNIND))  STOP "NNIND MISSING!"
IF (.NOT. ALLOCATED(NNDIST)) STOP "NNDIST MISSING!"

! Interpolate bottom depth
DO i=1,NX
   do j=1,NY
      ! select the three nearest points
      ig=(i-1)*NY+j
      interp=.FALSE.
      do k = 1,3
         if (NNdist(ig, k) < eps) then
             XOUT(i,j) = X(nnind(ig,k))
             interp    = .TRUE.
             EXIT
         end if
      enddo
     
      IF (.NOT. interp) THEN
         suminvw = 0.
         suminvwd= 0.
         do k = 1,3
            if (X(nnind(ig,k)) < NAN) then ! If the original value is valid
               suminvw = suminvw + (1./NNdist(ig,k))**2
               suminvwd= suminvwd+ (1./NNdist(ig,k))**2 * X(nnind(ig,k)) 
            endif
         enddo

         if (suminvw > 0.) then
            XOUT(i,j)= suminvwd/suminvw
         else
            XOUT(i,j)= 0.   ! A ghost point
         endif
      ENDIF
   enddo
ENDDO

return
END SUBROUTINE KNN3
END MODULE GRID_3D

