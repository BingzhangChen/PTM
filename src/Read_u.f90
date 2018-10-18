MODULE Read_u
USE NETCDF_IO
USE GRID_3D, only: numcells, NZ, mask, BOT, NX, NY
USE GRID_3D, only: NX1, NY1, NZ1, KNN3, DX, DY, HZ
USE GRID_3D, only: u,v,w,kap,u1,v1,w1,kap1,u2,v2,w2,kap2
USE GRID_3D, only: NNIND
IMPLICIT NONE

CONTAINS

SUBROUTINE interp_u(t, step) 
IMPLICIT NONE

! Current model time
integer, intent(in)   :: t


! Current step in offline SUNTANS data
integer, intent(in)   :: step

integer               :: i,j,k
integer               :: t1, t2

! Calculate two time points of step and step + 1
t1 = (step-1)*period
t2 =  step   *period

IF (t < t1 .or. t > t2) THEN
    STOP "Current time NOT within the time interval! Something wrong!"
ELSE
    do i = 1, NX-1
       do j = 1, NY
          do k = 1, NZ
             u(i,j,k) = (u1(i,j,k)*(t2 - t) + u2(i,j,k)*(t - t1))/(t2 - t1)
          enddo
       enddo
    enddo

    do i = 1, NX
       do j = 1, NY-1
          do k = 1, NZ
             v(i,j,k) = (v1(i,j,k)*(t2 - t) + v2(i,j,k)*(t - t1))/(t2 - t1)
          enddo
       enddo
    enddo

    do i = 1, NX
       do j = 1, NY
          do k = 1, NZ-1
            w(i,j,k)=  (w1(i,j,k)*(t2 - t) +   w2(i,j,k)*(t - t1))/(t2 - t1)
          kap(i,j,k)=(kap1(i,j,k)*(t2 - t) + kap2(i,j,k)*(t - t1))/(t2 - t1)
          enddo
       enddo
    enddo
END IF
RETURN
END SUBROUTINE interp_u

SUBROUTINE get_u(rec, u0, v0, w0, kap0) 
IMPLICIT NONE

! The record to be read from nc files (u.nc)
INTEGER, intent(in)   :: rec
REAL,    intent(out)  :: u0(0:NX,NY,NZ), v0(NX,0:NY,NZ), w0(NX,NY,0:NZ), kap0(NX,NY,0:NZ)
REAL,    ALLOCATABLE  :: cffw(:,:,:)
INTEGER               :: i,j,k,ig
INTEGER               :: u_varid, v_varid, w_varid, kap_varid
INTEGER               :: NCID

! Scratch matrix to store temporary data read from nc file
REAL                  :: cff(numcells,NZ) = 0d0
REAL                  :: cff1(numcells)   = 0d0
REAL                  :: cff2(NX, NY)     = 0d0

! Filename for u and w
CHARACTER(LEN=10)     :: ufile     = 'u.nc'
CHARACTER(LEN=10)     :: kappafile = 'kappat.nc'

! Inquire whether Flow.nc exists
INQUIRE(FILE=FLOW_FNAME, EXIST=FLOWFILEEXISTS)

IF (FLOWFILEEXISTS) THEN
    IF (REC .EQ. 1) WRITE(6,'(A30, A8)') 'Read current velocities from', trim(Flow_Fname)
    CALL check(NF90_OPEN(FLOW_FNAME, NF90_NOWRITE, NCID) )

    ALLOCATE(cffw(NX,NY,NZ))

    ! Get the varid of the data variable, based on its name.
    CALL CHECK(NF90_INQ_VARID(NCID, u_name, u_varid) )

    ! Read the data.
    ! step is the NO. of the record
    CALL CHECK(NF90_GET_VAR(ncid, u_varid, cffw, start = (/1, 1, 1, rec/), &
                                                 count = (/NX,NY,NZ,1/)  ))
    u0(0, :,:) = 0.
    u0(NX,:,:) = 0.

    do i = 1, NX-1
       u0(i,:,:)=(cffw(i,:,:) + cffw(i+1,:,:))/2.
    enddo

    ! reset cffw
    cffw(:,:,:)=0d0

    CALL CHECK(NF90_INQ_VARID(NCID, v_name, v_varid) )
    CALL CHECK(NF90_GET_VAR(NCID, v_varid, cffw, start = (/1, 1, 1, rec/), &
                                                 count = (/NX,NY,NZ,1/)  ))

    v0(:,0, :) = 0.
    v0(:,NY,:) = 0.
    DO i = 1, NY-1
       v0(:,i,:)=(cffw(:,i,:) + cffw(:,i+1,:))/2.
    ENDDO

    DEALLOCATE(cffw)

    ALLOCATE(cffw(NX, NY, NZ1))
    cffw(:,:,:)=0d0

    ! Read w
    CALL CHECK(NF90_INQ_VARID(NCID, W_NAME, W_VARID) )
    CALL CHECK(NF90_GET_VAR(NCID, W_VARID, cffw, start = (/1, 1, 1, rec/), &
                                                 count = (/NX,NY,NZ1,1/)  ))

    do i = 0, NZ
       w0(:,:,i)=cffw(:,:,i+1)
    enddo

    cffw(:,:,:)=0d0

    ! Read kap
    CALL CHECK(NF90_INQ_VARID(NCID, KAP_NAME, KAP_VARID) )
    CALL CHECK(NF90_GET_VAR(NCID, KAP_VARID, cffw,start = (/1, 1, 1, rec/), &
                                                  count = (/NX,NY,NZ1,1/)  ))

    do i = 0, NZ
       kap0(:,:,i)=cffw(:,:,i+1)
    enddo

    DEALLOCATE(cffw)
    ! Close the file, freeing all resources.
    CALL CHECK(NF90_CLOSE(NCID) )
ELSE
    ! Read the external u,  w, kappa to set up initial condition
    CALL NCREAD_3D(ufile, 'u', numcells, NZ, rec, cff)
    
    u0(0, :,:) = 0.
    u0(NX,:,:) = 0.
    
    DO j = 1,NZ
       do ig = 1, numcells
          cff1(ig)=cff(ig,j)
       enddo
    
       ! Interpolate horizontally
       call KNN3(cff1, cff2)
    
       ! Interpolate cff2 ==> u0
       do k = 1, NY
          do i = 1, NX-1
             u0(i,k,NZ-j+1) = (cff2(i,k) + cff2(i+1,k))/2.
          enddo
       enddo
    
    ENDDO
    
    CALL NCREAD_3D(ufile, 'v', numcells, NZ, rec, cff)
    
    v0(:,0, :) = 0.
    v0(:,NY,:) = 0.
    
    DO j = 1, NZ
       do ig = 1, numcells
          cff1(ig)=cff(ig,j)
       enddo
    
       ! Interpolate on the surface
       call KNN3(cff1, cff2)
    
       ! Interpolate cff2 ==> v0
       do k = 1, NX
          do i = 1, NY-1
             v0(k,i,NZ-j+1) = (cff2(k,i) + cff2(k,i+1))/2.
          enddo
       enddo
    ENDDO
    
    CALL NCREAD_3D(ufile, 'w', numcells, NZ, rec, cff)
    
    w0(:,:, 0) = 0.
    w0(:,:,NZ) = 0.
    
    ALLOCATE(cffw(NX,NY,NZ))
    cffw(:,:,:)=0.
    DO j = 1, NZ
    
       do ig = 1, numcells
          cff1(ig)=cff(ig,j)
       enddo
    
       ! Interpolate on the surface
       call KNN3(cff1, cff2)
    
       ! Interpolate cff2 ==> cffw
       do k = 1, NX
          do i = 1, NY
             cffw(k,i,j) = cff2(k,i)
          enddo
       enddo
    ENDDO
    
    ! Interpolate cffw ==> w0
    DO i = 1, NX
       do j = 1, NY
          do k = 1,NZ-1
             w0(i,j,NZ-k) = (cffw(i,j,k)*Hz(NZ-k) + cffw(i,j,k+1)*Hz(NZ-k+1)) &
                          / (Hz(NZ-k)+Hz(NZ-k+1))
          enddo
       enddo
    ENDDO
    
    CALL NCREAD_3D(kappafile, 'kappat', NX, NZ, rec, cff)
    
    kap0(:,:, 0) = 0.
    kap0(:,:,NZ) = 0.
    
    cffw(:,:,:)=0.
    
    DO j = 1, NZ
       do ig = 1, numcells
          cff1(ig)=cff(ig,j)
       enddo
    
       ! Interpolate on the surface
       call KNN3(cff1, cff2)
    
       ! Interpolate cff2 ==> cffw
       do k = 1, NX
          do i = 1, NY
             cffw(k,i,j) = cff2(k,i)
          enddo
       enddo
    ENDDO
    
    ! Interpolate cffw ==> kap0
    DO i = 1, NX
       do j = 1, NY
          do k = 1,NZ-1
             kap0(i,j,NZ-k) = (cffw(i,j,k)*Hz(NZ-k)   &
                          + cffw(i,j,k+1)*Hz(NZ-k+1)) &
                          / (Hz(NZ-k)+Hz(NZ-k+1))
          enddo
       enddo
    ENDDO
    DEALLOCATE(cffw)
ENDIF
return
END SUBROUTINE get_u
END MODULE Read_u
