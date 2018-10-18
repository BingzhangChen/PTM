! Lagrange particle tracking on one particle
! !INTERFACE:
SUBROUTINE LAGRANGE(nlev,zlev,mask,nuh,w,zi,zp)
USE PARAMS_MOD, only: dt
USE Grid_3D,    only: eps
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
LOGICAL,PARAMETER  :: visc_corr  =.TRUE.
LOGICAL,PARAMETER  :: bounce_back=.FALSE.
!EOP
!-----------------------------------------------------------------------
!BOC

!Judge whether the input is correct
IF (zp < zlev(0) .or. zp > zlev(nlev) .or. zi < 1 .or. zi > nlev) THEN
   stop "The particle out of the domain!"
ENDIF

dt_inv      = 1./dt
rnd_var_inv = 1./rnd_var

CALL RANDOM_NUMBER(rnd)

rnd = 2.*rnd-1.

DO i = 1,nlev
   
   ! Depth of each grid
   dz(i) = zlev(i)-zlev(i-1)

   ! gradient of Kv
   dzn(i)= (nuh(i)-nuh(i-1))/dz(i)
END DO


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
      if (bounce_back) then
         zloc = 2.*zlev(nlev) - zloc 
      else
         zloc = zlev(nlev)-eps
      endif
      zi   = nlev
    ! Judge whether after bouncing the particle overshoots through a grid. 
      if (zloc < zlev(nlev-1)) then 
          stop "The particle overshoots! Reduce time step!"
      endif
   elseif (zloc < zlev(0)) then     
      if (bounce_back) then
         zloc = 2d0*zlev(0) - zloc      
      else
         zloc = zlev(0) +eps
      endif
      zi   = 1
      if (zloc > zlev(1)) then 
          stop "The particle overshoots! Reduce time step!"
      endif
   endif      

   ! To judge whether the particle enters a fake grid (mask = 0)
   ! If yes, then bounce back
   IF (step > 0.) THEN
      if (zloc > zlev(zi)) then
         if (mask(zi + 1) == 0) then
             if (bounce_back) then
                zloc = 2d0*zlev(zi) - zloc
             else
                zloc = zlev(zi) - eps
             endif
        ! Judge whether after bouncing the particle overshoots through a grid. 
             if (zloc < zlev(zi-1)) then 
                 STOP "The particle overshoots! Reduce time step!"
             endif
         endif
      endif
   ELSEIF (step < 0.) THEN
      if (zloc < zlev(zi-1)) then
         if (mask(zi - 1) == 0) then
             if (bounce_back) then
                zloc = 2d0*zlev(zi-1) - zloc
             else
                zloc = zlev(zi-1) + eps
             endif

        ! Judge whether after bouncing the particle overshoots through a grid. 
             if (zloc > zlev(zi)) then 
                 STOP "The particle overshoots! Reduce time step!"
             endif
         endif
      endif
   ENDIF

   step = zloc-zp     ! The distance traveled

! The following is to find the index for the grid where the particle currently resides
   IF (STEP > 0.) THEN ! SEARCH NEW INDEX ABOVE OLD INDEX
      do i=zi, nlev
         if (zlev(i) .gt. zloc) then
             zi = i
             EXIT
         endif
      end do
   ELSE                ! SEARCH NEW INDEX BELOW OLD INDEX
      do i=zi,1,-1
         if (zlev(i-1) .lt. zloc) then
             zi=i
             EXIT
         endif
      end do
   END IF
ELSE
   i   =zi
   zloc=zp
END IF

rat  = (zloc-zlev(i-1))/dz(i)
visc = rat*nuh(i)+(1.-rat)*nuh(i-1)  ! interpolate Kv

IF (visc.lt.visc_back) visc=visc_back   !Kv background

zp_old = zp

step = dt*( SQRT(2d0*rnd_var_inv*dt_inv*visc) * rnd + w(i) + dzn(i) )

zp = zp + step

! Judge whether the particle overshoots through a grid. 
! If yes, the time step needs to be shortened
IF (i .lt. nlev .and. i .gt. 1) THEN
  if (zp > zlev(i+1) .or. zp < zlev(i-2)) then
    stop "The particle overshoots! Reduce time step!"
  endif
ENDIF

! First judge whether the particle jumps out of the domain
IF (zp > zlev(nlev)) THEN
   if (bounce_back) then
      zp = 2.*zlev(nlev) - zp 
   else
      zp = zlev(nlev) - eps
   endif

 ! Judge whether after bouncing the particle overshoots through a grid. 
   if (zp < zlev(nlev-1)) then 
       stop "The particle overshoots! Reduce time step!"
   endif
ELSEIF (zp < zlev(0)) THEN     
   if (bounce_back) then
      zp = 2d0*zlev(0) - zp      
   else
      zp = zlev(0) + eps
   endif
   if (zp > zlev(1)) then 
       stop "The particle overshoots! Reduce time step!"
   endif
ENDIF      

! To judge whether the particle enters a fake grid (mask = 0)
! If yes, then bounce back
if (step > 0.) then
   if (zp > zlev(i)) then
      if (mask(i + 1) == 0) then
          if (bounce_back) then
            zp = 2.*zlev(i) - zp
          else
            zp = zlev(i) - eps
          endif
     ! Judge whether after bouncing the particle overshoots through a grid. 
          if (zp < zlev(i-1)) then 
              stop "The particle overshoots! Reduce time step!"
          endif
      endif
   endif
elseif (step < 0.) then
   if (zp < zlev(i-1)) then
      if (mask(i - 1) == 0) then
          if (bounce_back) then
             zp = 2.*zlev(i-1) - zp
          else
             zp = zlev(i-1) + eps
          endif

     ! Judge whether after bouncing the particle overshoots through a grid. 
          if (zp > zlev(i)) then 
              stop "The particle overshoots! Reduce time step!"
          endif
      endif
   endif
endif

!Compute new zi
do i=1,nlev
   if (zlev(i) .gt. zp .and. zlev(i-1) .le. zp) then
       zi = i
       EXIT
   endif
end do

return
END SUBROUTINE LAGRANGE


