MODULE PARAMS_MOD

! DEFINE VARIOUS MODEL SETTINGS
IMPLICIT NONE

! Number of time steps
INTEGER, PARAMETER :: NDays     =    10

! Number of particle groups
INTEGER, PARAMETER :: NG        =    4

! Number of particles of each group (to be assigned, and sum should == TNPAR)
INTEGER, PARAMETER :: N_PAR     =    500
 
! Total number of particles
INTEGER, PARAMETER :: TNPAR     =    NG*N_PAR

! swimming speed (m/s) of particles of each group
REAL               :: wb(NG)    =    1.

! Timestep (unit: seconds)
REAL,    PARAMETER :: dt        =    30.

INTEGER, PARAMETER :: uniform   = 1   
INTEGER, PARAMETER :: point     = 2   
INTEGER            :: initial_condition = point
LOGICAL, PARAMETER :: RESTART   = .FALSE.
character(LEN=10), parameter :: restartname = 'restart.nc'
END Module PARAMS_MOD
