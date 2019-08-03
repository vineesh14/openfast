! ==============================================================================
SUBROUTINE FVW_READ_WAKE_PARAM( pFVW )!, uFVW )

  USE FVW_Parm
  USE AeroDyn14_Types, Only: FVW_ParameterType, AD14_InputType
  USE MultTurb_Params, Only: NumTurbs, TurbDist, PerUinf, TurbLocs !mype, npes
!!KS -- had to add TurbLocs here 7.1.19
  USE FileManipulation, Only: WriteInitWake
  USE NWTC_Library

  IMPLICIT NONE
  !INCLUDE 'mpif.h'

  TYPE( FVW_ParameterType ), INTENT( INOUT ) :: pFVW ! Parameters
  INTEGER(IntKi)  :: ErrStat           ! Error status of the operation
  REAL( ReKi ) :: LMax, LTot, TotRad, Uinf, CUTOFF_dist, TurbDist_Single

  INTEGER :: CUTOFF_diff, I, ierr
  LOGICAL :: file_exists, ForTwoTurbs
  CHARACTER( * ), PARAMETER :: root = 'InputFiles/'

  REAL(ReKi):: yloc, zloc

  NumTurbs = 1 !!!KS -- HACKHACKHACK
  ALLOCATE( TurbLocs( NumTurbs, 2 ))

  OPEN(unit=100,file="MultTurbineInputs.txt")
  !Reading in number of turbines and how they are aligned
  READ(100,*) NumTurbs
  PRINT*, 'Number of Turbines:  ', NumTurbs

  READ( 100, * ) yloc, zloc

  TurbLocs( 1, 1 ) = yloc
  TurbLocs( 1, 2 ) = zloc

  !Set local parameters for ease of use.
  DtAero          = pFVW%DtAero
  TMax            = pFVW%TMax
  NElm            = pFVW%NElm
  IF ( .NOT. ALLOCATED( Chord )) ALLOCATE( Chord( NElm ), RElm( NElm ))
  RElm            = pFVW%RElm
  NumBl           = pFVW%NumBl
  Radius          = pFVW%Radius
  Chord           = pFVW%C
  HubHt              = pFVW%HubHt
  HubR            = pFVW%HubRad
  RotSpeed        = pFVW%RotSpeed  !!KS -- added 6.28.19

  FVW_AirfoilParm = pFVW%AirfoilParm
  FVW_AirfoilOut  = pFVW%AirfoilOut

  open( unit=23, file="Input_FVW_parameters.txt" )

  read( 23, * ) CUTOFF_dist
  read( 23, * ) Uinf
  read( 23, * ) PerUinf
  read( 23, * ) NumBS
  read( 23, * ) Root_cut
  read( 23, * ) eps
  read( 23, * ) nu
  read( 23, * ) near_deg
  read( 23, * ) Nelm_start
  read( 23, * ) RotSpeed_Est
  !read( 23, * ) ForTwoTurbs
  !IF (ForTwoTurbs .EQV. .TRUE. .AND. NumTurbs .LT. 2) THEN
  !   read( 23, * ) TurbDist_Single
  !END IF
  close( 23 )

  CUTOFF_dist = CUTOFF_dist*(Radius+HubR)*2.0_ReKi
  near_deg = near_deg / 180.00_ReKi * Pi_D
  !IF ( mype .EQ. 0 ) THEN
  PRINT *, 'RotSpeed is ', RotSpeed
  PRINT *, 'RotSpeed_Est is ', RotSpeed_Est
  PRINT *, 'dtaero is ', DtAero
  PRINT *, 'CUTOFF_dist is ', CUTOFF_dist
  !END IF
  delta_psi = DtAero * RotSpeed
  delta_psi_Est = DtAero * RotSpeed_Est

  CUTOFF_prim = CUTOFF_dist / ( PerUinf*Uinf * DtAero ) ! # of markers
       !NINT(TwoPi_D/delta_psi) is the # of markers per revolution
  PRINT*, 'delta_psi: ', delta_psi(1), '; delta_psi_Est: ', delta_psi_Est
  PRINT*, 'CUTOFF_prim: ', CUTOFF_prim
  Nj = Ceiling( TMax/DtAero )   !# of time steps
  Nj2 = 20 * NINT( TwoPi_D / delta_psi(1) ) + Nj
  NnearMax = NINT( near_deg / delta_psi(1) ) + 1
  ALLOCATE( CUTOFF_upinit( NumTurbs ), CUTOFF_upmax( NumTurbs ), CUTOFF_up( NumTurbs ), BC( NumTurbs ))
  BC = 10
  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  !IF ( mype .EQ. 0 .AND. NumTurbs .GT. 1 ) THEN

  !   DO I = 1, NumTurbs
  !         PRINT*, 'TurbDist: ', TurbDist(I), 'mype: ', mype
  !         CUTOFF_diff = TurbDist(I) / ( PerUinf*Uinf * DtAero )
  !         CUTOFF_upinit( I ) = CUTOFF_prim + CUTOFF_diff
  !         CUTOFF_upmax( I ) = CUTOFF_upinit( I )
  !         CUTOFF_up( I ) = CUTOFF_upinit( I )
  !   END DO
  !ELSE IF (mype .EQ. 0 .AND. ForTwoTurbs .EQV. .TRUE. ) THEN
  !   CUTOFF_diff = TurbDist_Single / ( PerUinf*Uinf * DtAero )
  !   CUTOFF_upinit(1) = CUTOFF_prim+CUTOFF_diff
  !   CUTOFF_upmax(1) = CUTOFF_upinit(1)
  !   CUTOFF_up(1) = CUTOFF_upinit(1)
  !ELSE
  CUTOFF_upinit(1) = CUTOFF_prim
  CUTOFF_upmax(1) = CUTOFF_upinit(1)
  CUTOFF_up(1) = CUTOFF_upinit(1)
  !END IF

  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  !CALL MPI_BCAST( CUTOFF_upinit, NumTurbs, MPI_INTEGER, 0, MPI_COMM_WORLD, ErrStat )
  !CALL MPI_BCAST( CUTOFF_upmax, NumTurbs, MPI_INTEGER, 0, MPI_COMM_WORLD, ErrStat )
  !CALL MPI_BCAST( CUTOFF_up, NumTurbs, MPI_INTEGER, 0, MPI_COMM_WORLD, ErrStat )

  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  !IF ( mype .EQ. 0 ) THEN
  PRINT*, 'NnearMax is ', NnearMax
  PRINT*, 'Nj is ', Nj
  PRINT*, 'Nj2 is ', Nj2
  PRINT*, 'CUTOFF_upmax is: ', CUTOFF_upmax      !Not the same for each turbine/processor
  PRINT*, 'CUTOFF_upinit is: ', CUTOFF_upinit
  !END IF
  !Check if InitialWake.txt exists; if it does, now inital wakes for each turbine will be written;
  !if it does not, current files will be used.
  !IF ( mype .EQ. 0 ) THEN
  INQUIRE(FILE="InitialWake.txt", EXIST=file_exists)
  IF (file_exists .EQV. .TRUE.) THEN
     CALL WriteInitWake( CUTOFF_upinit )
  END IF
  !END IF

END SUBROUTINE FVW_READ_WAKE_PARAM
! ==============================================================================

! ==============================================================================
!FIXME: is this an initialize once routine, or is this something that might change periodically????
SUBROUTINE FVW_INITIALIZE_WAKE(  )

  USE FVW_Parm
  USE FVW_vars
!FIXME: why do we need AD14 types here????????
  USE AeroDyn14_Types, Only: FVW_ParameterType, AD14_InputType
  USE FVW_ComputeWake, Only: BladeQuarterChordjm1, BladeQuarterChordjm2
  USE MultTurb_Params
!FIXME: remove inflowwind from module entirely
  USE InflowWind

  IMPLICIT NONE
  !INCLUDE 'mpif.h'

  INTEGER nindx, kindx, jindx, nindx2, kindx2, kindx3, nbs, j2, ierr
  CHARACTER( * ), PARAMETER :: root = 'InitFiles/'

  NumWakes = NumTurbs * NumBl
  CUTOFF_Allocate = MAXVAL(CUTOFF_upmax) + NINT( TwoPi_D / delta_psi_Est )


  WakeAgeLimit = NINT( TwoPi_D / delta_psi_Est ) + MAXVAL(CUTOFF_upmax)

PRINT*,' WakeAgeLimit: ', WakeAgeLimit

  !CALL MPI_BCAST( CUTOFF_Allocate, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  PRINT*, ' CUTOFF_Allocate is: ', CUTOFF_Allocate
  IF ( .NOT. ALLOCATED( FWake%rj ))  THEN
     ALLOCATE ( CUTOFF( NumTurbs ), FWake%rj( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%rjm1( 3, CUTOFF_Allocate, NumWakes ), FWake%rjm2( 3, CUTOFF_Allocate, NumWakes ), &
        & NWake%r_nearj( 3, NumBS+1, NnearMax, NumBl ), NWake%r_nearjm1( 3, NumBS+1, NnearMax, NumBl ), &
        & NWake%r_nearjm2( 3, NumBS+1, NnearMax, NumBl ), FWake%r_primej( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_primejm1( 3, CUTOFF_Allocate, NumWakes ), FWake%r_primejm2( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_primejm3( 3, CUTOFF_Allocate, NumWakes ), FWake%r_oldj( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_oldjm1( 3, CUTOFF_Allocate, NumWakes ), FWake%r_newj( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_oldjm2( 3, CUTOFF_Allocate, NumWakes ), FWake%r_oldjm3( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_newjm1( 3, CUTOFF_Allocate, NumWakes ), FWake%r_newjm2( 3, CUTOFF_Allocate, NumWakes ), &
        & FWake%r_newjm3( 3, CUTOFF_Allocate, NumWakes ), &
        & NWake%Gammablj( NumBS, NumBl ), FWake%Gammaj( CUTOFF_Allocate, NumWakes ), &
        & NWake%Gamma_nearj( NnearMax, NumBS+1, NumBl ), FWake%Gammajp1( CUTOFF_Allocate, NumWakes ), &
        & NWake%Gamma_nearjp1( NnearMax, NumBS+1, NumBl ), NWake%Gammabljm1( NumBS, NumBl ), &
        & FWake%Gammajm1( CUTOFF_Allocate, NumWakes ), NWake%Gamma_nearjm1( NnearMax, NumBS+1, NumBl ), &
        & FWake%VinducedFarWakej( 3, CUTOFF_Allocate, NumBl, WakeAgeLimit, NumWakes ), &
        & FWake%VinducedFarWakejm1( 3, CUTOFF_Allocate, NumBl, WakeAgeLimit, NumWakes ), &
        & FWake%VinducedFarWakeRj( 3, NumBS, NumWakes, WakeAgeLimit, NumWakes ), &
        & FWake%VinducedFarWakeRjm1( 3, NumBS, NumWakes, WakeAgeLimit, NumWakes ), &
        & BladeQuarterChordjm1( 3, NumBS+1, NumBl ), BladeQuarterChordjm2( 3, NumBS+1, NumBl ), &
        & loop( NumTurbs ))
  END IF

  FWake%rj = 0.00_ReKi; FWake%rjm1 = 0.00_ReKi; FWake%rjm2 = 0.00_ReKi; NWake%r_nearj = 0.00_ReKi; NWake%r_nearjm1 = 0.00_ReKi
  FWake%r_primej = 0.00_ReKi; FWake%r_primejm1 = 0.00_ReKi; FWake%r_primejm2 = 0.00_ReKi; FWake%r_primejm3 = 0.00_ReKi
  FWake%r_oldj = 0.00_ReKi; FWake%r_oldjm1 = 0.00_ReKi; FWake%r_oldjm2 = 0.00_ReKi; FWake%r_oldjm3 = 0.00_ReKi
  FWake%r_newj = 0.00_ReKi; FWake%r_newjm1 = 0.00_ReKi; FWake%r_newjm2 = 0.00_ReKi; FWake%r_newjm3 = 0.00_ReKi
  NWake%Gammablj = 0.00_ReKi; FWake%Gammaj = 0.00_ReKi; NWake%Gamma_nearj = 0.00_ReKi
  NWake%Gamma_nearjp1 = 0.00_ReKi; FWake%Gammajp1 = 0.00_ReKi; NWake%Gammabljm1 = 0.00_ReKi; FWake%Gammajm1 = 0.00_ReKi
  NWake%Gamma_nearjm1 = 0.00_ReKi; BladeQuarterChordjm1 = 0.00_ReKi; BladeQuarterChordjm2 = 0.00_ReKi
  FWake%VinducedFarWakeRj = 0.00_ReKi; FWake%VinducedFarWakeRjm1 = 0.00_ReKi
  FWake%VinducedFarWakej = 0.00_ReKi; FWake%VinducedFarWakejm1 = 0.00_ReKi

  OPEN( unit = 1000, file = ( root//TRIM( 'Turb1' )//'_r_primej.txt'   ))
  OPEN( unit = 1001, file = ( root//TRIM( 'Turb1' )//'_r_primejm1.txt' ))
  OPEN( unit = 1002, file = ( root//TRIM( 'Turb1' )//'_r_primejm2.txt' ))
  OPEN( unit = 1003, file = ( root//TRIM( 'Turb1' )//'_r_primejm3.txt' ))

  OPEN( unit = 1025, file = ( root//TRIM( 'Turb1' )//'_r_oldj.txt'   ))
  OPEN( unit = 1004, file = ( root//TRIM( 'Turb1' )//'_r_oldjm1.txt' ))
  OPEN( unit = 1005, file = ( root//TRIM( 'Turb1' )//'_r_oldjm2.txt' ))
  OPEN( unit = 1006, file = ( root//TRIM( 'Turb1' )//'_r_oldjm3.txt' ))

  OPEN( unit = 1007, file = ( root//TRIM( 'Turb1' )//'_r_newjm1.txt' ))
  OPEN( unit = 1008, file = ( root//TRIM( 'Turb1' )//'_r_newjm2.txt' ))
  OPEN( unit = 1009, file = ( root//TRIM( 'Turb1' )//'_r_newjm3.txt' ))

  OPEN( unit = 1010, file = ( root//TRIM( 'Turb1' )//'_rjm1.txt' ))
  OPEN( unit = 1011, file = ( root//TRIM( 'Turb1' )//'_rjm2.txt' ))

  OPEN( unit = 1012, file = ( root//TRIM( 'Turb1' )//'_r_nearjm1.txt' ))
  OPEN( unit = 1013, file = ( root//TRIM( 'Turb1' )//'_r_nearjm2.txt' ))

  OPEN( unit = 1017, file = ( root//TRIM( 'Turb1' )//'_Gammablj.txt'   ))
  OPEN( unit = 1018, file = ( root//TRIM( 'Turb1' )//'_Gammabljm1.txt' ))

  OPEN( unit = 1019, file = ( root//TRIM( 'Turb1' )//'_Gammaj.txt'   ))
  OPEN( unit = 1020, file = ( root//TRIM( 'Turb1' )//'_Gammajm1.txt' ))

  OPEN( unit = 1021, file = ( root//TRIM( 'Turb1' )//'_Gamma_nearj.txt'   ))
  OPEN( unit = 1022, file = ( root//TRIM( 'Turb1' )//'_Gamma_nearjm1.txt' ))

  OPEN( unit = 1023, file = ( root//TRIM( 'Turb1' )//'_BladeQuarterChordjm1.txt' ))
  OPEN( unit = 1024, file = ( root//TRIM( 'Turb1' )//'_BladeQuarterChordjm2.txt' ))

  PRINT *, 'NB is', NumBl

  READ ( 1000, * ) ((( FWake%r_primej(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1001, * ) ((( FWake%r_primejm1( j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1002, * ) ((( FWake%r_primejm2( j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1003, * ) ((( FWake%r_primejm3( j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1025, * ) ((( FWake%r_oldj(     j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1004, * ) ((( FWake%r_oldjm1(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1005, * ) ((( FWake%r_oldjm2(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1006, * ) ((( FWake%r_oldjm3(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1007, * ) ((( FWake%r_newjm1(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1008, * ) ((( FWake%r_newjm2(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1009, * ) ((( FWake%r_newjm3(   j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1010, * ) ((( FWake%rjm1(       j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )
  READ ( 1011, * ) ((( FWake%rjm2(       j2, kindx, nindx ), j2=1,3), kindx=1,CUTOFF_upinit(1)), &
     & nindx=1,NumBl )

  READ ( 1012, * ) (((( NWake%r_nearjm1( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), &
     & kindx3 = 1, NnearMax ), nindx = 1, NumBl )
  READ ( 1013, * ) (((( NWake%r_nearjm2( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), &
     & kindx3 = 1, NnearMax ), nindx = 1, NumBl )

  READ ( 1017, * ) ((NWake%Gammablj(   nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )
  READ ( 1018, * ) ((NWake%Gammabljm1( nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )

  READ ( 1019, * ) ((FWake%Gammaj(   kindx, nindx ), kindx = 1, CUTOFF_upinit(1)), nindx = 1, NumBl )
  READ ( 1020, * ) ((FWake%Gammajm1( kindx, nindx ), kindx = 1, CUTOFF_upinit(1)), nindx = 1, NumBl )

  READ ( 1021, * ) (((NWake%Gamma_nearj(   kindx3,nbs,nindx),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )
  READ ( 1022, * ) (((NWake%Gamma_nearjm1( kindx3,nbs,nindx),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )

  READ ( 1023, * ) ((( BladeQuarterChordjm1( j2, nbs, nindx ), j2=1,3 ), nbs=1,NumBS+1 ), nindx=1,NumBl )
  READ ( 1024, * ) ((( BladeQuarterChordjm2( j2, nbs, nindx ), j2=1,3 ), nbs=1,NumBS+1 ), nindx=1,NumBl )

  CLOSE( 1001 ); CLOSE( 1002 ); CLOSE( 1003 ); CLOSE( 1004 ); CLOSE( 1005 ); CLOSE( 1006 )
  CLOSE( 1007 ); CLOSE( 1008 ); CLOSE( 1009 ); CLOSE( 1010 ); CLOSE( 1011 ); CLOSE( 1012 )
  CLOSE( 1013 ); CLOSE( 1017 ); CLOSE( 1018 )
  CLOSE( 1019 ); CLOSE( 1020 ); CLOSE( 1021 ); CLOSE( 1022 ); CLOSE( 1023 ); CLOSE( 1024 )
  CLOSE( 1000 ); CLOSE( 1025 )

  ! ****************
  !IF ( mype .EQ. 0 ) THEN
  PRINT*, 'NumWakes is: ', NumWakes
  !    IF (NumTurbs .GT. 1 ) THEN
  CUTOFF = CUTOFF_upmax
  !   ELSE
  !      CUTOFF( 1 ) = CUTOFF_upmax(1)
  !   END IF

     PRINT *, 'CUTOFF values are ', CUTOFF
  !END IF  ! mype
  !CALL MPI_BCAST( CUTOFF, NumTurbs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  ! ****************
END SUBROUTINE FVW_INITIALIZE_WAKE
! ==============================================================================





! ==============================================================================
SUBROUTINE FVWtest(J, IBlade, Initial, p, u, m, W, CLFW, VINDFW, Time )

  USE FVW_Types
  USE FVW_Parm, Only: Time_Real!, RotSpeed   !KS -- removed 6.28.19
  USE NWTC_Library

  IMPLICIT NONE

  TYPE( FVW_ParameterType ), INTENT( INOUT ) :: p
  TYPE( FVW_InputType ),     INTENT( IN    ) :: u
  TYPE( FVW_MiscVarType ),   INTENT( INOUT ) :: m

  REAL(ReKi),                INTENT(   OUT ) :: W, CLFW, VINDFW(3)
  REAL(ReKi),                INTENT( IN    ) :: Time

  LOGICAL,                   INTENT( IN    ) :: Initial

  INTEGER,                   INTENT( IN    ) :: J
  INTEGER,                   INTENT( IN    ) :: IBlade

  INTEGER :: InitVal, a


!FIXME: If Initial is only on the very first call to the routine, then the first section can be moved into Init.
!        However, if it is set once per timestep, then we should figure out how to store all the FVW_INITIALIZE_WAKE
!        data so we dont need to read files every timestep (that is hugely expensive)
  Time_Real = Time
!  PRINT*, 'RotSpeed in FVWtest: ', p%RotSpeed
  IF ( Initial .AND. J .EQ. 1 .AND. IBLADE .EQ. 1 ) THEN
     InitVal=1
     CALL FVW_READ_WAKE_PARAM( p )
!FIXME: is the initialize wake call necessary if all the data is moved to miscvars and other types accordingly?
     CALL FVW_INITIALIZE_WAKE(  )
!    PRINT*, 'a'
        CALL FVW_COMPUTE_WAKE( p, p%FVWTurbineComponents, u, m, p%FVW_Wind )
!    PRINT*, 'b'
  ELSE
!   PRINT*, 'c'
     InitVal=0

     IF (J .EQ. 1 .AND. IBLADE .EQ. 1) THEN
        m%AofA   = 0.0_ReKi
        m%W2FVW  = 0.0_ReKi
        CALL FVW_COMPUTE_WAKE( p, p%FVWTurbineComponents, u, m, p%FVW_Wind )
     ENDIF
  ENDIF

!FIXME: logic question here, why only return a single value when we have calculated the whole set???
  W      = m%W2FVW(J,IBLADE)
  CLFW   = m%CLFVW(J,Iblade)
  VINDFW = m%VINDFVW(:, J, IBlade)

END SUBROUTINE FVWtest
