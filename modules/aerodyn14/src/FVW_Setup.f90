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
  IF ( .NOT. ALLOCATED( Chord )) ALLOCATE( Chord( NElm ), RElm( NElm ), RNodes( NElm ))
  RNodes          = pFVW%RNodes
  RElm            = pFVW%RElm
  NumBl           = pFVW%NumBl
  Radius          = pFVW%Radius
  Chord           = pFVW%C
  HH              = pFVW%HH
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
SUBROUTINE FVW_INITIALIZE_WAKE(  )

  USE FVW_Parm
  USE FVW_vars
  USE AeroDyn14_Types, Only: FVW_ParameterType, AD14_InputType
  USE FVW_ComputeWake, Only: BladeQuarterChordjm1, BladeQuarterChordjm2
  USE MultTurb_Params
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
        & AofA( NELM, NumBl ), W2FVW( NELM, NumBl ), CLFVW( NELM, NumBl ), VINDFVW( 3, NELM, NumBl ), &
        & FVW_Velocity( 3, NELM, NumBl ), loop( NumTurbs ))
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

  !IF (.NOT. ALLOCATED( GCoord%rj )) THEN
  !   ALLOCATE( GCoord%rj( 3, CUTOFF_Allocate, NumWakes ), GCoord%rjm1( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%rjm2( 3, CUTOFF_Allocate, NumWakes ), GCoord%r_primej( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_primejm1( 3, CUTOFF_Allocate, NumWakes ),GCoord%r_primejm2( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_primejm3( 3, CUTOFF_Allocate, NumWakes ), GCoord%r_oldj( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_oldjm1( 3, CUTOFF_Allocate, NumWakes ), GCoord%r_newj( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_oldjm2( 3, CUTOFF_Allocate, NumWakes ), GCoord%r_oldjm3( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_newjm1( 3, CUTOFF_Allocate, NumWakes ), GCoord%r_newjm2( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%r_newjm3( 3, CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%VinducedFarWakej( 3, CUTOFF_Allocate, NumBl, WakeAgeLimit, NumWakes ), &
  !      & GCoord%VinducedFarWakejm1( 3, CUTOFF_Allocate, NumBl, WakeAgeLimit, NumWakes ), &
  !      & GCoord%Gammaj( CUTOFF_Allocate, NumWakes ), GCoord%Gammajp1( CUTOFF_Allocate, NumWakes ), &
  !      & GCoord%Gammajm1( CUTOFF_Allocate, NumWakes ), &

  !      & GCoord%r_nearj( 3, NumBS+1, NnearMax, NumWakes ), GCoord%Gammablj( NumBS, NumWakes ), &
  !      & GCoord%r_nearjm1( 3, NumBS+1, NnearMax, NumWakes ), GCoord%Gammabljm1( NumBS, NumWakes ), &
  !      & GCoord%r_nearjm2( 3, NumBS+1, NnearMax, NumWakes ), &
  !      & GCoord%VinducedFarWakeRj( 3, NumBS, NumWakes, WakeAgeLimit, NumWakes ), &
  !      & GCoord%VinducedFarWakeRjm1( 3, NumBS, NumWakes, WakeAgeLimit, NumWakes ), &
  !      & GCoord%Gamma_nearj( NnearMax, NumBS+1, NumWakes ), &
  !      & GCoord%Gamma_nearjp1( NnearMax, NumBS+1, NumWakes ), &
  !      & GCoord%Gamma_nearjm1( NnearMax, NumBS+1, NumWakes ), &

  !      & GCoord%BladeThreeQuarterChordj( 3, NumBS, NumBl ), GCoord%BladeQuarterChordj( 3, NumBS+1, NumWakes ), &
  !      & GCoord%BladeQuarterChordjm1( 3, NumBS+1, NumWakes ), GCoord%BladeQuarterChordjm2( 3, NumBS+1, NumWakes ), &
  !      & GCoord%BladeLoc2j( 3, NumBS, NumBl ), GCoord%BladeLoc2j_Real( 3, NumBS, NumBl ))
  !END IF

  !GCoord%rj = 0.00_ReKi; GCoord%rjm1 = 0.00_ReKi; GCoord%rjm2 = 0.00_ReKi; GCoord%r_nearj = 0.00_ReKi; GCoord%r_nearjm1 = 0.00_ReKi
  !GCoord%r_primej = 0.00_ReKi; GCoord%r_primejm1 = 0.00_ReKi; GCoord%r_primejm2 = 0.00_ReKi; GCoord%r_primejm3 = 0.00_ReKi
  !GCoord%r_oldj = 0.00_ReKi; GCoord%r_oldjm1 = 0.00_ReKi; GCoord%r_oldjm2 = 0.00_ReKi; GCoord%r_oldjm3 = 0.00_ReKi
  !GCoord%r_newj = 0.00_ReKi; GCoord%r_newjm1 = 0.00_ReKi; GCoord%r_newjm2 = 0.00_ReKi; GCoord%r_newjm3 = 0.00_ReKi
  !GCoord%VinducedFarWakeRj = 0.00_ReKi; GCoord%VinducedFarWakeRjm1 = 0.00_ReKi
  !GCoord%VinducedFarWakej = 0.00_ReKi; GCoord%VinducedFarWakejm1 = 0.00_ReKi
  !GCoord%Gammablj = 0.00_ReKi; GCoord%Gammaj = 0.00_ReKi; GCoord%Gamma_nearj = 0.00_ReKi
  !GCoord%Gamma_nearjp1 = 0.00_ReKi; GCoord%Gammajp1 = 0.00_ReKi; GCoord%Gammabljm1 = 0.00_ReKi; GCoord%Gammajm1 = 0.00_ReKi
  !GCoord%Gamma_nearjm1 = 0.00_ReKi; GCoord%BladeThreeQuarterChordj= 0.00_ReKi
  !GCoord%BladeQuarterChordj = 0.00_ReKi; GCoord%BladeQuarterChordjm1 = 0.00_ReKi; GCoord%BladeQuarterChordjm2 = 0.00_ReKi
  !GCoord%BladeLoc2j = 0.00_ReKi; GCoord%BladeLoc2j_Real = 0.00_ReKi

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
SUBROUTINE FVW_Terminate()

  USE FVW_vars
  USE FVW_Parm
  USE FVW_ComputeWake
  USE MultTurb_Params

!FVW_vars
  IF ( ALLOCATED( FWake%r_primej                )) DEALLOCATE( FWake%r_primej                )
  IF ( ALLOCATED( FWake%r_primejm1              )) DEALLOCATE( FWake%r_primejm1              )
  IF ( ALLOCATED( FWake%r_primejm2              )) DEALLOCATE( FWake%r_primejm2              )
  IF ( ALLOCATED( FWake%r_primejm3              )) DEALLOCATE( FWake%r_primejm3              )
  IF ( ALLOCATED( FWake%r_oldj                  )) DEALLOCATE( FWake%r_oldj                  )
  IF ( ALLOCATED( FWake%r_oldjm1                )) DEALLOCATE( FWake%r_oldjm1                )
  IF ( ALLOCATED( FWake%r_oldjm2                )) DEALLOCATE( FWake%r_oldjm2                )
  IF ( ALLOCATED( FWake%r_oldjm3                )) DEALLOCATE( FWake%r_oldjm3                )
  IF ( ALLOCATED( FWake%r_newj                  )) DEALLOCATE( FWake%r_newj                  )
  IF ( ALLOCATED( FWake%r_newjm1                )) DEALLOCATE( FWake%r_newjm1                )
  IF ( ALLOCATED( FWake%r_newjm2                )) DEALLOCATE( FWake%r_newjm2                )
  IF ( ALLOCATED( FWake%r_newjm3                )) DEALLOCATE( FWake%r_newjm3                )
  IF ( ALLOCATED( FWake%rj                      )) DEALLOCATE( FWake%rj                      )
  IF ( ALLOCATED( FWake%rjm1                    )) DEALLOCATE( FWake%rjm1                    )
  IF ( ALLOCATED( FWake%rjm2                    )) DEALLOCATE( FWake%rjm2                    )
  IF ( ALLOCATED( NWake%r_nearj                 )) DEALLOCATE( NWake%r_nearj                 )
  IF ( ALLOCATED( NWake%r_nearjm1               )) DEALLOCATE( NWake%r_nearjm1               )
  IF ( ALLOCATED( NWake%r_nearjm2               )) DEALLOCATE( NWake%r_nearjm2               )
  IF ( ALLOCATED( FWake%Gammaj                  )) DEALLOCATE( FWake%Gammaj                  )
  IF ( ALLOCATED( FWake%Gammajm1                )) DEALLOCATE( FWake%Gammajm1                )
  IF ( ALLOCATED( FWake%Gammajp1                )) DEALLOCATE( FWake%Gammajp1                )
  IF ( ALLOCATED( NWake%Gammablj                )) DEALLOCATE( NWake%Gammablj                )
  IF ( ALLOCATED( NWake%Gammabljm1              )) DEALLOCATE( NWake%Gammabljm1              )
  IF ( ALLOCATED( NWake%Gamma_nearj             )) DEALLOCATE( NWake%Gamma_nearj             )
  IF ( ALLOCATED( NWake%Gamma_nearjm1           )) DEALLOCATE( NWake%Gamma_nearjm1           )
  IF ( ALLOCATED( NWake%Gamma_nearjp1           )) DEALLOCATE( NWake%Gamma_nearjp1           )
  IF ( ALLOCATED( BladeQuarterChordjm1          )) DEALLOCATE( BladeQuarterChordjm1          )
  IF ( ALLOCATED( BladeQuarterChordjm2          )) DEALLOCATE( BladeQuarterChordjm2          )
  IF ( ALLOCATED( AofA                          )) DEALLOCATE( AofA                          )
  IF ( ALLOCATED( W2FVW                         )) DEALLOCATE( W2FVW                         )
  IF ( ALLOCATED( CLFVW                         )) DEALLOCATE( CLFVW                         )
  IF ( ALLOCATED( GammaBEM                      )) DEALLOCATE( GammaBEM                      )
  IF ( ALLOCATED( PhiOLD                        )) DEALLOCATE( PhiOLD                        )
  IF ( ALLOCATED( VINDFVW                       )) DEALLOCATE( VINDFVW                       )
  IF ( ALLOCATED( FVW_Velocity                  )) DEALLOCATE( FVW_Velocity                  )

!KS -- for Parallel Code

  ! ****************
  !IF ( ALLOCATED( GCoord%BladeThreeQuarterChordj ))   DEALLOCATE( GCoord%BladeThreeQuarterChordj )
  !IF ( ALLOCATED( GCoord%BladeQuarterChordj      ))   DEALLOCATE( GCoord%BladeQuarterChordj      )
  !IF ( ALLOCATED( GCoord%BladeQuarterChordjm1    ))   DEALLOCATE( GCoord%BladeQuarterChordjm1    )
  !IF ( ALLOCATED( GCoord%BladeQuarterChordjm2    ))   DEALLOCATE( GCoord%BladeQuarterChordjm2    )

  !IF ( ALLOCATED( GCoord%VinducedFarWakejm1   ))   DEALLOCATE( GCoord%VinducedFarWakejm1   )
  !IF ( ALLOCATED( GCoord%VinducedFarWakeRjm1     ))   DEALLOCATE( GCoord%VinducedFarWakeRjm1     )
  !IF ( ALLOCATED( GCoord%VinducedFarWakej     ))   DEALLOCATE( GCoord%VinducedFarWakej     )
  !IF ( ALLOCATED( GCoord%VinducedFarWakeRj       ))   DEALLOCATE( GCoord%VinducedFarWakeRj       )
  !IF ( ALLOCATED( GCoord%r_primej                ))   DEALLOCATE( GCoord%r_primej                )
  !IF ( ALLOCATED( GCoord%r_primejm1              ))   DEALLOCATE( GCoord%r_primejm1              )
  !IF ( ALLOCATED( GCoord%r_primejm2              ))   DEALLOCATE( GCoord%r_primejm2              )
  !IF ( ALLOCATED( GCoord%r_primejm3              ))   DEALLOCATE( GCoord%r_primejm3              )
  !IF ( ALLOCATED( GCoord%r_oldj                  ))   DEALLOCATE( GCoord%r_oldj                  )
  !IF ( ALLOCATED( GCoord%r_oldjm1                ))   DEALLOCATE( GCoord%r_oldjm1                )
  !IF ( ALLOCATED( GCoord%r_oldjm2                ))   DEALLOCATE( GCoord%r_oldjm2                )
  !IF ( ALLOCATED( GCoord%r_oldjm3                ))   DEALLOCATE( GCoord%r_oldjm3                )
  !IF ( ALLOCATED( GCoord%r_newj                  ))   DEALLOCATE( GCoord%r_newj                  )
  !IF ( ALLOCATED( GCoord%r_newjm1                ))   DEALLOCATE( GCoord%r_newjm1                )
  !IF ( ALLOCATED( GCoord%r_newjm2                ))   DEALLOCATE( GCoord%r_newjm2                )
  !IF ( ALLOCATED( GCoord%r_newjm3                ))   DEALLOCATE( GCoord%r_newjm3                )
  !IF ( ALLOCATED( GCoord%rj                      ))   DEALLOCATE( GCoord%rj                      )
  !IF ( ALLOCATED( GCoord%rjm1                    ))   DEALLOCATE( GCoord%rjm1                    )
  !IF ( ALLOCATED( GCoord%rjm2                    ))   DEALLOCATE( GCoord%rjm2                    )
  !IF ( ALLOCATED( GCoord%r_nearj                 ))   DEALLOCATE( GCoord%r_nearj                 )
  !IF ( ALLOCATED( GCoord%r_nearjm1               ))   DEALLOCATE( GCoord%r_nearjm1               )
  !IF ( ALLOCATED( GCoord%r_nearjm2               ))   DEALLOCATE( GCoord%r_nearjm2               )
  !IF ( ALLOCATED( GCoord%Gammaj                  ))   DEALLOCATE( GCoord%Gammaj                  )
  !IF ( ALLOCATED( GCoord%Gammajm1                ))   DEALLOCATE( GCoord%Gammajm1                )
  !IF ( ALLOCATED( GCoord%Gammajp1                ))   DEALLOCATE( GCoord%Gammajp1                )
  !IF ( ALLOCATED( GCoord%Gammablj                ))   DEALLOCATE( GCoord%Gammablj                )
  !IF ( ALLOCATED( GCoord%Gammabljm1              ))   DEALLOCATE( GCoord%Gammabljm1              )
  !IF ( ALLOCATED( GCoord%Gamma_nearj             ))   DEALLOCATE( GCoord%Gamma_nearj             )
  !IF ( ALLOCATED( GCoord%Gamma_nearjm1           ))   DEALLOCATE( GCoord%Gamma_nearjm1           )
  !IF ( ALLOCATED( GCoord%Gamma_nearjp1           ))   DEALLOCATE( GCoord%Gamma_nearjp1           )
  ! ****************

  IF ( ALLOCATED( Chord )) DEALLOCATE( Chord )
  IF ( ALLOCATED( RElm  )) DEALLOCATE( RElm  )

!FVW_params

!No allocatable parameters

!FVW_COMPUTE_WAKE

  IF ( ALLOCATED( C1                      )) DEALLOCATE( C1                      )
  IF ( ALLOCATED( Velsec                  )) DEALLOCATE( Velsec                  )
  IF ( ALLOCATED( Velsec2                 )) DEALLOCATE( Velsec2                 )
  IF ( ALLOCATED( C2                      )) DEALLOCATE( C2                      )
  IF ( ALLOCATED( velstorej               )) DEALLOCATE( velstorej               )
  IF ( ALLOCATED( a_of_a_storej           )) DEALLOCATE( a_of_a_storej           )
  IF ( ALLOCATED( a_of_a_effective        )) DEALLOCATE( a_of_a_effective        )

  IF ( ALLOCATED( BladeTanVectj           )) DEALLOCATE( BladeTanVectj           )
  IF ( ALLOCATED( BladeQuarterChordj      )) DEALLOCATE( BladeQuarterChordj      )
  IF ( ALLOCATED( BladeLocj               )) DEALLOCATE( BladeLocj               )
  IF ( ALLOCATED( BladeNormVect2j         )) DEALLOCATE( BladeNormVect2j         )
  IF ( ALLOCATED( BladeTanVect2j          )) DEALLOCATE( BladeTanVect2j          )
  IF ( ALLOCATED( Vind_storej             )) DEALLOCATE( Vind_storej             )
  IF ( ALLOCATED( BladeLoc2j              )) DEALLOCATE( BladeLoc2j              )
  IF ( ALLOCATED( BladeThreeQuarterChordj )) DEALLOCATE( BladeThreeQuarterChordj )
  IF ( ALLOCATED( VinducedNWFinal         )) DEALLOCATE( VinducedNWFinal         )
  IF ( ALLOCATED( VinducedNWtest          )) DEALLOCATE( VinducedNWtest          )

  IF ( ALLOCATED( VinducedFW1             )) DEALLOCATE( VinducedFW1             )
  IF ( ALLOCATED( VinducedNW1             )) DEALLOCATE( VinducedNW1             )
  IF ( ALLOCATED( VinducedBC1             )) DEALLOCATE( VinducedBC1             )
  IF ( ALLOCATED( VinducedTot1            )) DEALLOCATE( VinducedTot1            )
  IF ( ALLOCATED( VinducedTot2            )) DEALLOCATE( VinducedTot2            )
  IF ( ALLOCATED( VinducedTot1b           )) DEALLOCATE( VinducedTot1b           )
  IF ( ALLOCATED( VinducedTot2b           )) DEALLOCATE( VinducedTot2b           )

  IF ( ALLOCATED( rpcyl1 )) DEALLOCATE( rpcyl1 )
  IF ( ALLOCATED( rpcyl2 )) DEALLOCATE( rpcyl2 )
  IF ( ALLOCATED( rpcyl3 )) DEALLOCATE( rpcyl3 )
  IF ( ALLOCATED( rpcyl4 )) DEALLOCATE( rpcyl4 )
  IF ( ALLOCATED( rpcyl5 )) DEALLOCATE( rpcyl5 )
  IF ( ALLOCATED( rpcyl6 )) DEALLOCATE( rpcyl6 )
  IF ( ALLOCATED( rpcyl7 )) DEALLOCATE( rpcyl7 )
  IF ( ALLOCATED( rpcyl8 )) DEALLOCATE( rpcyl8 )
  IF ( ALLOCATED( rncyl1 )) DEALLOCATE( rncyl1 )
  IF ( ALLOCATED( rncyl2 )) DEALLOCATE( rncyl2 )
  IF ( ALLOCATED( rncyl3 )) DEALLOCATE( rncyl3 )
  IF ( ALLOCATED( rncyl4 )) DEALLOCATE( rncyl4 )
  IF ( ALLOCATED( rncyl5 )) DEALLOCATE( rncyl5 )
  IF ( ALLOCATED( rncyl6 )) DEALLOCATE( rncyl6 )
  IF ( ALLOCATED( rncyl7 )) DEALLOCATE( rncyl7 )
  IF ( ALLOCATED( rncyl8 )) DEALLOCATE( rncyl8 )
  IF ( ALLOCATED( rocyl1 )) DEALLOCATE( rocyl1 )
  IF ( ALLOCATED( rocyl2 )) DEALLOCATE( rocyl2 )
  IF ( ALLOCATED( rocyl3 )) DEALLOCATE( rocyl3 )
  IF ( ALLOCATED( rocyl4 )) DEALLOCATE( rocyl4 )
  IF ( ALLOCATED( rocyl5 )) DEALLOCATE( rocyl5 )
  IF ( ALLOCATED( rocyl6 )) DEALLOCATE( rocyl6 )
  IF ( ALLOCATED( rocyl7 )) DEALLOCATE( rocyl7 )
  IF ( ALLOCATED( rocyl8 )) DEALLOCATE( rocyl8 )

  IF ( ALLOCATED( CalcedVinf      )) DEALLOCATE( CalcedVinf      )
  IF ( ALLOCATED( VindTotal       )) DEALLOCATE( VindTotal       )

  IF ( ALLOCATED( Vaxial2j        )) DEALLOCATE( Vaxial2j        )
  IF ( ALLOCATED( VNElem2j        )) DEALLOCATE( VNElem2j        )
  IF ( ALLOCATED( Vaxialj         )) DEALLOCATE( Vaxialj         )
  IF ( ALLOCATED( VTT             )) DEALLOCATE( VTT             )
  IF ( ALLOCATED( VNElement       )) DEALLOCATE( VNElement       )

  IF ( ALLOCATED( BladeLoc2j_Real )) DEALLOCATE( BladeLoc2j_Real )
  IF ( ALLOCATED( r_oldj_Real     )) DEALLOCATE( r_oldj_Real     )
  IF ( ALLOCATED( r_primej_Real   )) DEALLOCATE( r_primej_Real   )

  IF ( ALLOCATED( CUTOFF          )) DEALLOCATE( CUTOFF          )
  IF ( ALLOCATED( CUTOFF_upinit   )) DEALLOCATE( CUTOFF_upinit   )
  IF ( ALLOCATED( CUTOFF_up       )) DEALLOCATE( CUTOFF_up       )
  IF ( ALLOCATED( CUTOFF_upmax    )) DEALLOCATE( CUTOFF_upmax    )

END SUBROUTINE FVW_Terminate
! ==============================================================================
