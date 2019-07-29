

SUBROUTINE FVW_COMPUTE_WAKE( TurbineComponents, InputMarkers, Wind_FVW )


  USE AeroDyn14_Types
  USE AD14AeroConf_Types
  USE FVW_vars
  USE FVW_Parm
  USE FVW_ComputeWake
  USE NWTC_Library
  USE MultTurb_Params
  USE InflowWind
  USE FileManipulation, Only: OutputFinalWake

  IMPLICIT NONE

  !INCLUDE 'mpif.h'

  TYPE( AD14AeroConf_InputType ),   INTENT( IN    ) :: TurbineComponents
  TYPE( MeshType ), DIMENSION(NumBl), INTENT( IN    ) :: InputMarkers
  TYPE( FVW_WindType ), INTENT( INOUT ) :: Wind_FVW

  INTEGER( IntKi ) :: nbs, j, k, n, indx, ErrStat, q, ierr, Tlow, Thigh, StartCount, ProcNum, kindx
  CHARACTER(124) :: ErrorMsg

  REAL( ReKi ) :: Cap_Gamma, RMSval, zloc, Period
  REAL( ReKi ), DIMENSION( 3 )            :: tmpvector
  REAL( ReKi ), DIMENSION( NumBS )        :: Rnumbs, xnumbs
  REAL( ReKi ), DIMENSION( NumBS + 1 )    :: Rnumbsp1, xnumbsp1, tmp1
  REAL( ReKi ), DIMENSION( 3, NumBS )     :: VTotal
  REAL( ReKi ), DIMENSION( NumBS, NumBl ) :: velstorej2, cl_storej


!These are temporary variables used for collecting/distributing global wake info before/after pred/corr step
  INTEGER( IntKi ) :: Size_Near, Size_NearV, Size_NearGBl, Size_NearG, Size_NearB, Size_Far, Size_FarV, Size_FarG, llow, lhigh, glow, ghigh
  REAL(ReKi), DIMENSION( :, : ),          ALLOCATABLE :: Gammaj_tmp
  REAL(ReKi), DIMENSION( :, : ),          ALLOCATABLE :: Gammablj_tmp
  REAL(ReKi), DIMENSION( :, :, : ),       ALLOCATABLE :: Gamma_nearj_tmp
  REAL(ReKi), DIMENSION( :, :, : ),       ALLOCATABLE :: r_primej_tmp
  REAL(ReKi), DIMENSION( :, :, :, : ),    ALLOCATABLE :: r_nearj_tmp
  REAL(ReKi), DIMENSION( :, :, : ),       ALLOCATABLE :: BladeQuarterChordj_tmp

  LOGICAL :: OutWake
  REAL :: WhichTurb
  tmpvector=0.0_ReKi
PRINT*, "Before SetupWake"
  CALL SetupWake( )
PRINT*, "After SetupWake"

  ! Set all the blade locations and blade tangent and normal vectors
  CALL BladeLocations( Vaxial2j, VNElem2j )
PRINT*, "After BladeLocations"
  !sets the BC for FW
  FWake%r_newj( :, 1, : ) = FWake%r_primej( :, 1, : )
  FWake%rj(     :, 1, : ) = FWake%r_newj(   :, 1, : )

  !for j = 1, computes the angle of attack of the blades given the rotational
  !and wind speeds

  loop = 1
  counter = 0
  DO WHILE ( ANY( loop .EQ. 1 ))
     Size_Near    = NumBl * NnearMax * ( NumBS+1 ) * 3  ! for location of near wake markers
!For # of wakes ^ For Near wake^^^^^^^^^       ^For wake position

     Size_NearG   = NumBl * ( NumBS+1 ) * NnearMax  ! for circulation of near wake
     Size_Far = NumBl * CUTOFF_Allocate * 3
     Size_NearV   = NumBl * WakeAgeLimit * NumBl * NumBS * 3  ! for velocity of near wake markers
!For # of wakes ^   For Far wake^^     ^For wake position
     Size_FarG = NumBl * CUTOFF_Allocate

     !!!This block loops over the local turbine blades and computes the induced velocity on each blade section (nbs) due to the far wake. This induced velocity (VTotal) is then used to compute the ciculation on the blades (Gammaj) as well as near wake marker locations (r_nearj). While this computation is for local turbine blades, the Vinduced considers the markers from EVERY turbine, which is why MPI_BCAST is necessary.
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
     !CALL TransformToGlobal( )
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
     !IF( NumTurbs .GT. 1 ) THEN
     !   DO ProcNum = 0, NumTurbs-1
     !      StartCount = ProcNum*NumBl+1
           !CALL MPI_BCAST( GCoord%rjm1( 1, 1, StartCount ), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
           !CALL MPI_BCAST( GCoord%rjm2( 1, 1, StartCount ), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
           !CALL MPI_BCAST( GCoord%Gammajm1( 1, StartCount ), Size_FarG, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !   END DO
     !END IF
     DO n = 1, NumBl      !Loop over # of blades
        VTotal    = 0.0_ReKi
        VindTotal = 0.0_ReKi

        WhichTurb = REAL(n-0.01)/REAL(NumBl)
        q = n + FLOOR( WhichTurb )*NumBl
        DO nbs = Num_start, NumBS      !Loops over blade segments
           CALL Vinduced3( FWake%rjm1, FWake%Gammajm1, BladeThreeQuarterChordj( :, nbs, n ), &
              & FWake%rjm2, q, j, nbs ) !KS -- Have to keep rjm1 and rjm2 options
  !NWake%VinducedFarWakeRj is calculated in Vinduced3; I don't think I have to pass these values between nodes.
           !adds the induced velocity at the blade from the 3 Time steps of near wake

           VindTotal( 1, nbs ) = Sum( FWake%VinducedFarWakeRj( 1, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
           VindTotal( 2, nbs ) = Sum( FWake%VinducedFarWakeRj( 2, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
           VindTotal( 3, nbs ) = Sum( FWake%VinducedFarWakeRj( 3, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))

           tmpvector = BladeLoc2j_Real( :, nbs, n )

           CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )
           Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector
           CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )

           CalcedVinf = Wind_FVW%OutputData%VelocityUVW( :, 1 )
           CALL TRANSFORM_TO_FVW_COORDS( CalcedVinf )
           VTotal( :, nbs ) = CalcedVinf + VindTotal( :, nbs ) + Vaxial2j( :, nbs, n ) + &
              & VNElem2j( :, nbs, n )

           VN  = DOT_PRODUCT( VTotal( :, nbs ), BladeNormVect2j( :, nbs, n ) )
           VT  = DOT_PRODUCT( VTotal( :, nbs ), BladeTanVect2j(  :, nbs, n ) )

           a_of_a_effective( nbs, n ) = atan( VN / VT ) !****NEED TO CHANGE ****
           Velsec( nbs ) = sqrt( VT * VT + VN * VN )
        ENDDO  ! End of Loop over blade segments

        Velstorej( :, n ) = Velsec !****NEED TO CHANGE ****
        ! ************
        ! KS -- can keep this inside parallel portion, but need to use LOCAL indexing.
        ! ************
        tmp1 = NWake%Gamma_nearj( 1, :, n )
        CALL Calculate_Gamma1( n, VTotal, BladeTanVectj( :, :, n ), BladeNormVect2j( :, :, n ), &
                             & BladeQuarterChordj, BladeThreeQuarterChordj( :, :, n ), &
                             & Cap_Gamma, NWake%Gammablj( :, n ), NWake%r_nearjm1( :, :, :, n ), &
                             & NWake%r_nearj( :, :, :, n ), tmp1, zloc, &!, NWake%Gamma_nearj( 1, :, n ), zloc, &
                             & VinducedNWFinal( :, :, n ), Wind_FVW )
        FWake%Gammaj( 1, n ) = Cap_Gamma

     ENDDO

     !Update rnear values from Calculate_Gamma Results
     FWake%rj(       :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_primej( :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_oldj(   :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_newj(   :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )

     counter = counter + 1

    ! *******************
    ! KS -- Convert to GLOBAL COORDINATES
    ! *******************

     !CALL TransformToGlobal(  )

    ! *******************
    ! KS -- Merge wake locations; do this on ONE thread
    ! *******************

     Size_Near    = NumBl * NnearMax * ( NumBS+1 ) * 3  ! for location of near wake markers
!For # of wakes ^ For Near wake^^^^^^^^^       ^For wake position
     Size_NearV   = NumBl * WakeAgeLimit * NumBl * NumBS * 3  ! for velocity of near wake markers
     Size_NearB = 3 * ( NumBS+1 ) * NumBl  !For blades
     Size_NearGBl = NumBl * NumBS  ! For circulation of blades
     Size_NearG   = NumBl * ( NumBS+1 ) * NnearMax  ! for circulation of near wake

     Size_Far  = NumBl * CUTOFF_Allocate * 3  ! for location of far wake markers
!For # of wakes ^   For Far wake^^     ^For wake position
     Size_FarV = NumBl * WakeAgeLimit * NumBl * CUTOFF_Allocate * 3  ! for velocity of far wake markers
     Size_FarG = NumBl * CUTOFF_Allocate  ! for circulation of far wake markers

! Allocate temporary variable used for send/receive of global wake data
     IF ( NumTurbs .GT. 1 ) THEN
        IF ( .NOT. ALLOCATED( r_primej_tmp)) ALLOCATE( r_primej_tmp( 3, CUTOFF_Allocate, NumBl ), &
            & Gammaj_tmp( CUTOFF_Allocate, NumBl ), &
           & r_nearj_tmp( 3, NumBS+1, NnearMax, NumBl ), Gammablj_tmp( NumBS, NumBl ), &
           & Gamma_nearj_tmp( NnearMax, NumBS+1, NumBl ), &
           & BladeQuarterChordj_tmp( 3, NumBS+1, NumBl ))

        Tlow = NumBl*(NTurb-1)+1; Thigh = ( NTurb*NumBl )
        r_primej_tmp = FWake%r_primej( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_primej(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_primejm1( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_primejm1(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_primejm2( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_primejm2(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_primejm3( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_primejm3(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_oldj( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_oldj(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_oldjm1( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_oldjm1(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_oldjm2( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_oldjm2(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_oldjm3( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_oldjm3(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_newj( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_newj(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_newjm1( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_newjm1(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_newjm2( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_newjm2(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_primej_tmp = FWake%r_newjm3( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, &
        !            & GCoord%r_newjm3(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

        Gammaj_tmp = FWake%Gammaj( :, Tlow:Thigh )
        !CALL MPI_GATHER( Gammaj_tmp(1,1), Size_FarG, MPI_REAL8, &
        !            & GCoord%Gammaj(1,1), Size_FarG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

        r_nearj_tmp = NWake%r_nearj( :, :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, &
        !            & GCoord%r_nearj(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_nearj_tmp = NWake%r_nearjm1( :, :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, &
        !            & GCoord%r_nearjm1(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        r_nearj_tmp = NWake%r_nearjm2( :, :, :, Tlow:Thigh )
        !CALL MPI_GATHER( r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, &
        !            & GCoord%r_nearjm2(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        Gamma_nearj_tmp = NWake%Gamma_nearj( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( Gamma_nearj_tmp(1,1,1), Size_NearG, MPI_REAL8, &
        !            & GCoord%Gamma_nearj(1,1,1), Size_NearG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        Gamma_nearj_tmp = NWake%Gamma_nearjm1( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( Gamma_nearj_tmp(1,1,1), Size_NearG, MPI_REAL8, &
        !            & NWake%Gamma_nearjm1(1,1,1), Size_NearG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        Gammablj_tmp = NWake%Gammablj( :, Tlow:Thigh )
        !CALL MPI_GATHER( Gammablj_tmp(1,1), Size_NearGBl, MPI_REAL8, &
        !            & GCoord%Gammablj(1,1), Size_NearGBl, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        Gammablj_tmp = NWake%Gammabljm1( :, Tlow:Thigh )
        !CALL MPI_GATHER( Gammablj_tmp(1,1), Size_NearGBl, MPI_REAL8, &
        !            & GCoord%Gammabljm1(1,1), Size_NearGBl, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )

        BladeQuarterChordj_tmp = BladeQuarterChordj( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( BladeQuarterChordj_tmp(1,1,1), Size_NearB, MPI_REAL8, &
        !            & GCoord%BladeQuarterChordj(1,1,1), Size_NearB, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        BladeQuarterChordj_tmp = BladeQuarterChordjm1( :, :, Tlow:Thigh )
        !CALL MPI_GATHER( BladeQuarterChordj_tmp(1,1,1), Size_NearB, MPI_REAL8, &
        !            & GCoord%BladeQuarterChordjm1(1,1,1), Size_NearB, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

        DEALLOCATE( r_primej_tmp, Gammaj_tmp, r_nearj_tmp, Gammablj_tmp, Gamma_nearj_tmp, &
               & BladeQuarterChordj_tmp)

        !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
     END IF
  ! ****************
     !IF ( mype .EQ. 0 ) THEN
     DO k = NnearMax + 1, MAXVAL(CUTOFF)
        CALL Predictor( k )
        CALL Corrector( k )
     END DO
     !END IF  ! mype

    ! *******************
    ! KS -- Split wake locations; do this on ALL threads
    ! *******************

! Allocate temporary variable used for send/receive of global wake data
     IF (NumTurbs .GT. 1 ) THEN
        IF ( .NOT. ALLOCATED( r_primej_tmp)) ALLOCATE( r_primej_tmp( 3, CUTOFF_Allocate, NumBl ), &
            & Gammaj_tmp( CUTOFF_Allocate, NumBl ), &
           & r_nearj_tmp( 3, NumBS+1, NnearMax, NumBl ), Gammablj_tmp( NumBS, NumBl ), &
           & Gamma_nearj_tmp( NnearMax, NumBS+1, NumBl ))

        Tlow = NumBl*(NTurb-1)+1; Thigh = ( NTurb*NumBl ) !This is probably redundant.
        r_primej_tmp = FWake%r_primej( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_primej(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_primej( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_primejm1( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_primejm1(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_primejm1( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_primejm2( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_primejm2(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_primejm2( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_primejm3( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_primejm3(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_primejm3( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_oldj( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_oldj(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_oldj( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_oldjm1( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_oldjm1(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_oldjm1( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_oldjm2( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_oldjm2(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_oldjm2( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_oldjm3( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_oldjm3(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_oldjm3( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_newj( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_newj(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_newj( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_newjm1( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_newjm1(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_newjm1( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_newjm2( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_newjm2(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_newjm2( :, :, Tlow:Thigh ) = r_primej_tmp

        r_primej_tmp = FWake%r_newjm3( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_newjm3(1,1,1), Size_Far, MPI_REAL8, &
        !             & r_primej_tmp(1,1,1), Size_Far, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%r_newjm3( :, :, Tlow:Thigh ) = r_primej_tmp

        Gammaj_tmp = FWake%Gammaj( :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%Gammaj(1,1), Size_FarG, MPI_REAL8, &
        !             & Gammaj_tmp(1,1), Size_FarG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        FWake%Gammaj( :, Tlow:Thigh ) = Gammaj_tmp

        r_nearj_tmp = NWake%r_nearj( :, :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_nearj(1,1,1,1), Size_Near, MPI_REAL8, &
        !             & r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%r_nearj( :, :, :, Tlow:Thigh ) = r_nearj_tmp

        r_nearj_tmp = NWake%r_nearjm1( :, :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_nearjm1(1,1,1,1), Size_Near, MPI_REAL8, &
        !             & r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%r_nearjm1( :, :, :, Tlow:Thigh ) = r_nearj_tmp

        r_nearj_tmp = NWake%r_nearjm2( :, :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%r_nearjm2(1,1,1,1), Size_Near, MPI_REAL8, &
        !             & r_nearj_tmp(1,1,1,1), Size_Near, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%r_nearjm2( :, :, :, Tlow:Thigh ) = r_nearj_tmp

        Gamma_nearj_tmp = NWake%Gamma_nearj( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%Gamma_nearj(1,1,1), Size_NearG, MPI_REAL8, &
        !             & Gamma_nearj_tmp(1,1,1), Size_NearG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%Gamma_nearj( :, :, Tlow:Thigh ) = Gamma_nearj_tmp

        Gamma_nearj_tmp = NWake%Gamma_nearjm1( :, :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%Gamma_nearjm1(1,1,1), Size_NearG, MPI_REAL8, &
        !             & Gamma_nearj_tmp(1,1,1), Size_NearG, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%Gamma_nearjm1( :, :, Tlow:Thigh ) = Gamma_nearj_tmp

        Gammablj_tmp = NWake%Gammablj( :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%Gammablj(1,1), Size_NearGBl, MPI_REAL8, &
        !             & Gammablj_tmp(1,1), Size_NearGBl, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%Gammablj( :, Tlow:Thigh ) = Gammablj_tmp

        Gammablj_tmp = NWake%Gammabljm1( :, Tlow:Thigh )
        !CALL MPI_SCATTER( GCoord%Gammabljm1(1,1), Size_NearGBl, MPI_REAL8, &
        !             & Gammablj_tmp(1,1), Size_NearGBl, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
        NWake%Gammabljm1( :, Tlow:Thigh ) = Gammablj_tmp
        !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


        DEALLOCATE( r_primej_tmp, Gammaj_tmp, r_nearj_tmp, Gammablj_tmp,&
        & Gamma_nearj_tmp)
     END IF
    ! *******************
    ! KS -- Convert to LOCAL COORDINATES
    ! *******************
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
     !CALL TransformToLocal (  )
     !check for convergence      ***KS -- Maybe move this into rms subroutine
     RMSval = 0.0_ReKi
     CALL rms( FWake%r_newj, FWake%r_oldj, RMSval )
     IF ( RMSval .LT. eps ) THEN
        loop( NTurb ) = 0

        FWake%rj( :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
        FWake%rj( :, NnearMax+1:CUTOFF_up( NTurb ), : ) = FWake%r_newj(  :, NnearMax+1:CUTOFF_up( NTurb ), : )
     !KS -- 10.14.16  should last dimension be : or 1:NumBl???)
     ELSE
        FWake%r_oldj = FWake%r_newj
     ENDIF !RMSval

     !IF ( NumTurbs .GT. 1 ) THEN
        !DO ProcNum = 0, NumTurbs-1
           !CALL MPI_BCAST( loop( ProcNum+1 ), 1, MPI_INTEGER, ProcNum, MPI_COMM_WORLD, ierr )
        !END DO
     !END IF !NumTurbs
  ENDDO

  !MOVE THE CIRCULATION THROUGH THE VORTEX

  DO k = 2, CUTOFF_up( NTurb)
     FWake%Gammajp1( k, : ) = FWake%Gammaj( k-1, : )
  ENDDO

  DO k = 2, NnearMax
     NWake%Gamma_nearjp1( k, :, : ) = NWake%Gamma_nearj( k-1, :, : )
  ENDDO

    !CALL TransformToGlobal
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

!Need to collect appropriate rj, Gammaj, and rjm1 data here!!

     !IF( NumTurbs .GT. 1 ) THEN
     !   DO ProcNum = 0, NumTurbs-1
     !      StartCount = ProcNum*NumBl+1
     !      CALL MPI_BCAST( GCoord%rj( 1, 1, StartCount), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !      CALL MPI_BCAST( GCoord%rjm1( 1, 1, StartCount ), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !      CALL MPI_BCAST( GCoord%Gammaj( 1, StartCount ), Size_FarG, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )

     !      CALL MPI_BCAST( GCoord%r_nearj( 1, 1, 1, StartCount ), Size_Near, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr )
     !      CALL MPI_BCAST( GCoord%r_nearjm1( 1, 1, 1, StartCount ), Size_Near, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr )
     !      CALL MPI_BCAST( GCoord%Gamma_nearj( 1, 1, StartCount ), Size_NearG, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr )
     !   END DO
     !END IF
  DO n = 1, NumBl

!Need to collect appropriate rj, Gammaj, and rjm1 data here!!
     VTotal    = 0.0_ReKi
     VindTotal = 0.0_ReKi
     k = 1
     q = n! + mype*NumBl
     DO nbs = Num_start, NumBS
        CALL Vinduced3( FWake%rj, FWake%Gammaj, BladeThreeQuarterChordj( :, nbs, n ), FWake%rjm1, &
           & q, j, nbs )

        !adds the induced velocity at the blade from the far wake
        VindTotal( 1, nbs ) = Sum( FWake%VinducedFarWakeRj( 1, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
        VindTotal( 2, nbs ) = Sum( FWake%VinducedFarWakeRj( 2, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
        VindTotal( 3, nbs ) = Sum( FWake%VinducedFarWakeRj( 3, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))

        tmpvector = BladeLoc2j( :, nbs, n )

        CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )

        Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector
        CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )
        CalcedVinf = Wind_FVW%OutputData%VelocityUVW( :, 1 )
        CALL TRANSFORM_TO_FVW_COORDS( CalcedVinf )
        VTotal( :, nbs ) = CalcedVinf + VindTotal( :, nbs ) + Vaxial2j( :, nbs, n ) + &
           & VNElem2j( :, nbs, n )
     ENDDO

     CALL Calculate_Gamma1( n, VTotal, BladeTanVectj( :, :, n ), BladeNormVect2j( :, :, n ), &
                          & BladeQuarterChordj, BladeThreeQuarterChordj( :, :, n ), Cap_Gamma, &
                          & NWake%Gammablj( :, n ), NWake%r_nearjm1( :, :, :, n ), &
                          & NWake%r_nearj( :, :, :, n ), NWake%Gamma_nearj( 1, :, n ), zloc, &
                          & VinducedNWFinal( :, :, n ), Wind_FVW)

     !Update rnear values from Calculate_Gamma Results
     FWake%rj(       :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_primej( :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_oldj(   :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%r_newj(   :, 1:NnearMax, 1:NumBl ) = NWake%r_nearj( :, NumBS+1, :, : )
     FWake%Gammaj = FWake%Gammajm1
     FWake%Gammaj( 1, n ) = Cap_Gamma

!Transform just these values to Global
     !IF ( mype .EQ. 0 ) THEN
     llow = 1; lhigh = NumBl; glow = 1; ghigh = NumBl
     !ELSE
     !   llow = 1; lhigh = NumBl; ghigh = (mype+1)*NumBl; glow = ghigh-(NumBl-1)
     !END IF

     !GCoord%rj(         1, :, glow:ghigh ) = FWake%rj(         1, :, llow:lhigh )
     !GCoord%r_primej(   1, :, glow:ghigh ) = FWake%r_primej(   1, :, llow:lhigh )
     !GCoord%r_oldj(     1, :, glow:ghigh ) = FWake%r_oldj(     1, :, llow:lhigh )
     !GCoord%r_newj(     1, :, glow:ghigh ) = FWake%r_newj(     1, :, llow:lhigh )

     !GCoord%rj(         2, :, glow:ghigh ) = FWake%rj(         2, :, llow:lhigh ) + TurbLocs( NTurb, 1 )
     !GCoord%r_primej(   2, :, glow:ghigh ) = FWake%r_primej(   2, :, llow:lhigh ) + TurbLocs( NTurb, 1 )
     !GCoord%r_oldj(     2, :, glow:ghigh ) = FWake%r_oldj(     2, :, llow:lhigh ) + TurbLocs( NTurb, 1 )
     !GCoord%r_newj(     2, :, glow:ghigh ) = FWake%r_newj(     2, :, llow:lhigh ) + TurbLocs( NTurb, 1 )

     !GCoord%rj(         3, :, glow:ghigh ) = FWake%rj(         3, :, llow:lhigh ) + TurbLocs( NTurb, 2 )
     !GCoord%r_primej(   3, :, glow:ghigh ) = FWake%r_primej(   3, :, llow:lhigh ) + TurbLocs( NTurb, 2 )
     !GCoord%r_oldj(     3, :, glow:ghigh ) = FWake%r_oldj(     3, :, llow:lhigh ) + TurbLocs( NTurb, 2 )
     !GCoord%r_newj(     3, :, glow:ghigh ) = FWake%r_newj(     3, :, llow:lhigh ) + TurbLocs( NTurb, 2 )

     !GCoord%Gammaj(   :, glow:ghigh ) = FWake%Gammaj(   :, llow:lhigh )



  ENDDO
  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  !CALL TransformToLocal( )
  CALL UpdateAeroVals( )
  Period = (TwoPi_D/RotSpeed)
  INQUIRE( FILE="OutputWake.txt", EXIST=OutWake )
  IF ( (Time_Real .GT. (TMax-Period)) .AND. (IOutWake .EQ. 1) .AND. (OutWake .EQV. .TRUE.) .AND. (FWake%rjm1(2,1,1) .GT. -8.0_ReKi) .AND. (FWake%rjm1(2,1,1) .LT. 8.0_ReKi) .AND. ( FWake%rjm1(1,1,1) .GT. 0.0_ReKi )) THEN
     CALL OutputFinalWake( BladeQuarterChordjm1, BladeQuarterChordjm2 )
     IOutWake = 2
  END IF

  IF ( Time_Real .GE. TMax-(TMax/10.0_ReKi)) THEN
     TurbLength = Turblength+FWake%rjm1( 3, CUTOFF_up(1), 1 )
     counter2 = counter2+1
     IF ( Time_Real .GE. TMax-dtAero ) THEN
        PRINT*, '=============================='
        !PRINT*, 'For Turbine', mype+1
        PRINT*, 'The final wake length is: ', TurbLength/DBLE(counter2)

        !IF ( mype .EQ. 0 ) THEN
!        OPEN( unit = 7001, file = 'done.txt', STATUS = 'new' )
        !END IF
     END IF
  END IF
  !IF ( Time_Real .GE. TMax-(TMax/100.0_ReKi)) THEN
  CALL WriteWake( )
  !END IF
CONTAINS

!==========================================================================
SUBROUTINE RMS(r_new,r_old,RMSval)

  USE FVW_Parm, Only : CUTOFF_up, CUTOFF_allocate
  USE MultTurb_Params, Only: NumWakes, NTurb

  IMPLICIT NONE

  INTEGER                                       :: k,n,indx,counter
  REAL( ReKi )                                        :: RMSval, summation
  REAL( ReKi ), DIMENSION( 3, CUTOFF_allocate, NumWakes ) :: r_new, r_old

  summation = 0.0_ReKi
  counter = 0
  NTurb = 1  !! KS -- HACKHACKHACK 7.1.19
  DO n = 1, NumWakes
     DO k = 1, CUTOFF_up( NTurb )
        DO indx = 1, 3
           counter = counter + 1
           summation = summation + ( r_new( indx, k, n ) - r_old( indx, k, n ))**2.0_ReKi
        END DO
     END DO
  END DO

  RMSval = sqrt( summation ) / REAL(counter,ReKi)

END SUBROUTINE RMS
!==========================================================================
SUBROUTINE SetupWake( )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! This subroutine allocates arrays, sets up the initializes arrays to be zero and !!!
  !!!sets up the discretization values (delta_psi).                                   !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  IF ( .NOT. ALLOCATED( C1 )) ALLOCATE( C1( NumBS+1 ), Velsec( NumBS ), Velsec2( NumBS ), &
     & C2( NumBS ), velstorej( NumBS, NumBl ), a_of_a_storej( NumBS, NumBl ), &
     & a_of_a_effective( NumBS, NumBl ), BladeTanVectj( 3, NumBS+1, NumBl ), &
     & BladeQuarterChordj( 3, NumBS+1, NumBl ), BladeLocj( 3, NumBS+1, NumBl ), &
     & BladeNormVect2j( 3, NumBS, NumBl ), BladeTanVect2j( 3, NumBS, NumBl ), &
     & Vind_storej( 3, NumBS, NumBl ), BladeLoc2j( 3, NumBS, NumBl ), &
     & BladeThreeQuarterChordj( 3, NumBS, NumBl ), BladeLoc2j_Real( 3, NumBS, NumBl ), &
     & CalcedVinf( 3 ), VindTotal( 3, NumBS ), VinducedNWtest( 3 ), &
     & Vaxial2j( 3, NumBS, NumBl ), VNElem2j( 3, NumBS, NumBl ), Vaxialj( 3, NumBS+1, NumBl ), &
     & VTT( 3, NELM, NumBl ), VNElement( 3, NELM, NumBl ), VinducedNWFinal( 3, NumBS, NumWakes ))

  IF ( .NOT. ALLOCATED( VinducedFW1 )) ALLOCATE( VinducedFW1( 3 ), VinducedNW1( 3 ), VinducedBC1( 3 ), &
     & VinducedTot1( 3 ), VinducedTot2( 3 ), VinducedTot1b( 3 ), VinducedTot2b( 3 ), &
     & rpcyl1( 3 ), rpcyl2( 3 ), rpcyl3( 3 ), rpcyl4( 3 ), rpcyl5( 3 ), rpcyl6( 3 ), rpcyl7( 3 ), &
     & rpcyl8( 3 ), rncyl1( 3 ), rncyl2( 3 ), rncyl3( 3 ), rncyl4( 3 ), rncyl5( 3 ), rncyl6( 3 ), &
     & rncyl7( 3 ), rncyl8( 3 ))

  Omega = ( RotSpeed )
  Rad   = ( Radius )

  IF ( delta_psi( 1 ) .EQ. 0.0_ReKi ) THEN     !Sets initial delta_psi values
     delta_psi = DtAero * RotSpeed
  ELSE                                      !propogates delta_psi
     delta_psi( 3 ) = delta_psi( 2 )
     delta_psi( 2 ) = delta_psi( 1 )
     delta_psi( 1 ) = DtAero * RotSpeed
  ENDIF

100 FORMAT( 4F15.7 )

  !! sets up blade segments and lengths
  Num_start = INT( Root_cut * NumBS + 1 )
  dRad = Rad / REAL( NumBS,ReKi )
  dx = 1.0_ReKi / REAL( NumBS,ReKi )

  DO nbs = 1, NumBS
     Rnumbsp1( nbs ) = dRad * REAL( nbs - 1,ReKi )
     xnumbsp1( nbs ) = dx * REAL( nbs - 1,ReKi )

     Rnumbs( nbs ) = dRad * REAL( nbs - 1,ReKi ) + dRad / 2.0_ReKi
     xnumbs( nbs ) = dx * REAL( nbs - 1,ReKi ) + dx / 2.0_ReKi
  END DO ! nbs

  Rnumbsp1( NumBS+1 ) = dRad * REAL( NumBS,ReKi )
  xnumbsp1( NumBS+1 ) = dx * REAL( NumBS,ReKi )

  C2 = interpolation_array( RELM, Chord , Rnumbs,   NumBS, NELM )
  C1 = interpolation_array( RELM, Chord , Rnumbsp1, NumBS + 1, NELM )

  a_of_a_effective = 0.0_ReKi; a_of_a_storej = 0.0_ReKi; cl_storej = 0.0_ReKi; Vind_storej = 0.0_ReKi
  BladeLocj = 0.0_ReKi; BladeThreeQuarterChordj = 0.0_ReKi; BladeQuarterChordj = 0.0_ReKi
  BladeTanVectj = 0.0_ReKi; BladeNormVect2j = 0.0_ReKi; BladeTanVect2j = 0.0_ReKi; BladeLoc2j = 0.0_ReKi

  Vaxial2j = 0.0_ReKi; VNElem2j = 0.0_ReKi ; Vaxialj = 0.0_ReKi ; VNElement = 0.0_ReKi ; VTT = 0.0_ReKi
  Velstorej = 0.0_ReKi; Velstorej2 = 0.0_ReKi

  VN = 0.0_ReKi; VT = 0.0_ReKi

  VTotal = 0.0_ReKi; Velsec = 0.0_ReKi;  Cap_Gamma = 0.0_ReKi
  VinducedNWTest = 0.0_ReKi; VinducedNWFinal = 0.0_ReKi; CalcedVinf = 0.0_ReKi; VindTotal = 0.0_ReKi

  VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
  VinducedTot1 = 0.0_ReKi; VinducedTot2 = 0.0_ReKi; VinducedTot1b = 0.0_ReKi; VinducedTot2b = 0.0_ReKi

  !j = int(( 10.0_ReKi * TwoPi_D / RotSpeed + Time_Real ) / DtAero )!! KS changed
  !9.19.16 b/c Vinduced codes weren't computing Vind for all markers!! j started
  !ar 360, which was fine when #markers = 216, but that's the reason the fiament
  !length as 0 for farther waay markers and the filament length started being
  !computed as the code progressed in time.
 !  j = CUTOFF_upinit(mype+1) + Time_Real/DtAero  !10.7.16 Don't think this is right
  j = int(( 10._ReKi * TwoPi_D / RotSpeed_Est + Time_Real ) / DtAero ) !10.7.16  changed to this -- see pg. 154

  zloc = InputMarkers(1)%Position( 1, 1 )

END SUBROUTINE SetupWake
!==========================================================================

!==========================================================================
SUBROUTINE BladeLocations( VaxialBL, VNElemBL )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! This subroutine sets up the blade positions                                     !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE

  REAL( ReKi ), DIMENSION( 3, NumBS, NumBl ),  INTENT(   OUT ) :: VaxialBL, VNElemBL

  REAL( ReKi ), DIMENSION( 3 ) :: tmpvectorj, tmpvectorj2
  REAL( ReKi ),       DIMENSION( 3               ) :: V
  INTEGER :: jmax, indx1, indx2
  REAL( ReKi ) :: TMP_Vect( 3 ), dr
  REAL( ReKi ), DIMENSION(:), ALLOCATABLE :: tmp, tmp2
  REAL( ReKi ), DIMENSION(NELM) :: tmp3

  ALLOCATE(tmp(NELM-nelm_start),tmp2(NELM-nelm_start))

  tmpvectorj = 0.0_ReKi; tmpvectorj2 = 0.0_ReKi
  DO IBlade = 1, NumBl
     DO indx = 1, 3
        tmp = RELM( nelm_start:NELM )
        tmp2 = REAL(InputMarkers( IBlade )%Orientation( 2, indx, nelm_start:NELM),ReKi)
        BladeTanVectj( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& REAL(InputMarkers( IBlade )%Orientation( 2, indx, nelm_start:NELM ),ReKi), &
           & Rnumbsp1, NumBS+1, NELM - nelm_start+1 )
        BladeTanVect2j( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& REAL(InputMarkers( IBlade )%Orientation( 2, indx, nelm_start:NELM ),ReKi), &
           & Rnumbs, numbs, NELM - nelm_start+1 )
        tmp2 = REAL(InputMarkers( IBlade )%Orientation( 1, indx, nelm_start:NELM),ReKi)
        BladeNormVect2j( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& REAL(InputMarkers( IBlade )%Orientation( 1, indx, nelm_start:NELM ),ReKi), &
           & Rnumbs, numbs, NELM - nelm_start+1 )
        tmp2 = InputMarkers( IBlade )%Position( indx, nelm_start:NELM)
        BladeQuarterChordj( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !&, (InputMarkers( IBlade )%Position( indx, nelm_start:NELM)),
           & Rnumbsp1, NumBS+1, NELM - nelm_start+1 )
        BladeThreeQuarterChordj( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& (InputMarkers( IBlade )%Position( indx, nelm_start:NELM )),
           & Rnumbs, NumBS, NELM - nelm_start+1 )
        IF ( indx .EQ. 1 ) Then
           BladeQuarterChordj(      indx, :, IBlade ) = 0.0_ReKi
           BladeThreeQuarterChordj( indx, :, IBlade ) = 0.0_ReKi
        END IF

        DO IElement = 1, NumBS
           BladeThreeQuarterChordj( indx, IElement, IBlade ) = BladeThreeQuarterChordj( indx, &
              & IElement, IBlade ) + 0.5_ReKi * C2( IElement ) * BladeTanVect2j( indx, IElement, IBlade )
        END DO ! IElement

        BladeLocj( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& (InputMarkers( IBlade )%Position( indx, nelm_start:NELM)),
           & Rnumbsp1, NumBS+1, NELM - nelm_start+1 )
        BladeLoc2j( indx, :, IBlade ) = interpolation_array( tmp, tmp2, &
           !& (InputMarkers( IBlade )%Position( indx, nelm_start:NELM)),
           & Rnumbs, NumBS, NELM - nelm_start+1 )
     END DO !indx

     DO IElement = 1, NumBS + 1
        CALL TRANSFORM_TO_FVW_COORDS( BladeQuarterChordj( :, IElement, IBlade ))
        BladeQuarterChordj( 1, IElement, IBlade ) = BladeQuarterChordj( 1, IElement, IBlade ) - HH

        CALL TRANSFORM_TO_FVW_COORDS( BladeLocj( :, IElement, IBlade ))
        BladeLocj( 1, IElement, IBlade ) = BladeLocj( 1, IElement, IBlade ) - HH
        BladeLocj( 3, IElement, IBlade ) = 0.0_ReKi
        CALL TRANSFORM_TO_FVW_COORDS( BladeTanVectj( 1:3, IElement, IBlade ))
        IF ( IElement .LE. NumBS ) THEN
           CALL TRANSFORM_TO_FVW_COORDS( BladeNormVect2j( 1:3, IElement, IBlade ))

           CALL TRANSFORM_TO_FVW_COORDS( BladeTanVect2j( 1:3, IElement, IBlade ))

           CALL TRANSFORM_TO_FVW_COORDS( BladeLoc2j( :, IElement, IBlade ))
           BladeLoc2j( 1, IElement, IBlade ) = BladeLoc2j( 1, IElement, IBlade ) - HH
           Bladeloc2j( 3, IElement, IBlade ) = 0.0_ReKi
           CALL TRANSFORM_TO_FVW_COORDS( BladeThreeQuarterChordj( :, IElement, IBlade ))
           BladeThreeQuarterChordj( 1, IElement, IBlade ) = BladeThreeQuarterChordj( 1, &
              & IElement, IBlade ) - HH
        END IF
     END DO ! IElement
  END DO ! IBlade

  DO IBlade = 1, NumBl
     DO IElement = 1, NELM
        PitNow = 0.0_ReKi; SPitch = 0.0_ReKi; CPitch = 0.0_ReKi
        PitNow = - 1.0_ReKi * ATAN2( -1.0_ReKi  * DOT_PRODUCT( TurbineComponents%Blade( IBlade )%&
           & Orientation( 1, : ), InputMarkers( IBlade )%Orientation( 2, :, IElement)), &
           & DOT_PRODUCT( TurbineComponents%Blade( IBlade )%Orientation( 1, : ), &
           & InputMarkers( IBlade )%Orientation( 1, :, IElement)))
        SPitch = SIN( PitNow )
        CPitch = COS( PitNow )
        tmpVectorj = 0.0_ReKi; tmpVectorj2 = 0.0_ReKi
        tmpVectorj = -1.0_ReKi * SPitch * InputMarkers( IBlade )%Orientation( 1, :, IElement) &
           & + CPitch * InputMarkers( IBlade )%Orientation( 2, :, IElement)
        tmpVectorj2 = -InputMarkers( IBlade )%TranslationVel( :, IElement )

        VTT( 1:3, IElement, Iblade ) = DOT_PRODUCT( tmpVectorj, tmpVectorj2) * tmpVectorj

        CALL TRANSFORM_TO_FVW_COORDS( VTT( 1:3, IElement, Iblade ))
        tmpVectorj = 0.0_ReKi; tmpVectorj2 = 0.0_ReKi
        tmpVectorj = CPitch * InputMarkers( IBlade )%Orientation( 1, :, IElement) + &
            & SPitch * InputMarkers( IBlade )%Orientation( 2, :, IElement)
        tmpVectorj2 = InputMarkers( IBlade )%TranslationVel( :, IElement )
        VNElement( 1:3, IElement, Iblade ) = -1.0_ReKi * DOT_PRODUCT( tmpVectorj, tmpVectorj2)* tmpVectorj

       CALL TRANSFORM_TO_FVW_COORDS( VNElement( :, IElement, Iblade ))
     END DO ! IElement

     DO indx = 1, 3
        tmp3 = VTT( indx, :, Iblade )
        VaxialBL( indx, :, IBlade ) = interpolation_array( RELM, &
           & tmp3, Rnumbs, NumBS, NELM )
        tmp3 = VNElement( indx, :, Iblade)
        VNElemBL( indx, :, IBlade ) = interpolation_array( RELM, &
           & tmp3, Rnumbs, NumBS, NELM )
     END DO ! indx
  END DO ! IBlade

  BladeLoc2j_Real = BladeLoc2j

END SUBROUTINE BladeLocations
!==========================================================================

!==========================================================================
SUBROUTINE Predictor( m )

  ! *******************************************************************************************
  ! * Predictor portion of predictor/corrector scheme. Results in predicted wake              *
  ! * position using inital guess values:                                                     *
  ! * r'(ψ+Δψ, ζ+Δζ) = 1/47*{r'(ψ+Δψ,ζ)-3r'(ψ,ζ+Δζ)+45r'(ψ,ζ)+3r'(ψ-Δψ,ζ+Δζ)+3r'(ψ-Δψ,ζ)-     *
  ! *                  r'(ψ-2Δψ,ζ+Δζ)-r'(ψ-2Δψ,ζ)+(48Δψ/4Ω)[4V∞+Vind(r^(n-1)(ψ,ζ))+           *
  ! *                  Vind(r^(n-1)(ψ+Δψ,ζ))+Vind(r^(n-1)(ψ,ζ+Δζ))+Vind(r^(n-1)(ψ+Δψ,ζ+Δζ))]} *
  ! *                                                                                         *
  ! *  where rpcyl1 = r'(ψ,ζ)                                                                 *
  ! *        rpcyl2 = r'(ψ+Δψ,ζ)                                                              *
  ! *        rpcyl3 = r'(ψ-Δψ,ζ)                                                              *
  ! *        rpcyl4 = r'(ψ-Δψ,ζ+Δζ)                                                           *
  ! *        rpcyl5 = r'(ψ,ζ+Δζ)                                                              *
  ! *        rpcyl6 = r'(ψ-2Δψ,ζ)                                                             *
  ! *        rpcyl7 = r'(ψ-2Δψ,ζ+Δζ)                                                          *
  ! *                                                                                         *
  ! *        rocyl1 = GCoord%r_oldjm1(   :, m-1, n )
  ! *        rocyl2 = GCoord%r_oldj(     :, m-1, n )
  ! *        rocyl3 = GCoord%r_oldjm2(   :, m-1, n )
  ! *        rocyl4 = GCoord%r_oldjm2(   :, m,   n )
  ! *        rocyl5 = GCoord%r_oldjm1(   :, m,   n )
  ! *        rocyl6 = GCoord%r_oldjm3(   :, m-1, n )
  ! *        rocyl7 = GCoord%r_oldjm3(   :, m,   n )
  ! *        rocyl8 = GCoord%r_oldj(     :, m,   n )
  ! *                                                                                         *
  ! *        rpcyl8 = r'(ψ+Δψ, ζ+Δζ)                                                          *
  ! *                                                                                         *
  ! *        CalcedVinf    = V∞
  ! *        VinducedTot1  = Vind(r^(n-1)(ψ+Δψ,ζ+Δζ))
  ! *        VinducedTot2  = Vind(r^(n-1)(ψ,ζ))
  ! *        VinducedTot1b = Vind(r^(n-1)(ψ+Δψ,ζ))
  ! *        VinducedTot2b = Vind(r^(n-1)(ψ,ζ+Δζ))
  ! *******************************************************************************************


  IMPLICIT NONE

  INTEGER :: m, a
  REAL :: WhichTurb
  REAL( ReKi ) :: phi, invDeltaPsi

  DO n = 1, NumWakes
     WhichTurb = REAL(n-0.01)/REAL(NumBl)
     q = n - FLOOR( WhichTurb )*NumBl
     IF ( m .GT. CUTOFF_up( CEILING( WhichTurb ))) THEN
        GO TO 100 !7.13.15 see pg. 29 of notebook
     END IF !7.13.15 see pg. 29 of notebook

     CALL VinducedNW( NWake%r_nearj, NWake%Gamma_nearj, FWake%r_oldj( :, m, n ), VinducedNW1, &
        & NWake%r_nearjm1, NumWakes )

     CALL VinducedBC( BladeQuarterChordj, NWake%Gammablj, FWake%r_oldj( :, m, n ), VinducedBC1 )

     CALL Vinduced2OLD( FWake%r_oldj, FWake%Gammaj, FWake%r_oldj( :, m, n ), &
        & FWake%r_oldjm1, n, j, m )
                         !VinducedFarWakeOld values come from Vinduced2OLD
     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakej( 1, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakej( 2, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakej( 3, m, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot1 = VinducedFW1 + VinducedNW1 + VinducedBC1

     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearjm1, NWake%Gamma_nearjm1, FWake%r_oldjm1( :, m-1, n ), VinducedNW1, &
        & NWake%r_nearjm2, NumWakes )

     CALL VinducedBC( BladeQuarterChordjm1, NWake%Gammabljm1, FWake%r_oldjm1( :, m-1, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakejm1( 1, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakejm1( 2, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakejm1( 3, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot2 = VinducedFW1 + VinducedNW1 + VinducedBC1


     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearj, NWake%Gamma_nearj, FWake%r_oldj( :, m-1, n ), VinducedNW1, &
        & NWake%r_nearjm1, NumWakes )

     CALL VinducedBC( BladeQuarterChordj, NWake%Gammablj, FWake%r_oldj( :, m-1, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakej( 1, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakej( 2, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakej( 3, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedTot1b = VinducedFW1 + VinducedNW1 + VinducedBC1


     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearjm1, NWake%Gamma_nearjm1, FWake%r_oldjm1( :, m, n ), VinducedNW1, &
        & NWake%r_nearjm2, NumWakes )

     CALL VinducedBC( BladeQuarterChordjm1, NWake%Gammabljm1, FWake%r_oldjm1( :, m, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakejm1( 1, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakejm1( 2, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakejm1( 3, m, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot2b = VinducedFW1 + VinducedNW1 + VinducedBC1

     tmpvector = FWake%r_oldj( :, m, n )
     CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )

     Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector
     CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )
     CalcedVinf = Wind_FVW%OutputData%VelocityUVW( :, 1 )

     CALL TRANSFORM_TO_FVW_COORDS( CalcedVinf )

     rpcyl1 = FWake%r_primejm1( :, m-1, n )
     rpcyl2 = FWake%r_primej(   :, m-1, n )
     rpcyl3 = FWake%r_primejm2( :, m-1, n )
     rpcyl4 = FWake%r_primejm2( :, m,   n )
     rpcyl5 = FWake%r_primejm1( :, m,   n )
     rpcyl6 = FWake%r_primejm3( :, m-1, n )
     rpcyl7 = FWake%r_primejm3( :, m,   n )

           !predictor

!     rpcyl8_cop = 1.d0 / 47.d0 * ( 45.d0 * rpcyl1 + rpcyl2 + 3.d0 * rpcyl3 + 3.d0 * rpcyl4 - 3.d0 * &
!        & rpcyl5 - rpcyl6 - rpcyl7 ) + 48.d0 / 47.d0 * delta_psi / Omega * ( CalcedVinf + 0.25d0 * &
!        & VinducedTot1 + 0.25d0 * VinducedTot2 + 0.25d0 * VinducedTot1b + 0.25d0 * VinducedTot2b )

  !!!!!!!!!!!NEW PREDICTOR SCHEME!!!  Added 5.8.17
     phi = 46.0_ReKi*delta_psi(1) + 4.0_ReKi*delta_psi(2) - 2.0_ReKi*delta_psi(3)
     invDeltaPsi = 1.0_ReKi/(2.0_ReKi*delta_psi(1))

     rpcyl8 = (( CalcedVinf + 0.25_ReKi * VinducedTot1 + 0.25_ReKi * VinducedTot2 + 0.25_ReKi * &
            & VinducedTot1b + 0.25_ReKi * VinducedTot2b ) / Omega - ( -1.0_ReKi*invDeltaPsi + &
            & 23.0_ReKi / phi ) * rpcyl2 - ( invDeltaPsi - 21.0_ReKi / phi ) * rpcyl5 + &
            & ( invDeltaPsi + 21.0_ReKi / phi ) * rpcyl1 + ( 3.0_ReKi / phi ) * rpcyl4 + &
            & ( 3.0_ReKi / phi ) * rpcyl3 - ( 1.0_ReKi / phi ) * rpcyl7 - &
            & ( 1.0_ReKi / phi ) * rpcyl6 ) / ( invDeltaPsi + 23.0_ReKi / phi )


     FWake%r_primej( :, m, n ) = rpcyl8
100 a = a !7.13.15 see pg. 29 of notebook
  END DO ! n



END SUBROUTINE Predictor
!==========================================================================

!==========================================================================
SUBROUTINE Corrector( m )

  USE MultTurb_Params, Only: NTurb
  USE FileManipulation


  IMPLICIT NONE

  INTEGER :: m, a
  REAL :: WhichTurb
  REAL( ReKi ) :: phi, invDeltaPsi

  DO n = 1, NumWakes
     WhichTurb = REAL(n-0.01)/REAL(NumBl)
     q = n - FLOOR( WhichTurb )*NumBl
     IF ( m .GT. CUTOFF_up( CEILING( WhichTurb ))) THEN
        GO TO 110 !7.13.15 see pg. 29 of notebook
     END IF !7.13.15 see pg. 29 of notebook
     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL Vinduced2PRIME( FWake%r_primej, FWake%Gammaj, FWake%r_primej( :, m, n ), FWake%r_primejm1, &
        & n, j, m )
     CALL VinducedNW( NWake%r_nearj, NWake%Gamma_nearj, FWake%r_primej( :, m, n ), VinducedNW1, &
        & NWake%r_nearjm1, NumWakes )
     CALL VinducedBC( BladeQuarterChordj, NWake%Gammablj, FWake%r_primej( :, m, n ), VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakej( 1, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakej( 2, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakej( 3, m, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot1 = VinducedFW1 + VinducedNW1 + VinducedBC1


     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearjm1, NWake%Gamma_nearjm1, FWake%r_primejm1( :, m-1, n ), VinducedNW1, &
        & NWake%r_nearjm2, NumWakes )

     CALL VinducedBC( BladeQuarterChordjm1, NWake%Gammabljm1, FWake%r_primejm1( :, m-1, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakejm1( 1, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakejm1( 2, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakejm1( 3, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot2 = VinducedFW1 + VinducedNW1 + VinducedBC1


     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearj, NWake%Gamma_nearj, FWake%r_primej( :, m-1, n ), VinducedNW1, &
        & NWake%r_nearjm1, NumWakes )

     CALL VinducedBC( BladeQuarterChordj, NWake%Gammablj, FWake%r_primej( :, m-1, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakej( 1, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakej( 2, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakej( 3, m-1, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot1b = VinducedFW1 + VinducedNW1 + VinducedBC1


     VinducedFW1 = 0.0_ReKi; VinducedNW1 = 0.0_ReKi; VinducedBC1 = 0.0_ReKi
     CALL VinducedNW( NWake%r_nearjm1, NWake%Gamma_nearj, FWake%r_primejm1( :, m, n ), VinducedNW1, &
        & NWake%r_nearjm2, NumWakes )

     CALL VinducedBC( BladeQuarterChordjm1, NWake%Gammabljm1, FWake%r_primejm1( :, m, n ), &
        & VinducedBC1 )

     VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakejm1( 1, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakejm1( 2, m, q, 2:WakeAgeLimit, 1:NumWakes ))
     VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakejm1( 3, m, q, 2:WakeAgeLimit, 1:NumWakes ))

     VinducedTot2b = VinducedFW1 + VinducedNW1 + VinducedBC1


     tmpvector = FWake%r_primej( :, m, n )

     CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )

     Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector

     CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )
     CalcedVinf = Wind_FVW%OutputData%VelocityUVW( :, 1 )

     CALL TRANSFORM_TO_FVW_COORDS( CalcedVinf )

     rncyl1 = FWake%r_newjm1( :, m-1, n )
     rncyl2 = FWake%r_newj(   :, m-1, n )
     rncyl3 = FWake%r_newjm2( :, m-1, n )
     rncyl4 = FWake%r_newjm2( :, m,   n )
     rncyl5 = FWake%r_newjm1( :, m,   n )
     rncyl6 = FWake%r_newjm3( :, m-1, n )
     rncyl7 = FWake%r_newjm3( :, m,   n )

           !corrector

   !  rncyl8_cop = 1.d0 / 47.d0 * ( 45.d0 * rncyl1 + rncyl2 + 3.d0 * rncyl3 + 3.d0 * rncyl4 - 3.d0 * &
   !     & rncyl5 - rncyl6 - rncyl7 ) + 48.d0 / 47.d0 * delta_psi / Omega * ( CalcedVinf + 0.25d0 * &
   !     & VinducedTot3 + 0.25d0 * VinducedTot4 + 0.25d0 * VinducedTot3b + 0.25d0 * VinducedTot4b )

  !!!!!!!!!!!NEW CORRECTOR SCHEME!!!  Added 5.8.17
     phi = 46.0_ReKi * delta_psi( 1 ) + 4.0_ReKi * delta_psi( 2 ) - 2.0_ReKi * delta_psi( 3 )
     invDeltaPsi = 1.0_ReKi / ( 2.0_ReKi * delta_psi( 1 ))

     rncyl8 = (( CalcedVinf + 0.25_ReKi * VinducedTot1 + 0.25_ReKi * VinducedTot2 + &
              & 0.25_ReKi * VinducedTot1b + 0.25_ReKi * VinducedTot2b ) / Omega - ( -1.0_ReKi*invDeltaPsi + &
            & 23.0_ReKi / phi ) * rncyl2 - ( invDeltaPsi - 21.0_ReKi / phi ) * rncyl5 + &
            & ( invDeltaPsi + 21.0_ReKi / phi ) * rncyl1 + ( 3.0_ReKi / phi ) * rncyl4 + &
            & ( 3.0_ReKi / phi ) * rncyl3 - ( 1.0_ReKi / phi ) * rncyl7 - ( 1.0_ReKi / phi ) * rncyl6 ) / &
            & ( invDeltaPsi + 23.0_ReKi / phi )

     FWake%r_newj( :, m, n ) =  rncyl8

110 a = a !7.13.15 see pg. 29 of notebook

  END DO ! n

END SUBROUTINE Corrector
!==========================================================================

!==========================================================================
SUBROUTINE UpdateAeroVals

  USE FileManipulation, Only: WakeVelProfile

  IMPLICIT NONE

  INTEGER :: Size_Far, Size_NearV, Size_FarG
  REAL :: WhichTurb
  LOGICAL :: file_exists

  Size_NearV   = NumBl * WakeAgeLimit * NumBl * NumBS * 3  ! for velocity of near wake markers
  Size_Far  = NumBl * CUTOFF_Allocate * 3  ! for location of far wake markers
!For # of wakes ^   For Far wake^^     ^For wake position
  Size_FarG = NumBl * CUTOFF_Allocate  ! for circulation of far wake markers
  !CALL TransformToGlobal
  DO n = 1, NumBl
!Need to collect appropriate rj, Gammaj, and rjm1 data here!!
     !CALL TransformToGlobal
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
     !IF( NumTurbs .GT. 1 ) THEN
     !   DO ProcNum = 0, NumTurbs-1
     !      StartCount = ProcNum*NumBl+1
     !      CALL MPI_BCAST( GCoord%rj( 1, 1, StartCount ), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !      CALL MPI_BCAST( GCoord%rjm1( 1, 1, StartCount ), Size_Far, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !      CALL MPI_BCAST( GCoord%Gammaj( 1, StartCount ), Size_FarG, MPI_REAL8, ProcNum, MPI_COMM_WORLD, ierr  )
     !   END DO
     !END IF
     !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

     VTotal = 0.0_ReKi
     VindTotal = 0.0_ReKi
     k = 1

     WhichTurb = REAL(n-0.01)/REAL(NumBl)
     q = n + FLOOR( WhichTurb )*NumBl
     DO nbs = Num_start, NumBS
        CALL Vinduced3( FWake%rj, FWake%Gammaj, BladeThreeQuarterChordj( :, nbs, n  ), &
           & FWake%rjm1, q, j, nbs )

        !adds the induced velocity at the blade from the 3 Time steps of near wake
        VindTotal( 1, nbs ) = Sum( FWake%VinducedFarWakeRj( 1, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
        VindTotal( 2, nbs ) = Sum( FWake%VinducedFarWakeRj( 2, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))
        VindTotal( 3, nbs ) = Sum( FWake%VinducedFarWakeRj( 3, nbs, q, 2:WakeAgeLimit, 1:NumWakes ))

        tmpvector = BladeLoc2j( :, nbs, n )

        !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !CALL TransformToLocal( )
        VinducedNW1 = 0.0_ReKi
        CALL VinducedNW( NWake%r_nearj, NWake%Gamma_nearj, BladeThreeQuarterChordj( :, nbs, n ), &
           & VinducedNW1, NWake%r_nearjm1, NumBl )

        CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )

        Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector
        CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )
        CalcedVinf = Wind_FVW%OutputData%VelocityUVW( :, 1 )

        CALL TRANSFORM_TO_FVW_COORDS( CalcedVinf )

        VinducedNWFinal( :, nbs, n ) = VinducedNW1

        VTotal( :, nbs ) = CalcedVinf + VindTotal( :, nbs ) + Vaxial2j( :, nbs, n ) + &
           & VinducedNWFinal( :, nbs, n ) + VNElem2j( :, nbs, n )

        VN  = DOT_PRODUCT( VTotal( :, nbs ), BladeNormVect2j( :, nbs, n ) )
        VT  = DOT_PRODUCT( VTotal( :, nbs ), BladeTanVect2j(  :, nbs, n ) )

        a_of_a_storej( nbs, n ) = atan( VN / VT ) !****NEED TO CHANGE ****
        Velstorej2( nbs, n ) = VT * VT + VN * VN
        Velstorej( nbs, n ) = sqrt( Velstorej2( nbs, n ))  !KS -- I don't think this has to be an array   !10.14.15


        cl_storej( nbs, n ) = -2.0_ReKi * NWake%Gammablj( nbs, n ) / ( Velstorej( nbs, n ) * C2( nbs ))
        Vind_storej( :, nbs, n ) = VinducedNWFinal( :, nbs, n ) + VindTotal( :, nbs )
     END DO ! nbs (NumBS)
  ENDDO ! n (NumBl)

  FWake%VinducedFarWakejm1 = FWake%VinducedFarWakej
  FWake%VinducedFarWakeRjm1 = FWake%VinducedFarWakeRj
  FWake%VinducedFarWakej = 0.0_ReKi
  FWake%VinducedFarWakeRj = 0.0_ReKi

  FWake%r_oldjm3 = FWake%r_oldjm2; FWake%r_oldjm2 = FWake%r_oldjm1;
  FWake%r_oldjm1 = FWake%r_oldj
  FWake%r_oldj = FWake%r_newj; FWake%r_primej = FWake%r_newj
  FWake%r_newjm3 = FWake%r_newjm2; FWake%r_newjm2 = FWake%r_newjm1
  FWake%r_newjm1 = FWake%r_newj;
  FWake%r_newj = 0.0_ReKi
  FWake%r_primejm3 = FWake%r_primejm2
  FWake%r_primejm2 = FWake%r_primejm1; FWake%r_primejm1 = FWake%r_primej

  FWake%Gammajm1=FWake%Gammaj; FWake%Gammaj=FWake%Gammajp1; NWake%Gammabljm1=NWake%Gammablj
  NWake%Gammablj=0.0_ReKi; NWake%Gamma_nearjm1=NWake%Gamma_nearj; NWake%Gamma_nearj=NWake%Gamma_nearjp1
  NWake%Gamma_nearjp1=0.0_ReKi

  NWake%r_nearjm2 = NWake%r_nearjm1; NWake%r_nearjm1 = NWake%r_nearj; FWake%rjm2 = FWake%rjm1
  FWake%rjm1 = FWake%rj; FWake%rj = 0.0_ReKi

  BladeQuarterChordjm2 = BladeQuarterChordjm1
  BladeQuarterChordjm1 = BladeQuarterChordj
  BladeQuarterChordj = 0.0_ReKi

  !CALL TransformToGlobal( )
  !CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  DO n = 1, NumBl
     AofA(  :, n ) = interpolation_array( Rnumbs, a_of_a_storej( :, n ), (RELM), NELM, NumBS )
     W2FVW( :, n ) = interpolation_array( Rnumbs, Velstorej2(    :, n ), (RELM), NELM, NumBS )
     CLFVW( :, n ) = interpolation_array( Rnumbs, cl_storej(     :, n ), (RELM), NELM, NumBS )

     VINDFVW( 1, :, n ) = interpolation_array( Rnumbs,          Vind_storej( 3, :, n ), (RELM), &
        & NELM, NumBS )
     VINDFVW( 2, :, n ) = interpolation_array( Rnumbs, -1.0_ReKi * Vind_storej( 2, :, n ), (RELM), &
        & NELM, NumBS )
     VINDFVW( 3, :, n ) = interpolation_array( Rnumbs,          Vind_storej( 1, :, n ), (RELM), &
        & NELM, NumBS )
  END DO

  IF (Time_Real .GE. TMax-(TMax/10.0_ReKi)) THEN
     INQUIRE(FILE="InputFiles/WakePoints.txt", EXIST=file_exists)
     IF (file_exists .EQV. .TRUE.) THEN
        CALL WakeVelProfile( zloc, Wind_FVW, j, FWake%rjm1, FWake%Gammajm1, FWake%rjm2, &
                           & NWake%r_nearjm1, NWake%Gamma_nearjm1, NWake%r_nearjm2, &
                           & BladeQuarterChordjm1, NWake%Gammabljm1, BladeQuarterChordjm2 )
     END IF
  END IF


  !FWake%r_oldjm3 = 0.0_ReKi; FWake%r_oldjm2 = 0.0_ReKi;
  !FWake%r_oldjm1 = 0.0_ReKi
  !FWake%r_oldj = 0.0_ReKi; FWake%r_primej = 0.0_ReKi; FWake%r_newjm3 = 0.0_ReKi
  !FWake%r_newjm2 = 0.0_ReKi; FWake%r_newjm1 = 0.0_ReKi; FWake%r_primejm3 = 0.0_ReKi
  !FWake%r_primejm2 = 0.0_ReKi; FWake%r_primejm1 = 0.0_ReKi

  !FWake%Gammajm1= 0.0_ReKi; NWake%Gammabljm1 = 0.0_ReKi; NWake%Gammablj = 0.0_ReKi
  !NWake%Gamma_nearjm1 = 0.0_ReKi; NWake%Gamma_nearj = 0.0_ReKi; NWake%Gamma_nearjp1 = 0.0_ReKi

  !NWake%r_nearjm2 = 0.0_ReKi; NWake%r_nearjm1 = 0.0_ReKi; FWake%rjm2 = 0.0_ReKi
  !FWake%rj = 0.0_ReKi

  !BladeQuarterChordjm2 = 0.0_ReKi; BladeQuarterChordjm1 = 0.0_ReKi
  !BladeQuarterChordj = 0.0_ReKi

END SUBROUTINE UpdateAeroVals
!==========================================================================

!==========================================================================
SUBROUTINE WriteWake( )

  IMPLICIT NONE

  INTEGER :: kindx
  DO kindx = 1, CUTOFF(1)
        WRITE( 101+NTurb, 200 ) FWake%rjm1( :, kindx, 1 ), &
                       & FWake%Gammaj( kindx, 1), Sum( FWake%VinducedFarWakejm1( 1, kindx, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 2, kindx, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 3, kindx, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes ))
        WRITE( 201+NTurb, 200 ) FWake%rjm1( :, kindx, 2 ), &
                       & FWake%Gammaj( kindx, 2), Sum( FWake%VinducedFarWakejm1( 1, kindx, 2, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 2, kindx, 2, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 3, kindx, 2, NnearMax+1:WakeAgeLimit, 1:NumWakes ))
        IF ( NumBl .EQ. 3 ) THEN
           WRITE( 301+NTurb, 200 ) FWake%rjm1( :, kindx, 3 ), &
                       & FWake%Gammaj( kindx, 3), Sum( FWake%VinducedFarWakejm1( 1, kindx, 3, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 2, kindx, 3, NnearMax+1:WakeAgeLimit, 1:NumWakes )), Sum( FWake%VinducedFarWakejm1( 3, kindx, 3, NnearMax+1:WakeAgeLimit, 1:NumWakes ))
        END IF
  END DO

  FWake%rjm1 = 0.0_ReKi
  FWake%Gammaj =0.0_ReKi


200  FORMAT( 5F15.7 )
300  FORMAT( 12F15.7 )

END SUBROUTINE WriteWake
!==========================================================================

!==========================================================================
  FUNCTION interpolation_array( xvals, yvals, xi, arrayilen, arrayvlen )

    USE precision

    IMPLICIT NONE

    INTEGER arrayilen, arrayvlen, arindx, ilo
    REAL(ReKi), DIMENSION( arrayilen ) :: interpolation_array, xi
    REAL(ReKi), DIMENSION( arrayvlen ) :: xvals, yvals, tmp2, tmp3
    REAL(ReKi)                         :: tmp1
    ilo = 1

    DO arindx = 1, arrayilen
       IF ( xi( arindx ) .LT. xvals( 1 )) THEN
          interpolation_array( arindx ) = yvals( 1 ) + ( xi( arindx ) - xvals( 1 )) / &
              & ( xvals( 2 ) - xvals( 1 )) * ( yvals( 2 ) - yvals( 1 ))
       ELSE IF ( xi( arindx ) .GT. xvals( arrayvlen )) THEN
          interpolation_array( arindx ) = yvals( arrayvlen - 1 ) + ( xi( arindx ) - &
             & xvals( arrayvlen - 1 )) / ( xvals( arrayvlen ) - xvals( arrayvlen - 1 )) * &
             & ( yvals( arrayvlen ) - yvals( arrayvlen - 1 ))
       ELSE
          tmp1 = real( xi( arindx ), ReKi)
          tmp2 = real( xvals , ReKi)
          tmp3 = real( yvals , ReKi)
          !interpolation_array( arindx ) = ( InterpBinReal( real( xi( arindx )), &
             !& real( xvals ), real( yvals ), ilo, arrayvlen ))

          interpolation_array( arindx ) = ( InterpBinReal( tmp1, tmp2, tmp3, ilo, arrayvlen ))
       END IF
    END DO

  END FUNCTION interpolation_array
!==========================================================================

END SUBROUTINE FVW_COMPUTE_WAKE
