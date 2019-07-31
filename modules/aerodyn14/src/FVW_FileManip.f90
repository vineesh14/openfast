MODULE FileManipulation

  USE MultTurb_Params

  IMPLICIT NONE

  CONTAINS

  !************************************
  ! MultTurb subroutines for opening, writing to,
  ! and creating files
  !************************************
  ! Kelsey Shaler 8/28/14

!=================================================

!=================================================
SUBROUTINE WriteInitWake( CUTOFF_init  )

  USE NWTC_Library
  USE FVW_Subs
  USE FVW_Parm, Only: delta_psi_Est, CUTOFF_upinit, NnearMax, NumBS, NumBl, Radius, CUTOFF_upmax, &
                     & Num_Start, Root_cut, RElm, RotSpeed_Est
  USE MultTurb_Params, Only: Turbines, PerUinf

  IMPLICIT NONE

  INTEGER :: I, J, nindx, kindx, jindx, nindx2, kindx2, kindx3, nbs, j2, CUTOFF_Write, WakeAgeLimit
  INTEGER, DIMENSION( : ) :: CUTOFF_init
  REAL( ReKi ), ALLOCATABLE, DIMENSION( : ) :: BladeLocy, BladeLocx
  REAL( ReKi ), ALLOCATABLE, DIMENSION( :, : ) :: gammaj, gamma_bl
  REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, : ) :: gamma_near
  REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, : ) :: far_wake, BladeQuarterChord1, BladeQuarterChord2
  REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, : ) :: near_wake
  REAL( ReKi ) :: const, Theta, Theta_init, BladeSec, Hub, LocRad, Alpha, Beta
  REAL( ReKi ) :: BladeLoc, Uinf, GammaInit
  CHARACTER( * ), PARAMETER :: root = 'InitFiles/'

  NumTurbs = 1

  PRINT*, ''
  PRINT*, '*******************************'
  PRINT*, 'Writing new InitFiles'

  PRINT*, ''

    WakeAgeLimit = NINT( TwoPi_D / delta_psi_Est ) + CUTOFF_upmax(1)

  OPEN( unit = 3000, file = 'InitialWake.txt' )
  READ( 3000, * ) Hub
  READ( 3000, * ) Uinf
  READ( 3000, * ) GammaInit


  DO J = 1, NumTurbs
  OPEN( unit = 1000, file = ( root//TRIM( 'Turb1' )//'_r_primej.txt'   ), STATUS = 'new')
  OPEN( unit = 1001, file = ( root//TRIM( 'Turb1' )//'_r_primejm1.txt' ), STATUS = 'new')
  OPEN( unit = 1002, file = ( root//TRIM( 'Turb1' )//'_r_primejm2.txt' ), STATUS = 'new')
  OPEN( unit = 1003, file = ( root//TRIM( 'Turb1' )//'_r_primejm3.txt' ), STATUS = 'new')

  OPEN( unit = 1025, file = ( root//TRIM( 'Turb1' )//'_r_oldj.txt'   ), STATUS = 'new')
  OPEN( unit = 1004, file = ( root//TRIM( 'Turb1' )//'_r_oldjm1.txt' ), STATUS = 'new')
  OPEN( unit = 1005, file = ( root//TRIM( 'Turb1' )//'_r_oldjm2.txt' ), STATUS = 'new')
  OPEN( unit = 1006, file = ( root//TRIM( 'Turb1' )//'_r_oldjm3.txt' ), STATUS = 'new')

  OPEN( unit = 1007, file = ( root//TRIM( 'Turb1' )//'_r_newjm1.txt' ), STATUS = 'new')
  OPEN( unit = 1008, file = ( root//TRIM( 'Turb1' )//'_r_newjm2.txt' ), STATUS = 'new')
  OPEN( unit = 1009, file = ( root//TRIM( 'Turb1' )//'_r_newjm3.txt' ), STATUS = 'new')

  OPEN( unit = 1010, file = ( root//TRIM( 'Turb1' )//'_rjm1.txt' ), STATUS = 'new')
  OPEN( unit = 1011, file = ( root//TRIM( 'Turb1' )//'_rjm2.txt' ), STATUS = 'new')

  OPEN( unit = 1012, file = ( root//TRIM( 'Turb1' )//'_r_nearjm1.txt' ), STATUS = 'new')
  OPEN( unit = 1013, file = ( root//TRIM( 'Turb1' )//'_r_nearjm2.txt' ), STATUS = 'new')

  OPEN( unit = 1017, file = ( root//TRIM( 'Turb1' )//'_Gammablj.txt'   ), STATUS = 'new')
  OPEN( unit = 1018, file = ( root//TRIM( 'Turb1' )//'_Gammabljm1.txt' ), STATUS = 'new')

  OPEN( unit = 1019, file = ( root//TRIM( 'Turb1' )//'_Gammaj.txt'   ), STATUS = 'new')
  OPEN( unit = 1020, file = ( root//TRIM( 'Turb1' )//'_Gammajm1.txt' ), STATUS = 'new')

  OPEN( unit = 1021, file = ( root//TRIM( 'Turb1' )//'_Gamma_nearj.txt'   ), STATUS = 'new')
  OPEN( unit = 1022, file = ( root//TRIM( 'Turb1' )//'_Gamma_nearjm1.txt' ), STATUS = 'new')

  OPEN( unit = 1023, file = ( root//TRIM( 'Turb1' )//'_BladeQuarterChordjm1.txt' ), STATUS = 'new')
  OPEN( unit = 1024, file = ( root//TRIM( 'Turb1' )//'_BladeQuarterChordjm2.txt' ), STATUS = 'new')

!  Allocate & initialize values

  CUTOFF_Write = CUTOFF_upinit( J )
  ALLOCATE( far_wake( 3, CUTOFF_Write, NumBl ), near_wake( 3, NumBS+1, NnearMax, NumBl ), gammaj( CUTOFF_Write, NumBl ), gamma_bl( NumBS, NumBl ), gamma_near( NnearMax, NumBS+1, NumBl ), BladeQuarterChord1( 3, NumBS+1, NumBl ), BladeQuarterChord2( 3, NumBS+1, NumBl ), BladeLocx( NumBS ), BladeLocy( NumBS ))

  gammaj = GammaInit; gamma_near = 0.50_ReKi

! Set up parameters to be used
  const = Uinf*PerUinf/RotSpeed_Est

  BladeSec = Radius/NumBS
  LocRad = Radius+Hub !To get coordinates right for far wake
  Num_start = INT( Root_cut * NumBS + 1 )

! Initial Blade circulation
  gamma_bl = -5.00_ReKi; gamma_bl( 1:Num_start, : ) = 0.0_ReKi

! Initial near wake circulation
  gamma_near( 1:NnearMax, 14:NumBS+1, : ) = -0.50_ReKi

  BladeQuarterChord1 = 0.00_ReKi; BladeQuarterChord2 = 0.00_ReKi
  IF ( NumBl .EQ. 2 ) THEN
        BladeQuarterChord1( 1, 1:NumBS, 1 ) = RElm( : ) + Hub
        BladeQuarterChord1( 1, NumBS+1, 1 ) = LocRad

        BladeQuarterChord1( :, :, 2 ) = -BladeQuarterChord1( :, :, 1 )
  ELSE
     Theta = -10.00_ReKi*TwoPi_D/360.0_ReKi
     Alpha = 120.00_ReKi*TwoPi_D/360.0_ReKi
     Beta = 30.00_ReKi*TwoPi_D/360.0_ReKi !angle from x-axis to blade)
        BladeQuarterChord1( 1, 1:NumBS, 1 ) = RElm( : )*cos(Theta)
        BladeQuarterChord1( 1, NumBS+1, 1 ) = LocRad*cos(Theta)
        BladeQuarterChord1( 2, 1:NumBS, 1 ) = RElm( : )*sin(Theta)
        BladeQuarterChord1( 2, NumBS+1, 1 ) = LocRad*sin(Theta)

     ! Blade 2 (down and back)
    BladeLocy = RElm(:)*Sin( Theta+Alpha ); BladeLocx = RElm(:)*Cos( Theta+Alpha )
        BladeQuarterChord1( 1, 1:NumBS , 2 ) = BladeLocx(:)
        BladeQuarterChord1( 2, 1:NumBS , 2 ) = BladeLocy(:)
        BladeQuarterChord1( 1, NumBS+1, 2 ) = LocRad*Cos( Theta+Alpha )
        BladeQuarterChord1( 2, NumBS+1, 2 ) = LocRad*Sin( Theta+Alpha )

     ! Blade 3 (down and foward)
    BladeLocy = RElm(:)*Sin( Theta-Alpha ); BladeLocx = RElm(:)*Cos( Theta-Alpha )
        BladeQuarterChord1( 1, 1:NumBS , 3 ) = BladeLocx(:)
        BladeQuarterChord1( 2, 1:NumBS , 3 ) = BladeLocy(:)
        BladeQuarterChord1( 1, NumBS+1, 3 ) = LocRad*Cos( Theta-Alpha )
        BladeQuarterChord1( 2, NumBS+1, 3 ) = LocRad*Sin( Theta-Alpha )

  END IF

  BladeQuarterChord2 = BladeQuarterChord1

! Create helix for far wake
  far_wake = 0.00_ReKi
  IF (NumBl .EQ. 2 ) THEN
     Theta=0.00_ReKi
     DO I = 1, CUTOFF_Write
        far_wake( 1, I, 1 ) = LocRad * cos( -Theta )
        far_wake( 2, I, 1 ) = LocRad * sin( -Theta )
        far_wake( 3, I, 1 ) = const*Theta

        far_wake( 1, I, 2 ) = -LocRad*cos( -Theta )
        far_wake( 2 ,I, 2 ) = -LocRad*sin( -Theta )
        far_wake( 3, I, 2 ) = const*Theta

        Theta = Theta+delta_psi_Est
     END DO
  ELSE
     Theta_init=-10.00_ReKi*TwoPi_D/360.0_ReKi
     Theta=-10.00_ReKi*TwoPi_D/360.0_ReKi
     DO I = 1, CUTOFF_Write
     ! 1st blade (vertical)
        far_wake( 1, I, 1 ) = LocRad * cos( Theta )
        far_wake( 2, I, 1 ) = LocRad * sin( Theta )
        far_wake( 3, I, 1 ) = -const * (Theta+Theta_init)
     ! 2nd blade (down and back)
        far_wake( 1, I, 2 ) = LocRad*cos( Theta + Alpha )
        far_wake( 2 ,I, 2 ) = LocRad*sin( Theta + Alpha )
        far_wake( 3, I, 2 ) = -const * (Theta+Theta_init)
     ! 3rd blade (down and forward)
        far_wake( 1, I, 3 ) = LocRad*cos( Theta - Alpha )
        far_wake( 2 ,I, 3 ) = LocRad*sin( Theta - Alpha )
        far_wake( 3, I, 3 ) = -const * (Theta+Theta_init)

        Theta = Theta-delta_psi_Est
     END DO
  END IF
!Create initial near wake region
  near_wake = 0.00_ReKi

  IF (NumBl .EQ. 2 ) THEN
  near_wake( :, Num_start:NumBS+1, 1, : ) = BladeQuarterChord1( :, Num_start:NumBS+1, : )
     Theta = 0.00_ReKi
     DO I = 2, NnearMax
        Theta = Theta + delta_psi_Est
        near_wake( 1, Num_start:NumBS+1, I, 1 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*cos( -Theta )
        near_wake( 2, Num_start:NumBS+1, I, 1 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*sin( -Theta )
     END DO
     near_wake( 1, :, :, 2 ) = -near_wake( 1, :, :, 1 )
     near_wake( 2, :, :, 2 ) = -near_wake( 2, :, :, 1 )

     near_wake( 3, NumBS+1, :, : ) = far_wake(3, 1:NnearMax, : )
     DO I = Num_start, NumBS
        near_wake( 3, I, :, : ) = near_wake( 3, NumBS+1, :, : )
     END DO

     near_wake( :, 1:Num_start-1, :, : ) = 0.00_ReKi
  ELSE
     Theta = -10.00_ReKi*TwoPi_D/360.0_ReKi
     near_wake( :, Num_start:NumBS+1, 1, : ) = BladeQuarterChord1( :, Num_start:NumBS+1, : )
     DO I = 2, NnearMax
   !Converting polar to cartesian coordinates
        Theta = Theta-delta_psi_Est
        near_wake( 1, Num_start:NumBS+1, I, 1 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*cos( Theta )
        near_wake( 2, Num_start:NumBS+1, I, 1 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*sin( Theta )

        near_wake( 1, Num_start:NumBS+1, I, 2 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*cos( Theta+Alpha )
        near_wake( 2, Num_start:NumBS+1, I, 2 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*sin( Theta+Alpha )

        near_wake( 1, Num_start:NumBS+1, I, 3 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*cos( Theta-Alpha )
        near_wake( 2, Num_start:NumBS+1, I, 3 ) = BladeQuarterChord1( 1, Num_start:NumBS+1, 1 )*sin( Theta-Alpha )

     END DO
      near_wake( 3, NumBS+1, :, : ) = far_wake(3, 1:NnearMax, : )
     DO I = 1, NumBS
        near_wake( 3, I, :, : ) = near_wake( 3, NumBS+1, :, : )
     END DO

     near_wake( :, 1:Num_start-1, :, : ) = 0.00_ReKi
  END IF

! Write files for initial wake
  WRITE ( 1000, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1001, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1002, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1003, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1025, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1004, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1005, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1006, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1007, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1008, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1009, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1010, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )
  WRITE ( 1011, * ) ((( far_wake( j2, kindx, nindx ), j2=1,3),kindx=1, CUTOFF_write), nindx=1,NumBl )

  WRITE ( 1012, * ) (((( near_wake( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )
  WRITE ( 1013, * ) (((( near_wake( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )

  WRITE ( 1017, * ) ((gamma_bl( nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )
  WRITE ( 1018, * ) ((gamma_bl( nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )

  WRITE ( 1019, * ) ((gammaj( kindx, nindx ), kindx = 1, CUTOFF_write), nindx = 1, NumBl )
  WRITE ( 1020, * ) ((gammaj( kindx, nindx ), kindx = 1, CUTOFF_write), nindx = 1, NumBl )

  WRITE ( 1021, * ) (((gamma_near( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )
  WRITE ( 1022, * ) (((gamma_near( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )

  WRITE ( 1023, * ) ((( BladeQuarterChord1( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )
  WRITE ( 1024, * ) ((( BladeQuarterChord2( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )

  CLOSE( 1001 ); CLOSE( 1002 ); CLOSE( 1003 ); CLOSE( 1004 ); CLOSE( 1005 ); CLOSE( 1006 )
  CLOSE( 1007 ); CLOSE( 1008 ); CLOSE( 1009 ); CLOSE( 1010 ); CLOSE( 1011 ); CLOSE( 1012 )
  CLOSE( 1013 )
  CLOSE( 1017 ); CLOSE( 1018 )
  CLOSE( 1019 ); CLOSE( 1020 ); CLOSE( 1021 ); CLOSE( 1022 ); CLOSE( 1023 ); CLOSE( 1024 )
  CLOSE( 1000 ); CLOSE( 1025 )

  DEALLOCATE( far_wake, near_wake, gammaj, gamma_bl, gamma_near, BladeQuarterChord1, BladeQuarterChord2, BladeLocy, BladeLocx )

  END DO

  PRINT*, 'Finished new InitFiles'
  PRINT*, '*******************************'
  PRINT*, ''

END SUBROUTINE WriteInitWake
!=================================================

SUBROUTINE OutputFinalWake( BladeQuarterChord1, BladeQuarterChord2 )

  USE MultTurb_Params, Only: NWake, FWake
  USE FVW_Subs
  USE FVW_Parm, Only: CUTOFF_prim, CUTOFF_upmax, CUTOFF_upinit, NumBl, NumBS, NnearMax, I1, WakeAgeLimit

  IMPLICIT NONE

  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:) :: BladeQuarterChord1, BladeQuarterChord2
  INTEGER :: j2, kindx, nindx, nindx2, kindx2, kindx3, jindx, nbs, cut_low, cut_high
  CHARACTER(10) :: File1!, File2

  File1 = 'Turb1'!; File2 = 'Turb2' !HACK!!
  cut_low = 1; cut_high = CUTOFF_upmax(1)

  OPEN( unit = 1000, file = ( 'Turb1_r_primej.txt_new'   ), STATUS = 'new' )
  OPEN( unit = 1001, file = ( 'Turb1_r_primejm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1002, file = ( 'Turb1_r_primejm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1003, file = ( 'Turb1_r_primejm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1025, file = ( 'Turb1_r_oldj.txt_new'   ), STATUS = 'new' )
  OPEN( unit = 1004, file = ( 'Turb1_r_oldjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1005, file = ( 'Turb1_r_oldjm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1006, file = ( 'Turb1_r_oldjm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1007, file = ( 'Turb1_r_newjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1008, file = ( 'Turb1_r_newjm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1009, file = ( 'Turb1_r_newjm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1010, file = ( 'Turb1_rjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1011, file = ( 'Turb1_rjm2.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1012, file = ( 'Turb1_r_nearjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1013, file = ( 'Turb1_r_nearjm2.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1017, file = ( 'Turb1_Gammablj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 1018, file = ( 'Turb1_Gammabljm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1019, file = ( 'Turb1_Gammaj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 1020, file = ( 'Turb1_Gammajm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1021, file = ( 'Turb1_Gamma_nearj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 1022, file = ( 'Turb1_Gamma_nearjm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 1023, file = ( 'Turb1_BladeQuarterChordjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 1024, file = ( 'Turb1_BladeQuarterChordjm2.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2000, file = ( 'Turb2_r_primej.txt_new'   ), STATUS = 'new' )
  OPEN( unit = 2001, file = ( 'Turb2_r_primejm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2002, file = ( 'Turb2_r_primejm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2003, file = ( 'Turb2_r_primejm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2025, file = ( 'Turb2_r_oldj.txt_new'   ), STATUS = 'new' )
  OPEN( unit = 2004, file = ( 'Turb2_r_oldjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2005, file = ( 'Turb2_r_oldjm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2006, file = ( 'Turb2_r_oldjm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2007, file = ( 'Turb2_r_newjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2008, file = ( 'Turb2_r_newjm2.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2009, file = ( 'Turb2_r_newjm3.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2010, file = ( 'Turb2_rjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2011, file = ( 'Turb2_rjm2.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2012, file = ( 'Turb2_r_nearjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2013, file = ( 'Turb2_r_nearjm2.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2017, file = ( 'Turb2_Gammablj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 2018, file = ( 'Turb2_Gammabljm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2019, file = ( 'Turb2_Gammaj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 2020, file = ( 'Turb2_Gammajm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2021, file = ( 'Turb2_Gamma_nearj.txt_new' ),   STATUS = 'new' )
  OPEN( unit = 2022, file = ( 'Turb2_Gamma_nearjm1.txt_new' ), STATUS = 'new' )

  OPEN( unit = 2023, file = ( 'Turb2_BladeQuarterChordjm1.txt_new' ), STATUS = 'new' )
  OPEN( unit = 2024, file = ( 'Turb2_BladeQuarterChordjm2.txt_new' ), STATUS = 'new' )

  WRITE ( 1000, * ) ((( FWake%r_primej(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1001, * ) ((( FWake%r_primejm1( j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1002, * ) ((( FWake%r_primejm2( j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1003, * ) ((( FWake%r_primejm3( j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1025, * ) ((( FWake%r_oldj(     j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1004, * ) ((( FWake%r_oldjm1(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1005, * ) ((( FWake%r_oldjm2(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1006, * ) ((( FWake%r_oldjm3(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1007, * ) ((( FWake%r_newjm1(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1008, * ) ((( FWake%r_newjm2(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1009, * ) ((( FWake%r_newjm3(   j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1010, * ) ((( FWake%rjm1(       j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )
  WRITE ( 1011, * ) ((( FWake%rjm2(       j2, kindx, nindx ), j2=1,3),kindx=1,cut_high), nindx=1,NumBl )

  WRITE ( 1012, * ) (((( NWake%r_nearj( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )
  WRITE ( 1013, * ) (((( NWake%r_nearjm1( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )

  WRITE ( 1017, * ) ((NWake%Gammabljm1(   nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl ) !b/c Gammabl set to 0.0 in Update Aero   9.30.16
  WRITE ( 1018, * ) ((NWake%Gammabljm1( nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )

  WRITE ( 1019, * ) ((FWake%Gammajp1(   kindx, nindx ), kindx = 1, cut_high), nindx = 1, NumBl )
  WRITE ( 1020, * ) ((FWake%Gammaj( kindx, nindx ), kindx = 1, cut_high), nindx = 1, NumBl )

  WRITE ( 1021, * ) (((NWake%Gamma_nearjp1( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )
  WRITE ( 1022, * ) (((NWake%Gamma_nearj( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )

  WRITE ( 1023, * ) ((( BladeQuarterChord1( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )
  WRITE ( 1024, * ) ((( BladeQuarterChord2( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )

  WRITE ( 2000, * ) ((( FWake%r_primej(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2001, * ) ((( FWake%r_primejm1(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2002, * ) ((( FWake%r_primejm2( j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2003, * ) ((( FWake%r_primejm3( j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2025, * ) ((( FWake%r_oldj(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2004, * ) ((( FWake%r_oldjm1(     j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2005, * ) ((( FWake%r_oldjm2(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2006, * ) ((( FWake%r_oldjm3(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2007, * ) ((( FWake%r_newjm1(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2008, * ) ((( FWake%r_newjm2(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2009, * ) ((( FWake%r_newjm3(   j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2010, * ) ((( FWake%rjm1(       j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )
  WRITE ( 2011, * ) ((( FWake%rjm2(       j2, kindx, nindx ), j2=1,3),kindx=1,CUTOFF_prim), nindx=1,NumBl )

  WRITE ( 2012, * ) (((( NWake%r_nearj( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )
  WRITE ( 2013, * ) (((( NWake%r_nearjm1( j2, nbs, kindx3, nindx ) , j2 = 1, 3 ), nbs = 1, NumBS+1 ), kindx3 = 1, NnearMax ), nindx = 1, NumBl )

  WRITE ( 2017, * ) ((NWake%Gammabljm1(   nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl ) !b/c Gammabl set to 0.0 in Update Aero   9.30.16
  WRITE ( 2018, * ) ((NWake%Gammabljm1( nbs, nindx ), nbs = 1, NumBS ), nindx = 1, NumBl )

  WRITE ( 2019, * ) ((FWake%Gammajp1(   kindx, nindx ), kindx = 1, CUTOFF_prim), nindx = 1, NumBl )
  WRITE ( 2020, * ) ((FWake%Gammaj( kindx, nindx ), kindx = 1, CUTOFF_prim), nindx = 1, NumBl )

  WRITE ( 2021, * ) (((NWake%Gamma_nearjp1( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )
  WRITE ( 2022, * ) (((NWake%Gamma_nearj( kindx3,nbs,nindx ),kindx3=1,NnearMax),nbs=1,NumBS+1),nindx=1,NumBl )

  WRITE ( 2023, * ) ((( BladeQuarterChord1( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )
  WRITE ( 2024, * ) ((( BladeQuarterChord2( j2, nbs, nindx ), j2=1, 3 ), nbs=1, NumBS+1 ), nindx=1, NumBl )

  CLOSE( 1001 ); CLOSE( 1002 ); CLOSE( 1003 ); CLOSE( 1004 ); CLOSE( 1005 ); CLOSE( 1006 )
  CLOSE( 1007 ); CLOSE( 1008 ); CLOSE( 1009 ); CLOSE( 1010 ); CLOSE( 1011 ); CLOSE( 1012 )
  CLOSE( 1013 ); CLOSE( 1014 ); CLOSE( 1015 ); CLOSE( 1016 ); CLOSE( 1017 ); CLOSE( 1018 )
  CLOSE( 1019 ); CLOSE( 1020 ); CLOSE( 1021 ); CLOSE( 1022 ); CLOSE( 1023 ); CLOSE( 1024 )
  CLOSE( 1000 ); CLOSE( 1025 )

  CLOSE( 2001 ); CLOSE( 2002 ); CLOSE( 2003 ); CLOSE( 2004 ); CLOSE( 2005 ); CLOSE( 2006 )
  CLOSE( 2007 ); CLOSE( 2008 ); CLOSE( 2009 ); CLOSE( 2010 ); CLOSE( 2011 ); CLOSE( 2012 )
  CLOSE( 2013 ); CLOSE( 2014 ); CLOSE( 2015 ); CLOSE( 2016 ); CLOSE( 2017 ); CLOSE( 2018 )
  CLOSE( 2019 ); CLOSE( 2020 ); CLOSE( 2021 ); CLOSE( 2022 ); CLOSE( 2023 ); CLOSE( 2024 )
  CLOSE( 2000 ); CLOSE( 2025 )


END SUBROUTINE OutputFinalWake
!=================================================

!=================================================
SUBROUTINE OpenFiles

  !*****************************
  !This subroutine opens the appropriate files to which turbine wake information
  !will be written. The file names will be determines by which turbine is
  !currently being calculated
  !
  !      NOTE: Initialization files will still be opened/closed in FVW_INITIALIZE_WAKE.f90
  !
  !      Kelsey Shaler 8/27/2014
  !*****************************

  USE MultTurb_Params, Only: Turbines, NTurb
  USE FVW_Subs
  USE FVW_Parm, Only: NumBl

  IMPLICIT NONE

  INTEGER :: IUNIT

  !FileRoot2 = TRIM( Turbines( NTurb )%TurbNames )

  OPEN( unit = 101+NTurb, file = ( TRIM( 'Turb1' )//'_valuesrn1.txt'       ), STATUS = 'new' )
  OPEN( unit = 201+NTurb, file = ( TRIM( 'Turb1' )//'_valuesrn2.txt'       ), STATUS = 'new' )

  IF (NumBl .EQ. 3 ) Then
    OPEN( unit = 301+NTurb, file = ( TRIM( 'Turb1' )//'_valuesrn3.txt' ), STATUS = 'new' )
  END IF


END SUBROUTINE OpenFiles
!=================================================

!=================================================
SUBROUTINE WakeVelProfile(zloc, Wind_FVW, jold, rm1, Gammam1, rm2, r_nearm1, Gamma_nearm1, r_nearm2, &
                         & BladeQCm1, Gammablm1, BladeQCm2 )

  !*****************************
  ! This subroutine writes out the wake velocity profile at locations
  ! specified in 'InputFiles/WakePoints.txt'. The velocity profiles will
  ! be written to NumTurb_WakeVel_Loc#.txt
  !
  !      Kelsey Shaler 10.6.15
  !*****************************

  USE FVW_Parm
  USE FVW_Subs
  USE InflowWind
  USE AeroDyn14_Types, Only: FVW_WindType
  USE MultTurb_Params,    Only: NTurb, FWake, TurbLocs!GCoord
  USE FVW_ComputeWake
  USE InflowWind_Subs, Only: CalculateOutput

  IMPLICIT NONE

  INTEGER( IntKi ),                                          INTENT( IN    ) :: jold
  REAL( ReKi ),                                              INTENT( IN    ) :: zloc
  REAL( ReKi ), DIMENSION( NumBS, NumWakes                ), INTENT( IN    ) :: Gammablm1
  REAL( ReKi ), DIMENSION( CUTOFF_Allocate, NumWakes      ), INTENT( IN    ) :: Gammam1
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NumWakes           ), INTENT( IN    ) :: BladeQCm1, BladeQCm2
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes   ), INTENT( IN    ) :: rm1, rm2
  REAL( ReKi ), DIMENSION( NnearMax, NumBS+1, NumWakes    ), INTENT( IN    ) :: Gamma_nearm1
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NnearMax, NumWakes ), INTENT( IN    ) :: r_nearm1, r_nearm2
  TYPE( FVW_WindType ),                                      INTENT( INOUT ) :: Wind_FVW

  INTEGER( IntKi ),                                SAVE :: NumPtsz, NumPtsy, WakeCounter = 1
  REAL( ReKi ), DIMENSION( :, :, : ), ALLOCATABLE, SAVE :: WakeOutput, WakeVel, VinducedAvgTot

  INTEGER( IntKi )                            :: I, J, ErrStat, WriteUnit_Int, n, m, q
  REAL( ReKi )                                :: y_limit, ystep, yloc, Hub, VinducedAvgNorm, VinducedNorm, WhichTurb
  REAL( ReKi ), DIMENSION( 3 )                :: tmpvector, VinducedAvg
  REAL( ReKi ), DIMENSION( : ), ALLOCATABLE   :: zloc_Wake
  CHARACTER(3)                                :: ErrorMsg, WriteUnit_Char
  CHARACTER( * ), PARAMETER                   :: root = 'InputFiles/'

  VinducedFW1 = 0.00_ReKi; VinducedNW1 = 0.00_ReKi; VinducedBC1 = 0.00_ReKi; VinducedTot1 = 0.00_ReKi
  tmpvector = 0.00_ReKi

! Reading file that contains points for outputting wake velocities
  IF ( WakeCounter .EQ. 1 ) THEN

     OPEN( unit = 1, file = root//'WakePoints.txt' )

     READ( 1, * ) NumPtsz
     READ( 1, * ) NumPtsy
     ALLOCATE( zloc_Wake( NumPtsz ), WakeOutput( 3, NumPtsy, NumPtsz ), WakeVel( 3, NumPtsy, NumPtsz ), &
             & VinducedAvgTot( 3, NumPtsy, NumPtsz ))
     VinducedAvgTot = 0.00_ReKi

     zloc_Wake = 0.00_ReKi; WakeOutput = 0.00_ReKi
     READ( 1, * ) (zloc_Wake(I), I = 1, NumPtsz )
     READ( 1, * ) Hub
     CLOSE( 1 )
     WakeOutput( 1, :, : ) = 0.00_ReKi    ! setting x=HH for all z-locs
PRINT*, 'Radius:',Radius
     y_limit = 2.0*Radius + Hub ; ystep = 2.0_ReKi * y_limit/dble(NumPtsy)
PRINT*, 'y_limit: ', y_limit, 'ystep: ', ystep
     DO I = 1, NumPtsz
        yloc = -y_limit
        DO J = 1, NumPtsy
           WakeOutput( 2, J, I ) = yloc!/Radius
           WakeOutput( 3, J, I ) = zloc_Wake( I ) + TurbLocs( 1, 2 )
           yloc = yloc + ystep
        END DO
        WriteUnit_Int = I+5030+NTurb
        Write( WriteUnit_Char, "(I0.3)") INT(zloc_Wake( I ))
        OPEN( unit = WriteUnit_Int, file = (TRIM( 'Turb1' )//'_z='//TRIM(WriteUnit_Char)//'_WakeVel.txt'))
        WRITE( WriteUnit_Int, * ), '  VARIABLES= "y-location (m)", "x-velocity", "y-velocity", "z-velocity", "Uinf-normalized z-velocity"'
     END DO

     DEALLOCATE( zloc_Wake )
  END IF

!Writing out wake velocities for FreestreamVel
  WakeVel = 0.00_ReKi; VinducedAvg = 0.00_ReKi

  DO I = 1, NumPtsz
     DO J = 1, NumPtsy
        tmpvector = WakeOutput( :, J, I )

        CALL TRANSFORM_TO_AERODYN_COORDS( tmpvector, zloc )

        Wind_FVW%InputData%PositionXYZ( :, 1 ) = tmpvector
        CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )

        WakeVel( :, J, I ) = Wind_FVW%OutputData%VelocityUVW( :, 1 )

        CALL TRANSFORM_TO_FVW_COORDS( WakeVel(:, J, I ))

     END DO
  END DO

  DO J = 1, NumPtsz
     WriteUnit_Int = J + 5000+NTurb
     DO I = 1, NumPtsy
        CALL Vinduced2PRIME( rm1, Gammam1, WakeOutput( :, I, J ), &
                                 & rm2, 1, jold, 1 ) !Vind due to far wake

        VinducedFW1( 1 ) = Sum( FWake%VinducedFarWakej( 1, 1, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes ))
        VinducedFW1( 2 ) = Sum( FWake%VinducedFarWakej( 2, 1, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes ))
        VinducedFW1( 3 ) = Sum( FWake%VinducedFarWakej( 3, 1, 1, NnearMax+1:WakeAgeLimit, 1:NumWakes ))

        CALL VinducedNW( r_nearm1, Gamma_nearm1, WakeOutput( :, I, J ), VinducedNW1, &
          & r_nearm2, NumWakes ) !Vind due to near wake

        CALL VinducedBC( BladeQCm1, Gammablm1, WakeOutput( :, I, J ), &
          & VinducedBC1 ) !Vind due to blades

        VinducedTot1 =  VinducedFW1 + VinducedNW1 + VinducedBC1 + WakeVel( :, I, J )
            WRITE( WriteUnit_Int, 100 ) WakeOutput( 2, I, J )/(Radius+Hub), VinducedTot1/WakeVel(:,I,J)
     END DO
  END DO
100  FORMAT( 16F12.7 )
  WakeCounter = WakeCounter+1

END SUBROUTINE WakeVelProfile
!=================================================

END MODULE FileManipulation
