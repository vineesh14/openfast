SUBROUTINE TRANSFORM_TO_AERODYN_COORDS(Vector,zloc)

  USE FVW_Parm, Only: HH
  USE precision

  IMPLICIT NONE

  REAL(ReKi), INTENT( IN    )               :: zloc
  REAL(ReKi), INTENT(   OUT ), DIMENSION(3) :: Vector
  REAL(ReKi), DIMENSION(:), ALLOCATABLE     :: TmpVector
  REAL(ReKi), DIMENSION(:,:), ALLOCATABLE   :: Transf

  ALLOCATE(TmpVector(3),Transf(3,3))

  Transf(1,:)=(/0.00_ReKi,0.00_ReKi,1.00_ReKi/)
  Transf(2,:)=(/0.00_ReKi,-1.00_ReKi,0.00_ReKi/)
  Transf(3,:)=(/1.00_ReKi,0.00_ReKi,0.00_ReKi/)
  TmpVector=MatMul(Transf,Vector)

  Vector=TmpVector

  DEALLOCATE(TmpVector,Transf)

  Vector(1)=Vector(1)+zloc
  Vector(3)=Vector(3)+HH
  Vector(1)=0.0_ReKi!  KS--I don't even understand why that was there.... 3.10.15

END SUBROUTINE TRANSFORM_TO_AERODYN_COORDS
! ==================================================================

! ==================================================================
SUBROUTINE TRANSFORM_TO_FVW_COORDS(Vector)

  USE precision

  IMPLICIT NONE

  REAL(ReKi), DIMENSION(3) :: Vector, TmpVector
  REAL(ReKi), DIMENSION(3,3) :: Transf


  Transf(1,:)=(/0.00_ReKi,0.00_ReKi,1.00_ReKi/)
  Transf(2,:)=(/0.00_ReKi,-1.00_ReKi,0.00_ReKi/)
  Transf(3,:)=(/1.00_ReKi,0.00_ReKi,0.00_ReKi/)
  TmpVector=MatMul(Transf,Vector)

  Vector=TmpVector
END SUBROUTINE TRANSFORM_TO_FVW_COORDS
! ==================================================================
! ==================================================================
!SUBROUTINE TransformToGlobal(  )
!
!  USE MultTurb_Params,  Only : FWake, NWake, GCoord, Turbines, NTurb, NumWakes, TurbLocs!,mype
!  USE FVW_Parm, Only: NumBl, CUTOFF_upmax, CUTOFF_upinit
!  USE FVW_ComputeWake, Only: BladeThreeQuarterChordj, BladeQuarterChordj, BladeQuarterChordjm1, BladeQuarterChordjm2, BladeLoc2j, BladeLoc2j_Real
!
!  IMPLICIT NONE
!
!  INTEGER :: llow, lhigh, glow, ghigh, lcut_low, lcut_high, gcut_low, gcut_high
!
!  NTurb = 1 !!!KS -- HACKHACKHACK 7.1.19
!  !IF ( mype .EQ. 0 ) THEN
!     llow = 1; lhigh = NumBl; glow = 1; ghigh = NumBl
!
!     lcut_low = 1; lcut_high = CUTOFF_upmax(1)
!     gcut_low = 1; gcut_high = CUTOFF_upmax(1)
!  !ELSE
!  !   llow = 1; lhigh = NumBl; ghigh = (mype+1)*NumBl; glow = ghigh-(NumBl-1)
!  !   lcut_low = 1; lcut_high = CUTOFF_upmax(mype+1)
!  !   gcut_low = 1; gcut_high = CUTOFF_upmax(mype+1)
!
!  !END IF
!
!  GCoord%rj(         1, :, glow:ghigh ) = FWake%rj(         1, :, llow:lhigh )
!  GCoord%rjm1(       1, :, glow:ghigh ) = FWake%rjm1(       1, :, llow:lhigh )
!  GCoord%rjm2(       1, :, glow:ghigh ) = FWake%rjm2(       1, :, llow:lhigh )
!  GCoord%r_primej(   1, :, glow:ghigh ) = FWake%r_primej(   1, :, llow:lhigh )
!  GCoord%r_primejm1( 1, :, glow:ghigh ) = FWake%r_primejm1( 1, :, llow:lhigh )
!  GCoord%r_primejm2( 1, :, glow:ghigh ) = FWake%r_primejm2( 1, :, llow:lhigh )
!  GCoord%r_primejm3( 1, :, glow:ghigh ) = FWake%r_primejm3( 1, :, llow:lhigh )
!  GCoord%r_oldj(     1, :, glow:ghigh ) = FWake%r_oldj(     1, :, llow:lhigh )
!  GCoord%r_oldjm1(   1, :, glow:ghigh ) = FWake%r_oldjm1(   1, :, llow:lhigh )
!  GCoord%r_oldjm2(   1, :, glow:ghigh ) = FWake%r_oldjm2(   1, :, llow:lhigh )
!  GCoord%r_oldjm3(   1, :, glow:ghigh ) = FWake%r_oldjm3(   1, :, llow:lhigh )
!  GCoord%r_newj(     1, :, glow:ghigh ) = FWake%r_newj(     1, :, llow:lhigh )
!  GCoord%r_newjm1(   1, :, glow:ghigh ) = FWake%r_newjm1(   1, :, llow:lhigh )
!  GCoord%r_newjm2(   1, :, glow:ghigh ) = FWake%r_newjm2(   1, :, llow:lhigh )
!  GCoord%r_newjm3(   1, :, glow:ghigh ) = FWake%r_newjm3(   1, :, llow:lhigh )
!
!  GCoord%rj(         2, :, glow:ghigh ) = FWake%rj(         2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%rjm1(       2, :, glow:ghigh ) = FWake%rjm1(       2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%rjm2(       2, :, glow:ghigh ) = FWake%rjm2(       2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_primej(   2, :, glow:ghigh ) = FWake%r_primej(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_primejm1( 2, :, glow:ghigh ) = FWake%r_primejm1( 2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_primejm2( 2, :, glow:ghigh ) = FWake%r_primejm2( 2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_primejm3( 2, :, glow:ghigh ) = FWake%r_primejm3( 2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_oldj(     2, :, glow:ghigh ) = FWake%r_oldj(     2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_oldjm1(   2, :, glow:ghigh ) = FWake%r_oldjm1(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_oldjm2(   2, :, glow:ghigh ) = FWake%r_oldjm2(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_oldjm3(   2, :, glow:ghigh ) = FWake%r_oldjm3(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_newj(     2, :, glow:ghigh ) = FWake%r_newj(     2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_newjm1(   2, :, glow:ghigh ) = FWake%r_newjm1(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_newjm2(   2, :, glow:ghigh ) = FWake%r_newjm2(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!  GCoord%r_newjm3(   2, :, glow:ghigh ) = FWake%r_newjm3(   2, :, llow:lhigh ) + TurbLocs(NTurb,1)
!
!  GCoord%rj(         3, :, glow:ghigh ) = FWake%rj(         3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%rjm1(       3, :, glow:ghigh ) = FWake%rjm1(       3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%rjm2(       3, :, glow:ghigh ) = FWake%rjm2(       3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_primej(   3, :, glow:ghigh ) = FWake%r_primej(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_primejm1( 3, :, glow:ghigh ) = FWake%r_primejm1( 3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_primejm2( 3, :, glow:ghigh ) = FWake%r_primejm2( 3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_primejm3( 3, :, glow:ghigh ) = FWake%r_primejm3( 3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_oldj(     3, :, glow:ghigh ) = FWake%r_oldj(     3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_oldjm1(   3, :, glow:ghigh ) = FWake%r_oldjm1(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_oldjm2(   3, :, glow:ghigh ) = FWake%r_oldjm2(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_oldjm3(   3, :, glow:ghigh ) = FWake%r_oldjm3(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_newj(     3, :, glow:ghigh ) = FWake%r_newj(     3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_newjm1(   3, :, glow:ghigh ) = FWake%r_newjm1(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_newjm2(   3, :, glow:ghigh ) = FWake%r_newjm2(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!  GCoord%r_newjm3(   3, :, glow:ghigh ) = FWake%r_newjm3(   3, :, llow:lhigh ) + TurbLocs(NTurb,2)
!
!  GCoord%r_nearj(  1, :, :, glow:ghigh ) = NWake%r_nearj(   1, :, :, : )
!  GCoord%r_nearjm1(1, :, :, glow:ghigh ) = NWake%r_nearjm1( 1, :, :, : )
!  GCoord%r_nearjm2(1, :, :, glow:ghigh ) = NWake%r_nearjm2( 1, :, :, : )
!
!  GCoord%r_nearj(  2, :, :, glow:ghigh ) = NWake%r_nearj(   2, :, :, : ) + TurbLocs(NTurb,1)
!  GCoord%r_nearjm1(2, :, :, glow:ghigh ) = NWake%r_nearjm1( 2, :, :, : ) + TurbLocs(NTurb,1)
!  GCoord%r_nearjm2(2, :, :, glow:ghigh ) = NWake%r_nearjm2( 2, :, :, : ) + TurbLocs(NTurb,1)
!
!  GCoord%r_nearj(  3, :, :, glow:ghigh ) = NWake%r_nearj(   3, :, :, : ) + TurbLocs(NTurb,2)
!  GCoord%r_nearjm1(3, :, :, glow:ghigh ) = NWake%r_nearjm1( 3, :, :, : ) + TurbLocs(NTurb,2)
!  GCoord%r_nearjm2(3, :, :, glow:ghigh ) = NWake%r_nearjm2( 3, :, :, : ) + TurbLocs(NTurb,2)
!
!  GCoord%BladeThreeQuarterChordj( 1, :, : ) = BladeThreeQuarterChordj( 1, :, : )
!  GCoord%BladeLoc2j_Real(         1, :, : ) = BladeLoc2j_Real(         1, :, : )
!  GCoord%BladeThreeQuarterChordj( 2, :, : ) = BladeThreeQuarterChordj( 2, :, : ) + TurbLocs(NTurb,1)
!  GCoord%BladeLoc2j_Real(         2, :, : ) = BladeLoc2j_Real(         2, :, : ) + TurbLocs(NTurb,1)
!  GCoord%BladeThreeQuarterChordj( 3, :, : ) = BladeThreeQuarterChordj( 3, :, : ) + TurbLocs(NTurb,2)
!  GCoord%BladeLoc2j_Real(         3, :, : ) = BladeLoc2j_Real(         3, :, : ) + TurbLocs(NTurb,2)
!
!  GCoord%BladeQuarterChordj(   1, :, glow:ghigh ) = BladeQuarterChordj(   1, :, : )
!  GCoord%BladeQuarterChordjm1( 1, :, glow:ghigh ) = BladeQuarterChordjm1( 1, :, : )
!  GCoord%BladeQuarterChordjm2( 1, :, glow:ghigh ) = BladeQuarterChordjm2( 1, :, : )
!
!  GCoord%BladeQuarterChordj(   2, :, glow:ghigh ) = BladeQuarterChordj(   2, :, : ) + TurbLocs(NTurb,1)
!  GCoord%BladeQuarterChordjm1( 2, :, glow:ghigh ) = BladeQuarterChordjm1( 2, :, : ) + TurbLocs(NTurb,1)
!  GCoord%BladeQuarterChordjm2( 2, :, glow:ghigh ) = BladeQuarterChordjm2( 2, :, : ) + TurbLocs(NTurb,1)
!
!  GCoord%BladeQuarterChordj(   3, :, glow:ghigh ) = BladeQuarterChordj(   3, :, : ) + TurbLocs(NTurb,2)
!  GCoord%BladeQuarterChordjm1( 3, :, glow:ghigh ) = BladeQuarterChordjm1( 3, :, : ) + TurbLocs(NTurb,2)
!  GCoord%BladeQuarterChordjm2( 3, :, glow:ghigh ) = BladeQuarterChordjm2( 3, :, : ) + TurbLocs(NTurb,2)
!
!  GCoord%Gammaj(   :, glow:ghigh ) = FWake%Gammaj(   :, llow:lhigh )
!  GCoord%Gammajm1( :, glow:ghigh ) = FWake%Gammajm1( :, llow:lhigh )
!  GCoord%Gammajp1( :, glow:ghigh ) = FWake%Gammajp1( :, llow:lhigh )
!
!  GCoord%Gammablj(      :,    glow:ghigh ) = NWake%Gammablj(      :,    : )
!  GCoord%Gammabljm1(    :,    glow:ghigh ) = NWake%Gammabljm1(    :,    : )
!  GCoord%Gamma_nearj(   :, :, glow:ghigh ) = NWake%Gamma_nearj(   :, :, : )
!  GCoord%Gamma_nearjm1( :, :, glow:ghigh ) = NWake%Gamma_nearjm1( :, :, : )
!  GCoord%Gamma_nearjp1( :, :, glow:ghigh ) = NWake%Gamma_nearjp1( :, :, : )
!
!END SUBROUTINE TransformToGlobal
!! ==================================================================
!! ==================================================================
!SUBROUTINE TransformToLocal(    )
!
!  USE MultTurb_Params,  Only : FWake, NWake, GCoord, Turbines, NTurb, NumWakes, TurbLocs!,mype
!  USE FVW_Parm, Only: NumBl, CUTOFF_upmax, CUTOFF_upinit
!  USE FVW_ComputeWake, Only: BladeThreeQuarterChordj, BladeQuarterChordj, BladeQuarterChordjm1, BladeQuarterChordjm2
!
!  IMPLICIT NONE
!
!  INTEGER :: glow, ghigh, llow, lhigh, lcut_low, lcut_high, gcut_low, gcut_high
!
!  !IF ( mype .EQ. 0 ) THEN
!     llow = 1; lhigh = NumBl; glow = 1; ghigh = NumBl
!     lcut_low = 1; lcut_high = CUTOFF_upmax(1)
!     gcut_low = 1; gcut_high = CUTOFF_upmax(1)
!  !ELSE
!  !   llow = 1; lhigh = NumBl; ghigh = (mype+1)*NumBl; glow = ghigh-(NumBl-1)
!  !   lcut_low = 1; lcut_high = CUTOFF_upmax(mype+1)
!  !   gcut_low = 1; gcut_high = CUTOFF_upmax(mype+1)
!  !END IF
!
!  FWake%rj(         1, :, llow:lhigh ) = GCoord%rj(         1, :, glow:ghigh )
!  FWake%rjm1(       1, :, llow:lhigh ) = GCoord%rjm1(       1, :, glow:ghigh )
!  FWake%rjm2(       1, :, llow:lhigh ) = GCoord%rjm2(       1, :, glow:ghigh )
!  FWake%r_primej(   1, :, llow:lhigh ) = GCoord%r_primej(   1, :, glow:ghigh )
!  FWake%r_primejm1( 1, :, llow:lhigh ) = GCoord%r_primejm1( 1, :, glow:ghigh )
!  FWake%r_primejm2( 1, :, llow:lhigh ) = GCoord%r_primejm2( 1, :, glow:ghigh )
!  FWake%r_primejm3( 1, :, llow:lhigh ) = GCoord%r_primejm3( 1, :, glow:ghigh )
!  FWake%r_oldj(     1, :, llow:lhigh ) = GCoord%r_oldj(     1, :, glow:ghigh )
!  FWake%r_oldjm1(   1, :, llow:lhigh ) = GCoord%r_oldjm1(   1, :, glow:ghigh )
!  FWake%r_oldjm2(   1, :, llow:lhigh ) = GCoord%r_oldjm2(   1, :, glow:ghigh )
!  FWake%r_oldjm3(   1, :, llow:lhigh ) = GCoord%r_oldjm3(   1, :, glow:ghigh )
!  FWake%r_newj(     1, :, llow:lhigh ) = GCoord%r_newj(     1, :, glow:ghigh )
!  FWake%r_newjm1(   1, :, llow:lhigh ) = GCoord%r_newjm1(   1, :, glow:ghigh )
!  FWake%r_newjm2(   1, :, llow:lhigh ) = GCoord%r_newjm2(   1, :, glow:ghigh )
!  FWake%r_newjm3(   1, :, llow:lhigh ) = GCoord%r_newjm3(   1, :, glow:ghigh )
!
!  FWake%rj(         2, :, llow:lhigh ) = GCoord%rj(         2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%rjm1(       2, :, llow:lhigh ) = GCoord%rjm1(       2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%rjm2(       2, :, llow:lhigh ) = GCoord%rjm2(       2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_primej(   2, :, llow:lhigh ) = GCoord%r_primej(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_primejm1( 2, :, llow:lhigh ) = GCoord%r_primejm1( 2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_primejm2( 2, :, llow:lhigh ) = GCoord%r_primejm2( 2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_primejm3( 2, :, llow:lhigh ) = GCoord%r_primejm3( 2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_oldj(     2, :, llow:lhigh ) = GCoord%r_oldj(     2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_oldjm1(   2, :, llow:lhigh ) = GCoord%r_oldjm1(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_oldjm2(   2, :, llow:lhigh ) = GCoord%r_oldjm2(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_oldjm3(   2, :, llow:lhigh ) = GCoord%r_oldjm3(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_newj(     2, :, llow:lhigh ) = GCoord%r_newj(     2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_newjm1(   2, :, llow:lhigh ) = GCoord%r_newjm1(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_newjm2(   2, :, llow:lhigh ) = GCoord%r_newjm2(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  FWake%r_newjm3(   2, :, llow:lhigh ) = GCoord%r_newjm3(   2, :, glow:ghigh ) - TurbLocs(NTurb,1)
!
!  FWake%rj(         3, :, llow:lhigh ) = GCoord%rj(         3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%rjm1(       3, :, llow:lhigh ) = GCoord%rjm1(       3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%rjm2(       3, :, llow:lhigh ) = GCoord%rjm2(       3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_primej(   3, :, llow:lhigh ) = GCoord%r_primej(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_primejm1( 3, :, llow:lhigh ) = GCoord%r_primejm1( 3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_primejm2( 3, :, llow:lhigh ) = GCoord%r_primejm2( 3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_primejm3( 3, :, llow:lhigh ) = GCoord%r_primejm3( 3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_oldj(     3, :, llow:lhigh ) = GCoord%r_oldj(     3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_oldjm1(   3, :, llow:lhigh ) = GCoord%r_oldjm1(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_oldjm2(   3, :, llow:lhigh ) = GCoord%r_oldjm2(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_oldjm3(   3, :, llow:lhigh ) = GCoord%r_oldjm3(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_newj(     3, :, llow:lhigh ) = GCoord%r_newj(     3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_newjm1(   3, :, llow:lhigh ) = GCoord%r_newjm1(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_newjm2(   3, :, llow:lhigh ) = GCoord%r_newjm2(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  FWake%r_newjm3(   3, :, llow:lhigh ) = GCoord%r_newjm3(   3, :, glow:ghigh ) - TurbLocs(NTurb,2)
!
!  NWake%r_nearj(  1, :, :, : ) = GCoord%r_nearj(   1, :, :, glow:ghigh )
!  NWake%r_nearjm1(1, :, :,  : ) = GCoord%r_nearjm1( 1, :, :, glow:ghigh )
!  NWake%r_nearjm2(1, :, :,  : ) = GCoord%r_nearjm2( 1, :, :, glow:ghigh )
!
!  NWake%r_nearj(  2, :, :,  : ) = GCoord%r_nearj(   2, :, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  NWake%r_nearjm1(2, :, :,  : ) = GCoord%r_nearjm1( 2, :, :, glow:ghigh ) - TurbLocs(NTurb,1)
!  NWake%r_nearjm2(2, :, :,  : ) = GCoord%r_nearjm2( 2, :, :, glow:ghigh ) - TurbLocs(NTurb,1)
!
!  NWake%r_nearj(  3, :, :,  : ) = GCoord%r_nearj(   3, :, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  NWake%r_nearjm1(3, :, :,  : ) = GCoord%r_nearjm1( 3, :, :, glow:ghigh ) - TurbLocs(NTurb,2)
!  NWake%r_nearjm2(3, :, :,  : ) = GCoord%r_nearjm2( 3, :, :, glow:ghigh ) - TurbLocs(NTurb,2)
!
!  FWake%Gammaj(   :, llow:lhigh ) = GCoord%Gammaj(   :, glow:ghigh )
!  FWake%Gammajm1( :, llow:lhigh ) = GCoord%Gammajm1( :, glow:ghigh )
!  FWake%Gammajp1( :, llow:lhigh ) = GCoord%Gammajp1( :, glow:ghigh )
!
!  NWake%Gammablj(      :,     : ) = GCoord%Gammablj(      :,    glow:ghigh )
!  NWake%Gammabljm1(    :,     : ) = GCoord%Gammabljm1(    :,    glow:ghigh )
!  NWake%Gamma_nearj(   :, :,  : ) = GCoord%Gamma_nearj(   :, :, glow:ghigh )
!  NWake%Gamma_nearjm1( :, :,  : ) = GCoord%Gamma_nearjm1( :, :, glow:ghigh )
!  NWake%Gamma_nearjp1( :, :,  : ) = GCoord%Gamma_nearjp1( :, :, glow:ghigh )
!
!END SUBROUTINE TransformToLocal
! ==================================================================
