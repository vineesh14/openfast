MODULE MultTurb_Params

!This subroutine contains variables used in the multiple turbine simulations but not in the single turbine simulations

  USE precision

IMPLICIT NONE

  LOGICAL :: MultTurbs
  INTEGER :: TurbAlignment, NTurb, NumWakes
  !INTEGER :: mype, npes
  CHARACTER*10 :: ITurb
  CHARACTER*30 :: WindFileNameKS!, FileRoot1, FileRoot2, FileRoot3

  REAL(ReKi) :: TurbDist1, RotRadius, PerUinf

  REAL( ReKi ), ALLOCATABLE, DIMENSION( :, : ) :: TurbLocs
  REAL( ReKi ), ALLOCATABLE, DIMENSION( :    ) :: TurbDist

  !For column relations
  INTEGER :: IOutWake=1

! =========  NearWake  =======
  TYPE  INearWake
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :    ) :: r_nearj, r_nearjm1, r_nearjm2
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :          ) :: Gammablj, Gammabljm1
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: Gamma_nearj, Gamma_nearjm1, Gamma_nearjp1
  END TYPE INearWake
! ============================
! =========  FarWake  =======
  TYPE  IFarWake
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: rj, rjm1, rjm2
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_primej, r_primejm1, r_primejm2, r_primejm3
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_oldj, r_oldjm1, r_oldjm2, r_oldjm3
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_newj, r_newjm1, r_newjm2, r_newjm3
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :          ) :: Gammaj, Gammajm1, Gammajp1
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :, : ) :: VinducedFarWakePrimej, VinducedFarWakePrimejm1
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :, : ) :: VinducedFarWakej, VinducedFarWakejm1
     REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :, : ) :: VinducedFarWakeRj, VinducedFarWakeRjm1
  END TYPE IFarWake
! ============================
! =========  GlobalCoords  =======
  !TYPE  GlobalCoords
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: rj, rjm1, rjm2
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_primej, r_primejm1, r_primejm2, r_primejm3
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_oldj, r_oldjm1, r_oldjm2, r_oldjm3
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: r_newj, r_newjm1, r_newjm2, r_newjm3
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :          ) :: Gammaj, Gammajm1, Gammajp1
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :, : ) :: VinducedFarWakej, VinducedFarWakejm1! VinducedFarWakeOldj

  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :    ) :: r_nearj, r_nearjm1, r_nearjm2
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :          ) :: Gammablj, Gammabljm1
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: Gamma_nearj, Gamma_nearjm1, Gamma_nearjp1
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :, :, : ) :: VinducedFarWakeRj, VinducedFarWakeRjm1
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: BladeThreeQuarterChordj, BladeQuarterChordj
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: BladeQuarterChordjm1, BladeQuarterChordjm2
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: BladeNormVect2j, BladeTanVect2j, BladeTanVectj
  !   REAL( ReKi ), ALLOCATABLE, DIMENSION( :, :, :       ) :: BladeLoc2j, BladeLoc2j_Real
  !END TYPE GlobalCoords
! ============================
! =========  IndvTurb  =======
  TYPE, PUBLIC :: IndvTurb
     CHARACTER*10 :: TurbNames
     REAL( ReKi ), DIMENSION(2) :: Location
     CHARACTER*30 :: WindFileName
     REAL( ReKi ) :: RowNum
  END TYPE IndvTurb
! ============================

  TYPE( IndvTurb         ), ALLOCATABLE, DIMENSION(:), SAVE :: Turbines
  TYPE( IFarWake         )              :: FWake
  TYPE( INearWake        )              :: NWake
  !TYPE( GlobalCoords     )              :: GCoord
END MODULE MultTurb_Params
!========================================================
