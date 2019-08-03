!FIXME: all of this should be put into a registry.  Nothing should have a 'SAVE' on it
!==========================================================================
!==========================================================================

MODULE FVW_Parm

  USE NWTC_Library 
  USE FVW_Types
  USE AD14AeroConf_Types
!FIXME: check if any of thes are inputs, or get changed.  Move to inputtype or miscvars accordingly.  Put rest in parameters
  INTEGER(IntKi) :: CUTOFF_prim, PerOverlap, CUTOFF_Allocate, NumBl, WakeAgeLimit
  INTEGER(IntKi) :: NumBS, Nj, Nj2, NnearMax, NElm, Nelm_start, Num_start, I1
  INTEGER(IntKi), ALLOCATABLE, DIMENSION(:), SAVE :: CUTOFF, CUTOFF_upinit, CUTOFF_upmax, CUTOFF_up, BC

  REAL(ReKi) :: Radius, HubHt, HubR, DtAero
  REAL(ReKi) :: Root_cut, eps, nu, near_deg, delta_psi_Est, Omega, Rad, dRad, TMax, RotSpeed, RotSpeed_Est
  REAL(DbKi) :: Time_Real
  REAL(ReKi), DIMENSION(:), ALLOCATABLE, SAVE  :: RELM
  REAL(ReKi), PARAMETER :: alpha_param=1.256430_ReKi, a1=0.00020_ReKi

  TYPE(AD14AeroConf_ParameterType) :: FVW_AirfoilParm
  TYPE(AD14AeroConf_MiscVarType) :: FVW_AirfoilOut
  REAL(ReKi), DIMENSION(:), ALLOCATABLE, SAVE :: FVW_CDD, FVW_CLL, FVW_CMM
  REAL(ReKi), DIMENSION(:), ALLOCATABLE, SAVE :: Chord
  REAL(ReKi), DIMENSION(3), SAVE :: delta_psi = 0.00_ReKi

END MODULE FVW_Parm

!==========================================================================
!==========================================================================

MODULE FVW_ComputeWake

USE Precision
!FIXME: perhaps create a new type specifically for the wake calculations.  Check if that makes sense given how everythin is handled internally.  Maybe put this type inside miscvars or another type?????
  INTEGER :: IElement, IBlade, nbsindx, counter, init, high, counter2=0

  REAL(ReKi), SAVE :: VN, VT, dx, SPitch, CPitch, Pitnow, TurbLength=0.00_ReKi

  REAL(ReKi), ALLOCATABLE, DIMENSION(:  ), SAVE     :: C1, Velsec, Velsec2, C2
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:), SAVE   :: velstorej
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:), SAVE   :: a_of_a_storej, a_of_a_effective

  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: BladeTanVectj, BladeQuarterChordj, BladeLocj
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: BladeQuarterChordjm1, BladeQuarterChordjm2
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: BladeNormVect2j, BladeTanVect2j, Vind_storej
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: BladeLoc2j, BladeThreeQuarterChordj, VinducedNWFinal

  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: VinducedFW1, VinducedNW1, VinducedBC1
  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: VinducedTot1, VinducedTot2, VinducedTot1b, VinducedTot2b
  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: VinducedNWtest

  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: rpcyl1, rpcyl2, rpcyl3, rpcyl4, rpcyl5, rpcyl6, rpcyl7, rpcyl8
  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: rncyl1, rncyl2, rncyl3, rncyl4, rncyl5, rncyl6, rncyl7, rncyl8
  REAL(ReKi), ALLOCATABLE, DIMENSION(:), SAVE     :: rocyl1, rocyl2, rocyl3, rocyl4, rocyl5, rocyl6, rocyl7, rocyl8

  REAL(ReKi), ALLOCATABLE, DIMENSION(:  ), SAVE   :: CalcedVinf
  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:), SAVE   :: VindTotal

  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: Vaxial2j, VNElem2j, Vaxialj, VTT, VNElement

  REAL(ReKi), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: BladeLoc2j_Real, r_oldj_Real, r_primej_Real

END MODULE FVW_ComputeWake

!==========================================================================
