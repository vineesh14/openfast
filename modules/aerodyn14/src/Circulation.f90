! =====================================================================================
SUBROUTINE Calculate_Gamma1( n, VTotal, BladeTanVect, normalvector, BladeLoc, ControlPoints, Ainv, Cap_Gamma, &
                           & Gammabl, VortexPointsJmin1, VortexPoints, Gamma_near, zloc, VinducedNWFinal, Wind_FVW )
!

      ! ********************************************************
      ! This subroutine computes the circulation on the blades
      ! using the Weissinger - L model and the flow tangency condition
      !
      !       -- Description added by Kelsey Shaler
      ! ********************************************************

  USE FVW_Parm
  USE AeroDyn14_Types, Only: FVW_WindType
  USE InflowWind
  USE MathOps, Only: Dot, Pinv


  IMPLICIT NONE

  INTEGER,                                     INTENT( IN    ) :: n
  REAL( ReKi ),                                      INTENT( IN    ) :: zloc
  REAL( ReKi ), DIMENSION( 3, NumBS ),               INTENT( IN    ) :: ControlPoints, normalvector, VTotal
  REAL( ReKi ), DIMENSION( 3, NumBS + 1 ),           INTENT( IN    ) :: BladeTanVect
  REAL( ReKi ), DIMENSION( 3, NumBS + 1, NumBl ),    INTENT( IN    ) :: BladeLoc
  REAL( ReKi ), DIMENSION( 3, NumBS + 1, NnearMax ), INTENT( IN    ) :: VortexPointsJmin1

  REAL( ReKi ), DIMENSION( 3, NumBS + 1, NnearMax ), INTENT( INOUT ) :: Vortexpoints
  TYPE( FVW_WindType ),                        INTENT( INOUT ) :: Wind_FVW

  REAL( ReKi ),                                      INTENT(   OUT ) :: Cap_Gamma
  REAL( ReKi ), DIMENSION( NumBS ),                  INTENT(   OUT ) :: Gammabl
  REAL( ReKi ), DIMENSION( NumBS + 1 ),              INTENT(   OUT ) :: Gamma_near
  REAL( ReKi ), DIMENSION( 3, NumBS ),               INTENT(   OUT ) :: VinducedNWFinal
  REAL( ReKi ), DIMENSION( NumBS, NumBS ),           INTENT(   OUT ) :: Ainv

  INTEGER                                 :: indx2, indx1, m, jmax, ErrStat, nbs
  REAL( ReKi )                                  :: dr, Cap_Gamma_max, Cap_Gamma_min
  REAL( ReKi ), DIMENSION( 3                  ) :: SumBS, V
  REAL( ReKi ), DIMENSION( NumBS              ) :: Rnumbs, B
  REAL( ReKi ), DIMENSION( NumBS+1            ) :: Rnumbsp1
  REAL( ReKi ), DIMENSION( NumBS, NumBS       ) :: A
  REAL( ReKi ), DIMENSION( 3, ( 2*NnearMax-1 )) :: BS
  REAL( ReKi ), DIMENSION( 3, NumBS, NumBS    ) :: A2
  CHARACTER(            124             ) :: ErrorMsg

  REAL( ReKi ) :: TMP_Vect( 3 )

  dRad = Rad / dble( NumBS )

   !Splitting up blade into sections
  DO nbs = 1, NumBS
     Rnumbsp1( nbs ) = dRad * dble( nbs - 1 )
     Rnumbs( nbs ) = dRad * dble( nbs - 1 ) + dRad / 2.0_ReKi
  END DO
  Rnumbsp1( NumBS+1 ) = dRad * dble( NumBS )

  SumBS = 0.0_ReKi; BS = 0.0_ReKi

  dr = Rad / dble( NumBS )   !need b / c changing this based on dRad later
  A = 0.0_ReKi; Ainv = 0.0_ReKi; B = 0.0_ReKi; Gammabl = 0.0_ReKi; Gamma_near = 0.0_ReKi

  jmax = NnearMax - 1   !Setting limit for near wake; only calcualting this for near wake

  Vortexpoints( :, :, 1 ) = BladeLoc( :, :, n )   !Set vortex points to corresponding points on blade

  DO indx2 = Num_start, NumBS + 1   !Constraining induced velocity to be normal to blade plane
    dRad=(indx2-1)*dr
     Vortexpoints( :, indx2, 2 ) = BladeLoc( :, indx2, n ) + BladeTanVect( :, indx2 ) * &
       & dRad * TAN( delta_psi(1) )
  END DO

  IF ( jmax .GE. 2 ) THEN   !Setting up rest of vortex points !Is this ever not true?
     DO indx1 = 3, jmax + 1
        DO indx2 = Num_start, NumBS + 1
           TMP_Vect( : ) = VortexpointsJmin1( :, indx2, indx1-1 )
           CALL TRANSFORM_TO_AERODYN_COORDS( TMP_Vect, zloc )

           Wind_FVW%InputData%PositionXYZ( :, 1 ) = TMP_Vect

           CALL InflowWind_CalcOutput( Time_Real, Wind_FVW%InputData, Wind_FVW%ParamData, Wind_FVW%ContData, &
              & Wind_FVW%DiscData, Wind_FVW%ConstrData, Wind_FVW%OtherData, Wind_FVW%OutputData, &
              & Wind_FVW%MiscData, ErrStat, ErrorMsg )
           V = Wind_FVW%OutputData%VelocityUVW( :, 1 )
           CALL TRANSFORM_TO_FVW_COORDS( V )
           Vortexpoints( :, indx2, indx1 ) = VortexpointsJmin1( :, indx2, indx1-1 ) + delta_psi(1) / Omega * V
        END DO ! NumBS+1
     END DO ! jmax+1
  END IF

  DO indx1 = Num_start, NumBS
     B( indx1 )  = DOT_PRODUCT( normalvector( :, indx1 ), VTotal( :, indx1 ) )
     DO indx2 = Num_start, NumBS
        Call BiotSavart( Vortexpoints( :, indx2, 1 ), Vortexpoints( :, indx2+1, 1 ), &
           & Controlpoints( :, indx1 ), BS( :, 1 ))
        DO m = 1, jmax
           Call BiotSavart( Vortexpoints( :, indx2, m ), Vortexpoints( :, indx2, m+1 ), &
           & Controlpoints( :, indx1 ), BS( :, m*2 ))
           Call BiotSavart( Vortexpoints( :, indx2+1, m ), Vortexpoints( :, indx2+1, m+1 ), &
           & Controlpoints( :, indx1 ), BS( :, m*2+1 ))
        END DO

        sumBS = 0.0_ReKi
        DO m = 2, (( 2 * NnearMax ) - 1 )
           sumBS = sumBS + ( -1.0_ReKi ) ** dble( m + 1 ) * BS( :, m )
        END DO
        A2( :, indx2, indx1 ) = sumBS

        m=1
        sumBS = sumBS + ( -1.0_ReKi ) **dble( m + 1 ) * BS( :, m)
        A( indx2, indx1 ) = DOT_PRODUCT( sumBS, normalvector( :, indx1 ) )
     END DO ! NumBS
  END DO ! NumBS

  PRINT*, 'Before PINV'

  CALL pinv( A, NumBS, Ainv )
  PRINT*, 'After PINV'

  Gammabl = matmul( Ainv, B )

  VinducedNWFinal( 1, : ) = -matmul( A2( 1, :, : ), Gammabl )
  VinducedNWFinal( 2, : ) = -matmul( A2( 2, :, : ), Gammabl )
  VinducedNWFinal( 3, : ) = -matmul( A2( 3, :, : ), Gammabl )

  Cap_Gamma_max = maxval( Gammabl( Num_start:NumBS ))
  Cap_Gamma_min = minval( Gammabl( Num_start:NumBS ))
  IF ( abs( Cap_Gamma_min ) .GT. abs( Cap_Gamma_max )) THEN
     Cap_Gamma = Cap_Gamma_min
  ELSE
     Cap_Gamma = Cap_Gamma_max
  END IF
  Gamma_near( Num_start ) = -Gammabl( Num_start )
  Gamma_near( NumBS + 1 ) = Gammabl( NumBS )
  DO indx1 = Num_start + 1, NumBS
     Gamma_near( indx1 ) = Gammabl( indx1 - 1 ) - Gammabl( indx1 )
  END DO

  Ainv = 0.0_ReKi
  VinducedNWFinal = 0.0_ReKi      !KS -- Why is this set to 0??

END SUBROUTINE Calculate_Gamma1
! =====================================================================================

! =====================================================================================
   SUBROUTINE CLCD_FVW( P,  O,  ErrStat, ErrMess, &
                    ALPHA, CLA, CDA, CMA, I )
!   SUBROUTINE CLCD( ALPHA, CLA, CDA, CMA, I, ErrStat )
 !  returns values of lift and drag coeffs.
 !   This subroutine interpolates airfoil coefficients
 !   from a table of airfoil data.  The table must consist
 !   of ALPHA, CL and CD over the entire range of angles
 !   that will be encountered.
 !
 ! VARIABLES:
 !    CLA      = Returned value of lift coefficient
 !    CDA      = Returned value of drag coeff
 !    CMA      = Returned value of pitching moment coeff
 !    ALPHA    = Angle of attack ( radians )
 !    AL       = Array containing the angle of attack
 !    CL       = Array containing the lift coeffs. at AL( I )
 !    CD       = Array containing the drag coeffs. at AL( I )
 !    CM       = Array containing the moment coeffs. at AL( I )
 !    I        = Airfoil ID for this element, equal to NFoil( J ), where J is the index identifying the blade element
 !    MulTabLoc = Multiple airfoil table location for this element
 !    MulTabMet = Array containing the multiple airfoil table metric

 !  !!!!!KS!!!!!!
 !  Took from AeroSubs.f90 and altered/put here for ease of data transfer.
 ! ******************************************************

   USE AeroDyn14_Types

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE( AirfoilParms ), INTENT( IN    ) :: p ! Parameters
   INTEGER,              INTENT( IN    ) :: I ! NFOIL( J )

   TYPE( Airfoil ),      INTENT( INOUT ) :: O !therState ! Initial other / optimization states
   REAL( ReKi ),         INTENT( INOUT ) :: ALPHA

   INTEGER,              INTENT(   OUT ) :: ErrStat
   CHARACTER( * ),       INTENT(   OUT ) :: ErrMess
   REAL( ReKi ),         INTENT(   OUT ) :: CDA, CLA, CMA

   ! Local Variables:
   INTEGER      :: N1, N1P1, N2, N2P1, NTAB
   REAL( ReKi ) :: CDA1, CDA2, CLA1, CLA2, CMA1, CMA2, P1, P2

  ErrStat = ErrID_None
  ErrMess = ""

  IF ( .NOT. ALLOCATED( P%NFoil )) THEN
     CDA = 0.; CLA = 0.; CMA = 0.
     ErrStat = ErrID_Fatal
     RETURN
  ELSE
     ErrStat = ErrID_None
  END IF

  NTAB = P%NLIFT( I )

  IF (( ALPHA < O%AL( I, 1 )) .OR. ( ALPHA > O%AL( I, NTAB )) )   THEN
!bjj: This error message isn't necessarially accurate:
     ErrMess = ' Angle of attack = ' // TRIM( Num2LStr( ALPHA * R2D )) // &
               ' deg is outside data table range. ' // & !Blade #' // TRIM( Int2LStr( IBLADE )) // &
               ' Airfoil ' // TRIM( Int2LStr( I )) // '.'
!                   ' element ' // TRIM( Int2LStr( J )) // '.' )

     ErrStat = ErrID_Fatal
     RETURN
  ENDIF

  ALPHA = MIN( MAX( ALPHA, O%AL( I, 1 )), O%AL( I, NTAB ))
  CALL LocateBin ( ALPHA, O%AL( I, 1:NTAB ), N1, NTAB )

  IF ( N1 == 0 ) THEN
     N1   = 1
     N1P1 = 2
     P1   = 0.0
  ELSEIF( N1 == NTAB ) THEN
     N1P1 = N1
     N1   = N1 - 1
     P1   = 1.0
  ELSE
     N1P1 = N1 + 1
     P1   = ( ALPHA - O%AL( I, N1 )) / ( O%AL( I, N1P1 ) - O%AL( I, N1 ))
  END IF

 ! If the element has multiple airfoil tables, do a 2-D linear interpolation
 !  for Cl and CD

  IF ( P%NTables( I ) > 1 ) THEN

     O%MulTabLoc = MIN( MAX( O%MulTabLoc, P%MulTabMet( I, 1 )), P%MulTabMet( I, P%NTables( I )) )
     CALL LocateBin ( O%MulTabLoc, P%MulTabMet( I, 1:P%NTables( I )), N2, P%NTables( I ))

     IF ( N2 == 0 ) THEN
        N2   = 1
        N2P1 = 2
        P2   = 0.0
     ELSE IF ( N2 == P%NTables( I )) THEN
        N2P1 = N2
        N2   = N2 - 1
        P2   = 1.0
     ELSE
        N2P1 = N2 + 1
        P2   = ( O%MulTabLoc - P%MulTabMet( I, N2 )) / ( P%MulTabMet( I, N2P1 ) - P%MulTabMet( I, N2 ))
     END IF

     CLA1 = O%CL( I, N1, N2 ) + P1 * ( O%CL( I, N1P1, N2 ) - O%CL( I, N1, N2 ))
     CDA1 = O%CD( I, N1, N2 ) + P1 * ( O%CD( I, N1P1, N2 ) - O%CD( I, N1, N2 ))
     CMA1 = O%CM( I, N1, N2 ) + P1 * ( O%CM( I, N1P1, N2 ) - O%CM( I, N1, N2 ))

     CLA2 = O%CL( I, N1, N2P1 ) + P1 * ( O%CL( I, N1P1, N2P1 ) - O%CL( I, N1, N2P1 ))
     CDA2 = O%CD( I, N1, N2P1 ) + P1 * ( O%CD( I, N1P1, N2P1 ) - O%CD( I, N1, N2P1 ))
     CMA2 = O%CM( I, N1, N2P1 ) + P1 * ( O%CM( I, N1P1, N2P1 ) - O%CM( I, N1, N2P1 ))

     CLA = CLA1 + P2 * ( CLA2 - CLA1 )
     CDA = CDA1 + P2 * ( CDA2 - CDA1 )
     CMA = CMA1 + P2 * ( CMA2 - CMA1 )

     ELSE

     CLA  = O%CL( I, N1, 1 ) + P1 * ( O%CL( I, N1P1, 1 ) - O%CL( I, N1, 1 ))
     CDA  = O%CD( I, N1, 1 ) + P1 * ( O%CD( I, N1P1, 1 ) - O%CD( I, N1, 1 ))
     CMA  = O%CM( I, N1, 1 ) + P1 * ( O%CM( I, N1P1, 1 ) - O%CM( I, N1, 1 ))

  ENDIF

RETURN
END SUBROUTINE CLCD_FVW
! =====================================================================================
