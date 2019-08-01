MODULE FVW_Subs


   USE NWTC_Library
   USE FVW_Types
   USE InflowWind

   IMPLICIT NONE



CONTAINS

! =====================================================================================
SUBROUTINE Calculate_Gamma1( p, n, VTotal, BladeTanVect, normalvector, BladeLoc, ControlPoints, Cap_Gamma, &
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


  IMPLICIT NONE

  type(FVW_ParameterType),                           intent( in    ) :: p              !< Parameters
  INTEGER,                                           INTENT( IN    ) :: n
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
  REAL( ReKi ), DIMENSION( NumBS, NumBS )                            :: Ainv

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

!FIXME: Err handling is not complete!
  INTEGER(IntKi)                      :: ErrStat2
  CHARACTER(ErrMsgLen)                :: ErrMsg2
 
 
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
           CALL TRANSFORM_TO_AERODYN_COORDS( p, TMP_Vect )

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

    ! Get inverse of A
  CALL Pinv( A, NumBS, Ainv, ErrStat2, ErrMsg2 )

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

  VinducedNWFinal = 0.0_ReKi      !KS -- Why is this set to 0??

CONTAINS
   !=================================================
   !> Calculate the inverse of the square matrix A using single value decomposition
   !! routines in the LAPACK library
   SUBROUTINE Pinv(A, NumBS, Ainv, ErrStat, ErrMsg)

      USE NWTC_Library
      USE NWTC_LAPACK

      IMPLICIT NONE

      INTEGER,                intent(in   )  :: NumBS
      REAL(ReKi),             intent(inout)  :: A(NumBS,NumBS)
      REAL(ReKi),             intent(  out)  :: Ainv(NumBS,NumBS)
      INTEGER(IntKi),         intent(  out)  :: ErrStat
      CHARACTER(ErrMsgLen),   intent(  out)  :: ErrMsg

      INTEGER(IntKi)                      :: lwork, lwmax
      INTEGER(IntKi)                      :: summation, i

      REAL( ReKi )                        :: tolerance
      REAL( ReKi ),  ALLOCATABLE          :: WORK(:)
      REAL( ReKi )                        :: S(NumBS), U(NumBS,NumBS), VT(NumBS,NumBS), S_mat(NumBS,NumBS)

      INTEGER(IntKi)                      :: ErrStat2
      CHARACTER(ErrMsgLen)                :: ErrMsg2
      CHARACTER(*),  PARAMETER            :: RoutineName='Pinv'

      ErrStat = ErrID_None
      ErrMsg  = ""


   !FIXME: To optimize, setup the work array as a miscvar, calculate optimal size at init and set
      !--------------------------
      !  Find size of work array.
      !--------------------------

         ! set the size of the work array to something that we know will be possible to use
      LWMAX = MIN(7*NumBS,1000)
      ALLOCATE(WORK(LWMAX))
      work = 0.0_ReKi

         ! Query the [d,s]gesvd LAPACK routines to find out the optimal size for the work array.
      LWORK = -1
      call LAPACK_gesvd('A', 'A', NumBS, NumBS, A, S, U, Vt, work, lwork, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! If LAPACK (or MKL) suggested a larger work array as optimal, reallocate work array.
      if (LWMAX < work(1) ) then
         LWMAX=work(1)
         deallocate(work)
         allocate(work(LWMAX))
         work = 0.0_ReKi
      endif

      !--------------------------
      !  Compute SVD.
      !--------------------------

      call LAPACK_gesvd('A', 'A', NumBS, NumBS, A, S, U, Vt, work, size(work), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! To speed up the calculations, find the tolerance and only calculate results
         ! for terms above the tolerance
      tolerance = NumBS*epsilon(maxval(S))

      summation=0
      DO i=1,NumBS
         IF (s(i) .GT. tolerance) THEN
            summation=summation+1;
         END IF
      END DO

         ! Set the diagonal elements of S_mat
      S_mat = 0.0_ReKi
      DO i = 1, NumBS
         IF (i .LE. summation)THEN
            S_mat(i,i)=1.0_ReKi/s(i)
         END IF
      END DO

         ! Calculate the inverse of A
      Ainv=transpose(matmul( matmul(U(:,1:summation),S_mat(1:summation,1:summation)), VT(1:summation,:)))

      DEALLOCATE(WORK)

   END SUBROUTINE Pinv
   !=================================================




END SUBROUTINE Calculate_Gamma1
! =====================================================================================






Subroutine BiotSavart ( rp1, rp2, rp3, BS )

  USE NWTC_Library

  IMPLICIT NONE

  REAL( ReKi ), DIMENSION( 3 ) :: rp1, rp2, rp3, BS
  REAL( ReKi ), DIMENSION( : ), ALLOCATABLE :: r1, crossr1r2, r2
  REAL( ReKi ) :: dotr1r2, normr1, normr2

  ALLOCATE( r1( 3 ), crossr1r2( 3 ), r2( 3 ))

  r1 = rp3 - rp1
  r2 = rp3 - rp2
  BS = 0.00_ReKi
  crossr1r2  = CROSS_PRODUCT( r1, r2 )
  normr1 = TwoNorm( r1 )
  normr2 = TwoNorm( r2 )
  dotr1r2  = DOT_PRODUCT( r1, r2 )

  IF ( abs( normr1 ) .GT. 0.00_ReKi .AND. abs( normr2 ) .GT. 0.00_ReKi .AND. &
    & abs( normr1 * normr2 + dotr1r2 ) .GT. 0.00_ReKi ) THEN
     BS = 1.00_ReKi / ( 4.00_ReKi * Pi_D ) * crossr1r2 * ( 1.00_ReKi / normr1 + 1.00_ReKi / normr2 ) * &
       & ( 1.00_ReKi / ( normr1 * normr2 + dotr1r2 ))
  END IF

  DEALLOCATE( r1, crossr1r2, r2 )

END Subroutine BiotSavart

!FIXME: make this a function, and remove the parameter input (can pass p%TransformFVWtoAD instead?)
SUBROUTINE TRANSFORM_TO_AERODYN_COORDS(p,Vector)!,zloc)

   USE FVW_Parm, Only: HH

   IMPLICIT NONE

   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   real(ReKi),                      intent(  out)  :: Vector(3)
 
   REAL(ReKi)                                      :: TmpVector(3)
   REAL(ReKi)                                      :: Transf(3,3)
 
   TmpVector=MatMul(p%TransformFVWtoAD,Vector)
 
   Vector=TmpVector
 
   Vector(3)=Vector(3)+HH
   Vector(1)=0.0_ReKi!  KS--I don't even understand why that was there.... 3.10.15

END SUBROUTINE TRANSFORM_TO_AERODYN_COORDS
! ==================================================================

! ==================================================================
SUBROUTINE TRANSFORM_TO_FVW_COORDS(Vector)


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
SUBROUTINE Vinduced2OLD( rblade, Gamma, rp, rblade2, n, jold, kold )

  USE NWTC_Library
  USE FVW_Parm
  USE MultTurb_Params,    Only: NumWakes, NTurb, FWake!GCoord

  IMPLICIT NONE

  INTEGER,                                           INTENT( IN    ) :: n, jold, kold
  REAL( ReKi ), DIMENSION( 3 ),                            INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: rblade, rblade2


  INTEGER :: i, kx, indx, q, limit, BC_loc

  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2
  REAL( ReKi ) :: delta, strain, len1, len2, zeta, rc0, Denom
  REAL( ReKi ) :: r1(3), r2(3), crossr1r2(3), rbladetemp(3)
  REAL( ReKi ) :: INTEGRAL, rc

  REAL( ReKi ) :: WhichTurb

  INTEGRAL = 0.00_ReKi; rc=0.00_ReKi

  WhichTurb = (REAL(n,ReKi) -0.01)/REAL(NumBl, ReKi)
  q = n - FLOOR( WhichTurb )*NumBl

  DO i = 1, NumWakes
     WhichTurb = (REAL(i,ReKi)-0.01)/REAL(NumBl, ReKi)

     limit = CUTOFF_up( CEILING( WhichTurb ))
     BC_loc = BC( CEILING( WhichTurb ) )
     DO kx = NnearMax, WakeAgeLimit-1
        IF ( kx .LT. limit ) THEN
           zeta = dble( kx ) * delta_psi(1)

           DO indx = 1, 3
              r2( indx ) = rp( indx ) - rblade( indx, kx,   i )
              r1( indx ) = rp( indx ) - rblade( indx, kx+1, i )
           END DO

           mag_r1 = TwoNorm( r1 )
           mag_r2 = TwoNorm( r2 )
           dotr1r2  = DOT_PRODUCT( r1, r2 )

           delta = 1.00_ReKi + a1 * ( abs( Gamma( kx, i ))) / nu
           rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade(  :, kx+1, i ) - rblade(  :, kx, i )
           len2 = TwoNorm( rbladetemp )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade2( :, kx,   i ) - rblade2( :, kx-1, i )
           len1 = TwoNorm( rbladetemp )

           strain = ( len2 - len1 ) / len1

           INTEGRAL = delta_psi(1) / ( 1.00_ReKi + strain )
           rc = sqrt(( rc0 * rc0 + 4.00_ReKi * alpha_param * delta * nu * zeta / &
              & Omega ) * INTEGRAL )

           IF ( kx .GE. 2 ) THEN
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           crossr1r2  = CROSS_PRODUCT( r1, r2 )

           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
              & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.00_ReKi * &
              & dotr1r2 ) **2.00_ReKi

           IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
              & (denom .NE. 0.00_ReKi)) THEN
              DO indx = 1, 3
                 FWake%VinducedFarWakej( indx, kold, q, kx+1, i ) = Gamma( kx, i ) / &
                    & ( 4.00_ReKi * Pi_D ) * crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * &
                    & ( 1.00_ReKi - dotr1r2 / ( mag_r1 * mag_r2 )) / sqrt(denom)
              END DO
           ELSE
              FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
           END IF
        ELSE IF ( jold .GE. ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc ) .AND. &
               & kx .LE. ( jold - ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc - limit )) .AND. &
               & kx .LT. ( limit  + NINT( TwoPi_D / delta_psi_Est ))) THEN
           CALL VinducedFWOLD( i, kx+1, q, kold )
        ELSE IF (kx .LE. WakeAgeLimit ) THEN
           FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
        END IF
     END DO
  END DO


END SUBROUTINE Vinduced2OLD


!-----------------------------------------------------------------------

SUBROUTINE Vinduced2PRIME( rblade, Gamma, rp, rblade2, n, jold, kold )

  USE NWTC_Library
  USE FVW_Parm
  USE MultTurb_Params,    Only: NumWakes, NTurb, FWake!GCoord

  IMPLICIT NONE

  INTEGER,                                           INTENT( IN    ) :: n, jold, kold
  REAL( ReKi ), DIMENSION( 3 ),                            INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: rblade, rblade2


  INTEGER :: i, kx, indx, q, limit, BC_loc

  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2
  REAL( ReKi ) :: delta, strain, len1, len2, zeta, rc0, denom, INTEGRAL, rc
  REAL( ReKi ) :: rbladetemp(3), r1(3), r2(3), crossr1r2(3)

  REAL( ReKi ) :: WhichTurb


  INTEGRAL = 0.00_ReKi; rc=0.00_ReKi

  WhichTurb = (REAL(n,ReKi)-0.01)/REAL(NumBl, ReKi)
  q = n - FLOOR( WhichTurb)*NumBl

  DO i = 1, NumWakes
     WhichTurb = (REAL(i,ReKi)-0.01)/REAL(NumBl, ReKi)

        limit = CUTOFF_up( CEILING( WhichTurb ))
        BC_loc = BC( CEILING( WhichTurb ) )

     DO kx = NnearMax, WakeAgeLimit-1
        IF ( kx .LT. limit ) THEN
           zeta = dble( kx ) * delta_psi(1)
           DO indx = 1, 3
              r2( indx ) = rp( indx ) - rblade( indx, kx,   i )
              r1( indx ) = rp( indx ) - rblade( indx, kx+1, i )
           END DO

           mag_r1 = TwoNorm( r1 )
           mag_r2 = TwoNorm( r2 )
           dotr1r2  = DOT_PRODUCT( r1, r2 )

           delta = 1.00_ReKi + a1 * ( abs( Gamma( kx, i ))) / nu
           rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade(  :, kx+1, i ) - rblade(  :, kx, i )
           len2 = TwoNorm( rbladetemp )
           rbladetemp = 0.00_ReKi
           rbladetemp = rblade2( :, kx,   i ) - rblade2( :, kx-1, i )
           len1 = TwoNorm( rbladetemp )

           strain = ( len2 - len1 ) / len1

           INTEGRAL = delta_psi(1) / ( 1.00_ReKi + strain )
           rc = sqrt(( rc0 * rc0 + 4.00_ReKi * alpha_param * delta * nu * zeta / &
              & Omega ) * INTEGRAL )

           IF ( kx .GE. 2 ) THEN
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           crossr1r2  = CROSS_PRODUCT( r1, r2 )

           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
              & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.00_ReKi * &
              & dotr1r2 ) **2.00_ReKi

           IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
              & ( denom .NE. 0.00_ReKi )) THEN

              DO indx = 1, 3
                 FWake%VinducedFarWakej( indx, kold, q, kx+1, i ) = Gamma( kx, i ) / &
                    & ( 4.00_ReKi * Pi_D ) * crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * &
                    & ( 1.00_ReKi - dotr1r2 / ( mag_r1 * mag_r2 )) / sqrt(denom)
              END DO
           ELSE
              FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
           END IF
        ELSE IF ( jold .GE. ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc ) .AND. &
          & kx .LE. ( jold - ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc - limit )) .AND. &
          & kx .LT. ( limit  + NINT( TwoPi_D / delta_psi_Est ))) THEN
           CALL VinducedFWPRIME( i, kx+1, q, kold )
        ELSE IF (kx .LE. WakeAgeLimit ) THEN
           FWake%VinducedFarWakej( :, kold, q, kx+1, i ) = 0.00_ReKi
        END IF
     END DO
  END DO

END SUBROUTINE Vinduced2PRIME


!-----------------------------------------------------------------------

SUBROUTINE Vinduced3( rblade, Gamma, rp, rblade2, n, jold, kold )

 !******************************************************************
 !* Computes induced velocity on NEAR WAKE due to the Far Wake
 !*
 !*
 !******************************************************************

  USE FVW_Parm
  USE NWTC_Library
  USE MultTurb_Params,  Only: NumWakes, NTurb, FWake!GCoord

  IMPLICIT NONE

  INTEGER, INTENT( IN ) :: n, jold, kold   !=Blade #
  REAL( ReKi ), DIMENSION( 3                            ), INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, CUTOFF_Allocate, NumWakes ), INTENT( IN    ) :: rblade, rblade2



  INTEGER i, kx, indx, limit, BC_loc, a
  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2, delta, strain
  REAL( ReKi ) :: len1, len2, zeta, rc0, denom, rc, INTEGRAL
  REAL( ReKi ) :: r1(3), r2(3), crossr1r2(3)
  REAL( ReKi ) :: rbladetemp(3)

  REAL( ReKi ) :: WhichTurb

  INTEGRAL = 0.0_ReKi; rc=0.0_ReKi
  DO i = 1, NumWakes
     WhichTurb = (REAL(i, ReKi) -0.01)/REAL(NumBl, ReKi)
        limit = CUTOFF_up( CEILING( WhichTurb ))
        BC_loc = BC( CEILING( WhichTurb ) )

     DO kx = NnearMax, WakeAgeLimit-1
        IF ( kx .LT. limit ) THEN
           r2( : ) = rp( : ) - rblade( :, kx,   i )
           r1( : ) = rp( : ) - rblade( :, kx+1, i )

           mag_r1 = TwoNorm( r1 )
           mag_r2 = TwoNorm( r2 )
           dotr1r2  = DOT_PRODUCT( r1, r2 )

           delta = 1.0_ReKi + a1 * ( abs( Gamma( kx, i ))) / nu
           rc0 = sqrt( 4.0_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )
           rbladetemp = 0.0_ReKi
           rbladetemp = rblade(  :, kx+1, i ) - rblade(  :, kx, i )
           len2 = TwoNorm( rbladetemp )
           rbladetemp = 0.0_ReKi
           rbladetemp = rblade2( :, kx,   i ) - rblade2( :, kx-1, i )
           len1 = TwoNorm( rbladetemp )

           zeta = dble( kx ) * delta_psi(1)

           strain = ( len2 - len1 ) / len1
           INTEGRAL = delta_psi(1) / ( 1.0_ReKi + strain )
           rc = sqrt(( rc0 * rc0 + 4.0_ReKi * alpha_param * delta * nu * zeta / &
              & Omega ) * INTEGRAL )

           IF ( kx .GE. 2 ) THEN
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           crossr1r2  = CROSS_PRODUCT( r1, r2 )

           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.0_ReKi + &
              & close_to_zero **4.0_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.0_ReKi * &
              & dotr1r2 ) **2.0_ReKi

           IF ( mag_r1 .GT. 0.0_ReKi .AND. mag_r2 .GT. 0.0_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.0_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.0_ReKi .OR. crossr1r2( 3 ) .NE. 0.0_ReKi ) .AND. &
              & ( denom .NE. 0.0_ReKi )) THEN

              DO indx=1, 3
                 FWake%VinducedFarWakeRj( indx, kold, n, kx+1, i ) = Gamma( kx, i ) / ( 4.0_ReKi * Pi_D ) * &
                    & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.0_ReKi - dotr1r2 / &
                    & ( mag_r1 * mag_r2 )) / sqrt( denom )
              END DO
           ELSE
              FWake%VinducedFarWakeRj( :, kold, n, kx+1, i ) = 0.0_ReKi
           END IF
        ELSE IF ( jold .GE. ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc ) .AND. &
           & kx .LE. ( jold - ( NINT( TwoPi_D / delta_psi_Est ) * BC_loc - limit)) .AND. &
           & kx .LT. ( limit + NINT( TwoPi_D / delta_psi_Est ))) THEN
           CALL VinducedFW3( i, kx+1, n, kold )
        ELSE
           FWake%VinducedFarWakeRj( :, kold, n, kx+1, i ) = 0.0_ReKi
        END IF
     END DO
  END DO

END SUBROUTINE Vinduced3


!-----------------------------------------------------------------------

SUBROUTINE VinducedBC( rblade, Gamma, rp, Vind )


  USE NWTC_Library
  USE FVW_Parm
  USE MultTurb_params,  Only: NumWakes

  IMPLICIT NONE

  REAL( ReKi ), DIMENSION( 3 ), INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION(    NumBS+1, NumWakes ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NumWakes ), INTENT( IN    ) :: rblade

  INTEGER  i, indx, nbs

  REAL( ReKi ) :: close_to_zero, mag_r1, mag_r2, dotr1r2
  REAL( ReKi ) :: delta, rc0, denom
  REAL( ReKi ), DIMENSION( 3 ) :: Vind, r1, r2, crossr1r2


  Vind = 0.00_ReKi

  DO i = 1, NumWakes
     DO nbs = Num_start, NumBS
        DO indx = 1, 3
           r1( indx ) = rp( indx ) - rblade( indx, nbs,   i )
           r2( indx ) = rp( indx ) - rblade( indx, nbs+1, i )
        END DO

        mag_r1 = TwoNorm( r1 )
        mag_r2 = TwoNorm( r2 )
        dotr1r2  = DOT_PRODUCT( r1, r2 )

        delta = 1.00_ReKi + a1 * ( abs( Gamma( nbs, i ))) / nu
        rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )

        close_to_zero = rc0

        crossr1r2  = CROSS_PRODUCT( r1, r2 )

        denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
           & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - &
           & 2.00_ReKi * dotr1r2 ) **2.00_ReKi

        IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR.&
           & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
           & (denom .NE. 0.00_ReKi )) THEN

           DO indx = 1, 3
              Vind( indx ) = Vind( indx ) + Gamma( nbs, i ) / ( 4.00_ReKi * Pi_D ) * &
                 & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.00_ReKi - dotr1r2 / &
                 & ( mag_r1 * mag_r2 )) / sqrt( denom )
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE VinducedBC


!-----------------------------------------------------------------------

SUBROUTINE VinducedFW3(n1, k1, n2, k2)

  USE NWTC_Library
  USE FVW_Parm,       Only: delta_psi_Est
  USE MultTurb_params, Only: FWake!GCoord

  IMPLICIT NONE

  INTEGER, INTENT( IN ) :: n1, n2, k1, k2
  FWake%VinducedFarWakeRj( :, k2, n2, k1, n1 ) = FWake%VinducedFarWakeRjm1( :, k2, n2, k1-1, n1 ) + &
     & FWake%VinducedFarWakeRj( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est ), n1 ) - &
     & FWake%VinducedFarWakeRjm1( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est )-1, n1 )

END SUBROUTINE VinducedFW3


!-----------------------------------------------------------------------

SUBROUTINE VinducedFWOLD( n1, k1, n2, k2 )

  USE NWTC_Library
  USE FVW_Parm, Only: delta_psi_Est
  USE MultTurb_params, Only: FWake!GCoord

  IMPLICIT NONE

  INTEGER n1, n2, k1, k2

  FWake%VinducedFarWakej( :, k2, n2, k1, n1 ) = FWake%VinducedFarWakejm1( :, k2, n2, k1-1, n1 ) + &
     & FWake%VinducedFarWakej( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est ), n1 )- &
     & FWake%VinducedFarWakejm1( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est )-1, n1 )


END SUBROUTINE VinducedFWOLD


!-----------------------------------------------------------------------

SUBROUTINE VinducedFWPRIME( n1, k1, n2, k2 )

  USE NWTC_Library
  USE FVW_Parm, Only: delta_psi_Est
  USE MultTurb_params, Only: FWake!GCoord

  IMPLICIT NONE

  INTEGER n1, n2, k1, k2

  FWake%VinducedFarWakej( :, k2, n2, k1, n1 ) = FWake%VinducedFarWakejm1( :, k2, n2, k1-1, n1 ) +&
     & FWake%VinducedFarWakej( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est ), n1 ) - &
     & FWake%VinducedFarWakej( :, k2, n2, k1-NINT( TwoPi_D/delta_psi_Est )-1, n1 )


END SUBROUTINE VinducedFWPRIME


!-----------------------------------------------------------------------

SUBROUTINE VinducedNW( rblade, Gamma, rp, Vind, rblade2, up )


 !******************************************************************
 !* Computes induced velocity on referance LAGRANGIAN MARKER # "m" due to all of the near wake points
 !*     --> include near wake away from blade up to NnearMax-1 and along blade from Num_Start to NumBS+1
 !*  *** When called in UpdateAeroVals, this instead computes the induced velocity ON THE BLADE
 !*                                     due to the near wake
 !*
 !* dV = F_c*(GAMMA_mu/4pi)*(r_1xr_2)*[1/mag(r_1)+1/mag(r_2)]*{1/[mag(r_1)*mag(r_2)+r_1DOTr_2]}
 !* F_c = mag(r)^2/sqrt[mag(r)^4+r_c^4]
 !* r_c( zeta, epsilon ) = sqrt{r_c0^2+(4*alpha*delta*mu*zeta/Omega)INTEGRAL[(1+epsilon)^-1*dzeta]}
 !*                        ^*                 diffusion          *^ ^*  effect of strain field  *^
 !* r_c0^2 = 4*alpha*delta*mu*zeta_0/Omega; zeta_0 = 30 degrees (location of near-wake roll-up)
 !*                                            ^^^^^hardcoded.
 !*
 !* where:
 !*      dV = Vind          <--- induced velocity in NEAR WAKE (Outputted to VinducedNWx)
 !*      r_blade1, r_blade2 <---distances between the point of interest and the end points of the vortex segment)
 !*                 and these points come from r_nearj, r_nearjm1, or r_nearjm2
 !*      GAMMA = Gamma; comes from Gamma_nearj or Gamma_nearjm1
 !*      INTEGRAL[1+epsilon)^-1*dzeta] = INTEGRAL
 !*            = rp <--- reference point
 !*                  from PREDICTOR: r_oldj(:,m,n), r_oldj(:,m-1,n), r_oldj(:,m,n), r_oldjm1(:,m-1,n)
 !*                       CORRECTOR: r_primej(:,m,n), r_primej(:,m-1,n), r_primejm1(:,m,n), r_primejm1(:,m-1,n)
 !*                       OTHER:     BladeThreeQuarterChord(:, nbs, n),
 !*
 !*
 !*
 !*
 !******************************************************************


  USE NWTC_Library
  USE MultTurb_Params, Only: NumWakes, NTurb
  USE FVW_Parm

  IMPLICIT NONE

  INTEGER,                                             INTENT( IN    ) :: up
  REAL( ReKi ), DIMENSION( 3 ),                              INTENT( IN    ) :: rp
  REAL( ReKi ), DIMENSION( NnearMax, NumBS+1, NumWakes    ), INTENT( IN    ) :: Gamma
  REAL( ReKi ), DIMENSION( 3, NumBS+1, NnearMax, NumWakes ), INTENT( IN    ) :: rblade, rblade2

  REAL( ReKi ), DIMENSION( 3 ),                              INTENT(   OUT ) :: Vind

  INTEGER :: i, kx, indx, nbs

  REAL( ReKi ) :: delta, strain, len1, len2, zeta, rc0, denom, close_to_zero, mag_r1, mag_r2, dotr1r2, rc, INTEGRAL

  REAL( ReKi ) :: r1(3), r2(3), crossr1r2(3)
  REAL( ReKi ) :: rbladetemp(3)

  INTEGRAL = 0.00_ReKi; rc = 0.00_ReKi; Vind = 0.00_ReKi

  DO i = 1, up      !   NumWakes for pred/corr and NumBl for UpdateAero & Circulation
     DO kx = 1, NnearMax - 1
        zeta = dble( kx ) * delta_psi(1)
        DO nbs = Num_start, NumBS+1
              r2( : ) = rp( : ) - rblade( :, nbs, kx,   i )
              r1( : ) = rp( : ) - rblade( :, nbs, kx+1, i )

           mag_r1 = TwoNorm( r1 )
           mag_r2 = TwoNorm( r2 )
           dotr1r2  = DOT_PRODUCT( r1, r2 )

           !Calculate the core radius for the Vind cut off distance
           delta = 1.00_ReKi + a1 * ( abs( Gamma( kx, nbs, i ))) / nu
           rc0 = sqrt( 4.00_ReKi * alpha_param * delta * nu * 30.0_ReKi * D2R_D / Omega )

           IF ( kx .GE. 2 ) THEN
              rbladetemp = 0.00_ReKi
              rbladetemp = rblade(  :, nbs, kx+1, i ) - rblade(  :,  nbs,  kx, i )
              len2 = TwoNorm( rbladetemp )
              rbladetemp = 0.00_ReKi
              rbladetemp = rblade2( :, nbs, kx,   i ) - rblade2( :, nbs, kx-1, i )
              len1 = TwoNorm( rbladetemp )
              strain = ( len2 - len1 ) / len1
              INTEGRAL = delta_psi(1) / ( 1.0_ReKi + strain )
              rc = sqrt(( rc0 * rc0 + 4.00_ReKi * alpha_param * delta * nu * zeta / &
                 & Omega ) * INTEGRAL )
              close_to_zero = rc
           ELSE
              close_to_zero = rc0
           END IF

           crossr1r2  = CROSS_PRODUCT( r1, r2 )

           denom = ( mag_r1 * mag_r1 * mag_r2 * mag_r2 - dotr1r2 * dotr1r2 ) **2.00_ReKi + &
              & close_to_zero **4.00_ReKi * ( mag_r1 * mag_r1 + mag_r2 * mag_r2 - 2.00_ReKi * dotr1r2 ) **2.00_ReKi

           IF ( mag_r1 .GT. 0.00_ReKi .AND. mag_r2 .GT. 0.00_ReKi .AND. ( crossr1r2( 1 ) .NE. 0.00_ReKi .OR. &
              & crossr1r2( 2 ) .NE. 0.00_ReKi .OR. crossr1r2( 3 ) .NE. 0.00_ReKi ) .AND. &
              & ( denom .NE. 0.00_ReKi )) THEN

              DO indx = 1, 3
                 Vind( indx ) = Vind( indx ) + Gamma( kx, nbs, i ) / ( 4.00_ReKi * Pi_D ) * &
                    & crossr1r2( indx ) * ( mag_r1 + mag_r2 ) * ( 1.00_ReKi - dotr1r2 / &
                    & (mag_r1 * mag_r2 )) / sqrt( denom )
              END DO
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE VinducedNW



END MODULE FVW_Subs
