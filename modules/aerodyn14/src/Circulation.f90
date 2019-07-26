! =====================================================================================
SUBROUTINE Calculate_Gamma1( n, VTotal, BladeTanVect, normalvector, BladeLoc, ControlPoints, Cap_Gamma, &
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

