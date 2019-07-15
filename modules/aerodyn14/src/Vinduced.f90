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

  REAL :: WhichTurb

  INTEGRAL = 0.00_ReKi; rc=0.00_ReKi

  WhichTurb = REAL(n-0.01)/REAL(NumBl)
  q = n - FLOOR( WhichTurb )*NumBl

  DO i = 1, NumWakes
     WhichTurb = REAL(i-0.01)/REAL(NumBl)

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

  REAL :: WhichTurb


  INTEGRAL = 0.00_ReKi; rc=0.00_ReKi

  WhichTurb = REAL(n-0.01)/REAL(NumBl)
  q = n - FLOOR( WhichTurb)*NumBl

  DO i = 1, NumWakes
     WhichTurb = REAL(i-0.01)/REAL(NumBl)

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

  REAL :: WhichTurb

  INTEGRAL = 0.0_ReKi; rc=0.0_ReKi
  DO i = 1, NumWakes
     WhichTurb = REAL(i-0.01)/REAL(NumBl)
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
