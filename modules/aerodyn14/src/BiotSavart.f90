Subroutine BiotSavart ( rp1, rp2, rp3, BS )

  USE NWTC_Library
  USE MathOps,  Only: cross, norm, dot

  IMPLICIT NONE

  REAL( ReKi ), DIMENSION( 3 ) :: rp1, rp2, rp3, BS
  REAL( ReKi ), DIMENSION( : ), ALLOCATABLE :: r1, crossr1r2, r2
  REAL( ReKi ) :: dotr1r2, normr1, normr2

  ALLOCATE( r1( 3 ), crossr1r2( 3 ), r2( 3 ))

  r1 = rp3 - rp1
  r2 = rp3 - rp2
  BS = 0.00_ReKi
  CALL cross( r1, r2, crossr1r2 )
  CALL norm( r1, normr1 )
  CALL norm( r2, normr2 )
  CALL dot( r1, r2, dotr1r2 )

  IF ( abs( normr1 ) .GT. 0.00_ReKi .AND. abs( normr2 ) .GT. 0.00_ReKi .AND. &
    & abs( normr1 * normr2 + dotr1r2 ) .GT. 0.00_ReKi ) THEN
     BS = 1.00_ReKi / ( 4.00_ReKi * Pi_D ) * crossr1r2 * ( 1.00_ReKi / normr1 + 1.00_ReKi / normr2 ) * &
       & ( 1.00_ReKi / ( normr1 * normr2 + dotr1r2 ))
  END IF

  DEALLOCATE( r1, crossr1r2, r2 )

END Subroutine BiotSavart
